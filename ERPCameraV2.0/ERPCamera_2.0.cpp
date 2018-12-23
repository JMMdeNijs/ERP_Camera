/*********************************************************************************************************
 
            Main Program Virtual ERP Camera Version 2.0
            ********************************************

This code captures a regular 2D image from a Equirectangular Picture (ERP) using a virtual camera positioned  
in the center of the "ERP sphere". The camera can rotate (pan-azimuth and tilt-inclination) and zoom in and 
out using the keyboard 


The code uses the GPU (OpenCL) to calculate and orientation (phi, theta) of each pixel and 
gather the RGB values from the ERP image. Menory demand is substantial because half the ERP image is 
loaded in the GPU memory.
 
Jan de Nijs Copyright (C), December 2018
 
*****************************************************************************************/

#include "erplib_2.0.h"

 
int main(int argc,char **argv) {
     
    //Get path and file image
	//std::string imDir(GetFile());
	std::string imDir("/home/jan/Projects/ERPCamera/Laon_Cathedral.jpg");
     
    //Load equirectangle image
    CImg<unsigned char> inputImageChannels(&imDir[0]); //In CImg, RGB are stored in channels 
    CImgDisplay vcDisplay;
	
	//Specification of the width and height of the virtual outputImage and ERP image 
	int const outputImageSize[2] = {1024, 576};
	int const inputImageSize[2] = {inputImageChannels.width(), inputImageChannels.height()};
	int const inputSize = inputImageChannels.width() * inputImageChannels.height();
		
	//Segmentation ERP image to reduce (half) memory demand global device memory
	//Only half of the image will be loaded into the memory
	//but this memory is refreshed if the camera image moves outside loaded window
	WindowSegmentation segmentationTheta;
	segmentationTheta.numberRows = inputImageSize[1]/2 + inputImageSize[1]%2;
	
	for (int ii=0; ii< HORIZONTAL_SEGMENT_OPTIONS ; ii++){
		segmentationTheta.index.start[ii] = ((inputImageSize[1] / 2) * ii) / (HORIZONTAL_SEGMENT_OPTIONS-1);
		segmentationTheta.index.end[ii]   = segmentationTheta.index.start[ii] + segmentationTheta.numberRows - 1;	
		segmentationTheta.angle.start[ii] = 90.0f * ( 1.0f - 2.0f*(float)segmentationTheta.index.start[ii] / (float)(inputImageSize[1] - 1));
		segmentationTheta.angle.end[ii]   = 90.0f * ( 1.0f - 2.0f*(float)segmentationTheta.index.end[ii]   / (float)(inputImageSize[1] - 1));
	}
	
	int segmentNo(HORIZONTAL_SEGMENT_OPTIONS/2); //segment theta +45 up to -45 is selected
	int   erpRowRange[2]   = {segmentationTheta.index.start[segmentNo], segmentationTheta.index.end[segmentNo]};
	float erpAngleRange[2] = {segmentationTheta.angle.start[segmentNo], segmentationTheta.angle.end[segmentNo]};
	
	unsigned char *outputImageVector;
	outputImageVector = new unsigned char[outputImageSize[0] * outputImageSize[1] * 4]; //image is stored as RGBA vectors
    CImg<unsigned char> outputImageChannels(outputImageSize[0], outputImageSize[1],1,3,0);
	
	//Timing & Clock  STD::CHRONO CAN CAUSE CORE DUMP!!!!!
	std::chrono::time_point<std::chrono::system_clock> clock_start, clock_end; 
	int clock_interval;
	char *textstring;
	std::string time_interval;
	
	//Colors to show information in display
	unsigned char const yellow[3] = {0, 255, 255}, black[3] = {0, 0, 0};
	
	//Camera parameters: pose and window and window utility parameters
	float theta0(0.), phi0(0.), thetaWindow(M_PI*39.375/180), phiWindow(M_PI*70.0/180);
	float thetaMax, thetaMin, thetaMaxPrevious, thetaMinPrevious;
	float const thetaWindowMax(M_PI*60/180);
	
	//Declaration camera sensor pose in Euclidean coordinates
	SensorPose sensor;
	float	originSensor[3], incWidth[3], incHeight[3];
	
	//Declaration incTheta and incPhi
	float const resolutionInputImage[2] = {2*M_PI/inputImageSize[0], M_PI/inputImageSize[1]};
	
    
	//Arrangements for parallelization code with OpenCL
	std::vector<cl::Platform> platforms;
	cl::Platform::get(&platforms);
	std::vector<cl::Device> devices;
	platforms[0].getDevices(CL_DEVICE_TYPE_GPU, &devices);
	cl::Context context(devices);
	cl::CommandQueue queue = cl::CommandQueue(context, devices[0]);
	//Create global memory for deviceCameraImage with RGB unsigned char channels
	cl::ImageFormat imageFormat = cl::ImageFormat(CL_RGBA, CL_UNSIGNED_INT8);
	cl::Image2D deviceOutputImage = cl::Image2D(context, CL_MEM_WRITE_ONLY, imageFormat, outputImageSize[0], outputImageSize[1]);
	//Create global memory for segment of the the (ERP) inputImageChannels
	cl::ImageFormat inputImageFormat = cl::ImageFormat(CL_R, CL_UNSIGNED_INT8);
	cl::Image2D deviceInputImageR = cl::Image2D(context, CL_MEM_READ_ONLY, inputImageFormat, inputImageSize[0], segmentationTheta.numberRows);
	cl::Image2D deviceInputImageG = cl::Image2D(context, CL_MEM_READ_ONLY, inputImageFormat, inputImageSize[0], segmentationTheta.numberRows);
	cl::Image2D deviceInputImageB = cl::Image2D(context, CL_MEM_READ_ONLY, inputImageFormat, inputImageSize[0], segmentationTheta.numberRows);
	//Create sampler
	cl::Sampler sampler = cl::Sampler(context, CL_FALSE, CL_ADDRESS_CLAMP, CL_FILTER_NEAREST);
	//Create bufferss for originSensor[3], incWidth[3], incHeight[3], resolutionInputImage[2] and inputImageSize[2]
	size_t vecSize = sizeof(float) * 3;
	cl::Buffer bufferOriginSensor = cl::Buffer(context, CL_MEM_READ_ONLY, vecSize);
	cl::Buffer bufferIncWidth = cl::Buffer(context, CL_MEM_READ_ONLY, vecSize);
	cl::Buffer bufferIncHeight = cl::Buffer(context, CL_MEM_READ_ONLY, vecSize);
	cl::Buffer bufferResolutionInputImage = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * 2 );
	cl::Buffer bufferInputImageSize = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int) * 2 );
	cl::Buffer bufferSegmentRows = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int) * 2 );
	cl::size_t<3> origin;
	origin[0] = 0; origin[1] = 0; origin[2] = 0;
	cl::size_t<3> regionIn, regionOut;
	regionIn[0] = inputImageSize[0]; regionIn[1] = segmentationTheta.numberRows; regionIn[2] = 1;
	regionOut[0] = outputImageSize[0]; regionOut[1] = outputImageSize[1]; regionOut[2] = 1;
	//Enqueue ERP image
	unsigned char *ptrSegment(inputImageChannels.data() + inputImageSize[0] * segmentationTheta.index.start[segmentNo]);
	queue.enqueueWriteImage(deviceInputImageR, CL_TRUE, origin, regionIn, 0, 0, ptrSegment);
	queue.enqueueWriteImage(deviceInputImageG, CL_TRUE, origin, regionIn, 0, 0, (ptrSegment + inputSize));
	queue.enqueueWriteImage(deviceInputImageB, CL_TRUE, origin, regionIn, 0, 0, (ptrSegment +  2 * inputSize));
	//load kernel and build 
	std::ifstream sourceFile("../virtual_camera_kernelv2.0.cl");
	std::string sourceCode(std::istreambuf_iterator<char>(sourceFile), (std::istreambuf_iterator<char>()));
	cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length() +1));
	cl::Program program = cl::Program(context, source);
	program.build(devices);
	cl_int* err;
	cl::Kernel virtual_camera_kernel(program, "virtual_camera_kernel", err);
	std::cout <<"Kernel valid?:  " << err << "  " << *err << std::endl;
	//load static kernel parameters
	queue.enqueueWriteBuffer(bufferResolutionInputImage, CL_TRUE, 0, sizeof(float) * 2 , resolutionInputImage);
	queue.enqueueWriteBuffer(bufferInputImageSize, CL_TRUE, 0, sizeof(int) * 2 , inputImageSize);
	//Set kernel arguments	
	virtual_camera_kernel.setArg(0, deviceOutputImage);
	virtual_camera_kernel.setArg(1, bufferOriginSensor);
	virtual_camera_kernel.setArg(2, bufferIncWidth);
	virtual_camera_kernel.setArg(3, bufferIncHeight);
	virtual_camera_kernel.setArg(4, bufferResolutionInputImage);
	virtual_camera_kernel.setArg(5, bufferInputImageSize);
	virtual_camera_kernel.setArg(6, bufferSegmentRows);
	virtual_camera_kernel.setArg(7, deviceInputImageR);
	virtual_camera_kernel.setArg(8, deviceInputImageG);
	virtual_camera_kernel.setArg(9, deviceInputImageB);
	virtual_camera_kernel.setArg(10, sampler);
		
	//Control loop to build virtual image and control the camera
    bool proceed(true);
	bool updateImageSegment(false);
    
	while (proceed) {     
        
		clock_start = std::chrono::system_clock::now(); // STD::CHRONO CAN CAUSE CORE DUMP
	
		//generate pose parameters of the rotated camera sensor     
        sensor = RotatedSensor(phiWindow, thetaWindow, theta0, phi0);
		for (int kk=0; kk<3; kk++) {
			originSensor[kk] = sensor.origin(kk);
			incWidth[kk]  = (sensor.rightBottom(kk) - sensor.origin(kk))/(outputImageSize[0]-1);
			incHeight[kk] = (sensor.leftTop(kk) - sensor.origin(kk))/(outputImageSize[1]-1);
		}
		
		//que dynamic parameters
		queue.enqueueWriteBuffer(bufferOriginSensor, CL_TRUE, 0, vecSize, originSensor);
		queue.enqueueWriteBuffer(bufferIncWidth, CL_TRUE, 0, vecSize, incWidth);
		queue.enqueueWriteBuffer(bufferIncHeight, CL_TRUE, 0, vecSize, incHeight);
		queue.enqueueWriteBuffer(bufferSegmentRows, CL_TRUE, 0, sizeof(int) * 2 , erpRowRange);
		
		//Execute kernel and read image
		cl::NDRange global(outputImageSize[0], outputImageSize[1]);
		cl::NDRange local(16,16);
		queue.enqueueNDRangeKernel(virtual_camera_kernel, cl::NullRange, global, local);
		queue.finish(); 
		queue.enqueueReadImage(deviceOutputImage, CL_TRUE, origin, regionOut, 0, 0, outputImageVector);
		queue.finish();

		//Demux image to CImg format
		Demux(outputImageVector, outputImageSize, outputImageChannels.data());
		
		//STD::CHRONO CAN CAUSE CORE DUMP !!!!!!!!!!!!!!!!!!!!!!!!!!
		clock_interval = ceil((std::chrono::system_clock::now() - clock_start).count()/1e6); 
		time_interval = "Time per frame : " + std::to_string(clock_interval) + "  ms"; 
		outputImageChannels.draw_text(10, 10, &time_interval[0], yellow, black, 1.0, 20);
		
		//Display image
        vcDisplay.display(outputImageChannels);
        vcDisplay.wait();
      
        // Camera control
		thetaMaxPrevious = theta0 + thetaWindow/2;
		thetaMinPrevious = theta0 - thetaWindow/2;
		
        if (vcDisplay.is_keyARROWRIGHT())    phi0 +=M_PI/360;
        if (vcDisplay.is_keyARROWLEFT())     phi0 -=M_PI/360;
        if (vcDisplay.is_keyARROWUP())       theta0 +=M_PI/360;   if (theta0 >  M_PI/2) theta0 =  M_PI/2;
        if (vcDisplay.is_keyARROWDOWN())     theta0 -=M_PI/360;   if (theta0 < -M_PI/2 + thetaWindow/2) theta0 = -M_PI/2 + thetaWindow/2;
        
		if (vcDisplay.is_keyI()) {
            thetaWindow *= 0.9875;
            phiWindow   *= 0.9875;
        }
        if (vcDisplay.is_keyO() && thetaWindow < thetaWindowMax ) {
            thetaWindow *= 1.0125;
            phiWindow   *= 1.0125;
        }
		
		thetaMax = theta0 + thetaWindow/2;
		thetaMin = theta0 - thetaWindow/2;
		
        if (vcDisplay.is_keyS()) proceed = false;
		
		// Loading new segment inputImage when camera reaches the edge of the loaded segment
		updateImageSegment = false;
		for (int ii=0; ii < HORIZONTAL_SEGMENT_OPTIONS; ii++){
			if (segmentNo != 0 && (thetaMaxPrevious < M_PI * segmentationTheta.angle.start[segmentNo]/180) && (thetaMax > M_PI * segmentationTheta.angle.start[segmentNo]/180)) {
				segmentNo--;
				updateImageSegment = true;
			}
			if (segmentNo != HORIZONTAL_SEGMENT_OPTIONS && thetaMinPrevious > M_PI * segmentationTheta.angle.end[segmentNo] / 180 && thetaMin <= M_PI * segmentationTheta.angle.end[segmentNo] / 180) {
				segmentNo++;
				updateImageSegment = true;
			}
		}
		
		if (updateImageSegment) {
			erpRowRange[0] = segmentationTheta.index.start[segmentNo];
			erpRowRange[1] = segmentationTheta.index.end[segmentNo];
			ptrSegment = inputImageChannels.data() + inputImageSize[0] * segmentationTheta.index.start[segmentNo];
			queue.enqueueWriteImage(deviceInputImageR, CL_TRUE, origin, regionIn, 0, 0, ptrSegment);
			queue.enqueueWriteImage(deviceInputImageG, CL_TRUE, origin, regionIn, 0, 0, (ptrSegment + inputSize));
			queue.enqueueWriteImage(deviceInputImageB, CL_TRUE, origin, regionIn, 0, 0, (ptrSegment +  2 * inputSize));
		}
		
    }
  return 0;
}
 


 
