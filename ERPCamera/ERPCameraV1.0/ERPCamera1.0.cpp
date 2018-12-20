/****************************************************************************************
This code captures a regular 2D image from a Equirectangular Picture (ERP) using a virtual camera. 

The ERP image is projected on a sphere. The virtual camera is positioned in the centre of the sphere 
and can be oriented in the vertical (inclination - theta) and horizontal (azimuth - PHI) directions using the Arrow Keys
The virtual camera can be zoomed in and out using the Keys i/I and o/O
The aspect ratio can be changed using the Keys c/C
The camera can be stopped using the Key s/S.

To accelerate the code, threads are used for parallel processing.

Copyright (C), December 2018, Jan de Nijs

*****************************************************************************************/

#include "../erplib_1.0.h"


int main(int argc,char **argv) {
    
	std::string imDir(GetFile());
	//string imDir ="D:/Projects/ERPCamera/Laon_Cathedral.jpg";
    
	//Load and display equirectangle image
    CImg<unsigned char> image(&imDir[0]);
    CImgDisplay main_disp, sec_disp;
    main_disp.display(image);

	//Declaration and initialization camera parameter field (struct) for threading
	Parameters camera;
	camera.erpSize[0] = image.width();
	camera.erpSize[1] = image.height();
    camera.theta0 =0.;
    camera.phi0 = 0.;
    camera.thetaWindow = M_PI*39.375/180;
    camera.phiWindow = M_PI*70.0/180;
	camera.camSize[0] = 1024; camera.camSize[1] = 576;
	
	//Timing & Clock
	std::chrono::time_point<std::chrono::system_clock> clock_start, clock_end; 
	int clock_interval;
	char *textstring;
	std::string time_interval;
	
	//Colors information in display
	unsigned char yellow[3] = {0, 255, 255}, black[3] = {0, 0, 0};
	
	//Show-hide info bar (help)
	bool show_hide(true);
	
	//Define segmentation for threads
	camera.rowThreadEntrances[0][0] = 0; camera.rowThreadEntrances[0][1] = camera.camSize[1]/THREAD_NUMBER + camera.camSize[1]%THREAD_NUMBER -1;
    for(int i=1; i<THREAD_NUMBER; i++){  
        camera.rowThreadEntrances[i][0] = camera.rowThreadEntrances[i-1][0] + camera.camSize[1]/THREAD_NUMBER;
        camera.rowThreadEntrances[i][1] = camera.rowThreadEntrances[i-1][1] + camera.camSize[1]/THREAD_NUMBER;
    }
        
	//Declaration camera image variables
    CImg<float> sensorXYZ_Rotated(camera.camSize[0], camera.camSize[1],1,3,0.0);
	CImg<float> sensorSphericalCoordinates(camera.camSize[0], camera.camSize[1],1,2,0.0);
	CImg<unsigned char> cameraImage(camera.camSize[0], camera.camSize[1],1,3,0.0);
	
	//Completion camera parameter field for threading
	camera.ptrRotated 	= sensorXYZ_Rotated.data();
	camera.ptrSpherical = sensorSphericalCoordinates.data();
	camera.ptrImage		= image.data();
	camera.ptrCamera	= cameraImage.data();
			
	std::vector<std::thread> threadList;

    bool proceed(true);
    		
    while (proceed) {
		
		clock_start = std::chrono::system_clock::now(); 
		
		//generate and display camera image		
		for (int i=0; i<THREAD_NUMBER; i++) 	threadList.push_back(std::thread(CameraIm, &camera, i));
        for (int i=0; i<THREAD_NUMBER; i++)		(threadList.at(i)).join();
        threadList.erase (threadList.begin(),threadList.end());	
		
		//show time interval in image
		clock_interval = ceil((std::chrono::system_clock::now() - clock_start).count()/1e6); 
		time_interval = "Time per frame : " + std::to_string(clock_interval) + "  ms";
		cameraImage.draw_text(10, 10, &time_interval[0], yellow, black, 1.0, 20);
		
		//Legend control and display
		if (show_hide) {
			cameraImage.draw_text(10, camera.camSize[1] - 30, "Hide/Show legend: H,   Zoom In/Out:  i/o ", yellow, black, 1.0, 20);
			cameraImage.draw_text(450, camera.camSize[1] - 30, "Camera rotation: Arrows     Change camera aspects:   C", yellow, black, 1.0, 20);
		}
		
		sec_disp.display(cameraImage);
		sec_disp.wait();
		
		// Camera control		
		if (sec_disp.is_keyARROWRIGHT()) 	camera.phi0 +=M_PI/360;
		if (sec_disp.is_keyARROWLEFT()) 	camera.phi0 -=M_PI/360;
		if (sec_disp.is_keyARROWUP()) 		camera.theta0 +=M_PI/360;	if (camera.theta0 >  M_PI/2) camera.theta0 =  M_PI/2;
		if (sec_disp.is_keyARROWDOWN()) 	camera.theta0 -=M_PI/360;	if (camera.theta0 < -M_PI/2 + camera.thetaWindow/2) camera.theta0 = -M_PI/2 + camera.thetaWindow/2;
		
		if (sec_disp.is_keyI()) {
			camera.thetaWindow 	*= 0.9875;
			camera.phiWindow	*= 0.9875;
		}
		
		if (sec_disp.is_keyO()) {
			camera.thetaWindow 	*= 1.0125;
			camera.phiWindow	*= 1.0125;
		}
		
		if (sec_disp.is_keyS()) proceed = false;

		if (sec_disp.is_keyC()) {
			ChangeCamData(camera);
			sec_disp.resize(camera.camSize[0], camera.camSize[1]);
					
			sensorXYZ_Rotated.assign(camera.camSize[0], camera.camSize[1],1,3,0.0);				camera.ptrRotated 	= sensorXYZ_Rotated.data();
			sensorSphericalCoordinates.assign(camera.camSize[0], camera.camSize[1],1,2,0.0);	camera.ptrSpherical = sensorSphericalCoordinates.data();
			cameraImage.assign(camera.camSize[0], camera.camSize[1],1,3,0.0);					camera.ptrCamera	= cameraImage.data();
		}		
		
		if (sec_disp.is_keyH()) show_hide = !show_hide;
    }
  return 0;
}



