/**************************************************************************************************************************
            DEFINITIONS OF FUNCTIONS of ERP LIBRARY Version 2.0
			
			
Copyright (C), December 2018, Jan de Nijs
***************************************************************************************************************************/
 
#include <erplib_2.0.h>


void Demux(unsigned char *ptrMuxedImage, int imageSize[2], unsigned char *ptrImR){
	//demuxes the RGBA vectors into RGB channels of a CImg image
	int sizeImage(imageSize[0]*imageSize[1]);	
	unsigned char *ptrImG(ptrImR + sizeImage), *ptrImB(ptrImG + sizeImage);
	for (int ii=0; ii < sizeImage; ii++){
		*ptrImR++ = *ptrMuxedImage++;
		*ptrImG++ = *ptrMuxedImage++;
		*ptrImB++ = *ptrMuxedImage++;
		*ptrMuxedImage++;
	}
}



std::string GetFile(void){
	std::string fileDir, fileName;
	bool found(false);
	while (!found) {
		std::cout << std::endl << "Give image name              : "; std::cin >> fileName;
		std::cout << std::endl << "Give path to image directory : "; std::cin >> fileDir;
		for(unsigned int i=0; i<fileDir.length(); i++)
			if (fileDir[i] ==92) fileDir[i]=47; //replaces '\'  by '/'
		fileDir+= ("/" + fileName);
		std::ifstream filecheck(fileDir);
		found = (bool)filecheck;
		found? std::cout << "Image found \n" : std::cout << "Image not found \n\n" ;
	}
	
	return fileDir;
}


SensorPose RotatedSensor(float phiWindow, float thetaWindow, float theta0, float phi0){
	const float xMax=tan(phiWindow/2), yMax = tan(thetaWindow/2), z0(1.0);
    SensorPose rotatedSensor;
	Vector3f r0,  r1, r2;
    float Inc_Width[3], Inc_Height[3];
    r0 << -xMax, -yMax, z0; //origin
    r1 <<  xMax, -yMax, z0; //right_bottom
    r2 << -xMax, +yMax, z0; //left_top
         
    //calculate rotation matrix from phi0 (panning), theta0 (tilt)
    Matrix3f    rM_theta, rM_phi, rMatrix;
    rM_theta     << 1, 0, 0, 0, cos(theta0), -sin(theta0), 0, sin(theta0), cos(theta0);
    rM_phi       << cos(phi0), 0, sin(phi0), 0, 1, 0, -sin(phi0), 0, cos(phi0);
    rMatrix    = rM_phi*rM_theta;
    
	//Calculate origin, leftTop en rightBottom in Euclidean coordinates
    rotatedSensor.origin = rMatrix * r0;
    rotatedSensor.leftTop = rMatrix * r2;
    rotatedSensor.rightBottom = rMatrix * r1;
	
	return rotatedSensor;
} 
   
 





