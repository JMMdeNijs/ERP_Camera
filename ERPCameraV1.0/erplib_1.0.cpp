

/**************************************************************************************************************************
            DEFINITIONS OF FUNCTIONS of ERP_LIBRARY
Jan de Nijs
December 2018
***************************************************************************************************************************/


#include <erplib_1.0.h>

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

void RotatedSensor(float *ptrRX, int imSize[2], int rowRange[2], float phiWindow, float thetaWindow, float phi0, float theta0){
	//RotatedSensor generates a grid/sensor located in the Z=1 plane and center at the (X,Y) = (0,0) and rotates the sensor over the angles phi0, theta0.
	const int siz(imSize[0]*imSize[1]);
	ptrRX 	+= rowRange[0]*imSize[0];      
    float *ptrRY(ptrRX +siz), *ptrRZ(ptrRY + siz);
    
	//Declaration of vectors defining 3 corners of the sensor 
	const float xMax=tan(phiWindow/2), yMax = tan(thetaWindow/2), z0(1.0);
	Vector3f r0,  r1, r2, rr0, rr1, rr2;
	float incWidth[3], incHeight[3];
	r0 << -xMax, -yMax, z0;
	r1 <<  xMax, -yMax, z0;
	r2 << -xMax,  yMax, z0;
		
    //calculate rotation matrix from phi0 (panning), theta0 (tilt)
    Matrix3f 	rM_theta, rM_phi, rMatrix;
    rM_theta 	<< 1, 0, 0, 0, cos(theta0), -sin(theta0), 0, sin(theta0), cos(theta0);
    rM_phi 		<< cos(phi0), 0, sin(phi0), 0, 1, 0, -sin(phi0), 0, cos(phi0);
    rMatrix 	= rM_phi*rM_theta;
	
	rr0 = rMatrix * r0;
	rr1 = rMatrix * r1;
	rr2 = rMatrix * r2;
	
	for (int i = 0; i<3; i++){
		incWidth[i]  = (rr1(i) - rr0(i))/(imSize[0]-1);
		incHeight[i] = (rr2(i) - rr0(i))/(imSize[1]-1);
	}
	
	float xBegin, yBegin, zBegin;
	for (int i=rowRange[0]; i <=rowRange[1]; i++){
		xBegin = rr0(0) + i*incHeight[0];
		yBegin = rr0(1) + i*incHeight[1];
		zBegin = rr0(2) + i*incHeight[2];
        for (int j=0; j < imSize[0]; j++){
            *ptrRX++ = xBegin + j*incWidth[0];
			*ptrRY++ = yBegin + j*incWidth[1];
			*ptrRZ++ = zBegin + j*incWidth[2];
			
        }
	}
}    


void Tran2SphericalCoordinates(float *ptrX, float *ptrPhi, int imSize[2], int rowRange[2]){
    //Tran2SphericalCoordinates calculates for each pixel of the sensor (as defined in Cartesian coordinates XYZ) the Spherical angles Phi and theta.  
	int siz = imSize[0]*imSize[1];
    ptrX	+= rowRange[0]*imSize[0];
	ptrPhi 	+= rowRange[0]*imSize[0];
	
	float *ptrY(ptrX + siz),  *ptrZ(ptrY + siz);
    float *ptrTheta(ptrPhi + siz);
    //iterate over pixel to calculate Phi and theta;
    for (int i=rowRange[0]; i <= rowRange[1]; i++)
        for (int j=0; j < imSize[0]; j++){
            *ptrTheta++ = -atan(*ptrY++/sqrt(*ptrX * *ptrX + *ptrZ * *ptrZ));
            *ptrPhi++ = atan2(*ptrX++, *ptrZ++);
        }
}

void CropImage(unsigned char *ptrCamImR, int imSize[2], int rowRange[2], unsigned char *ptrErpR, int erpSize[2], float *ptrPhi){
    //Given a EquiRectangular Picture (ERP), CropImage crops the ERP to the window/pixels defined by the sensor of the virtual camera given in Spherical coordinates 
    int siz = imSize[0]*imSize[1], erpSiz = erpSize[0] * erpSize[1];
    ptrCamImR   += rowRange[0]*imSize[0];
    ptrPhi      += rowRange[0]*imSize[0];
     
    unsigned char *ptrCamImG(ptrCamImR + siz),  *ptrCamImB(ptrCamImG + siz);
    unsigned char *ptrErpG(ptrErpR + erpSiz), *ptrErpB(ptrErpG + erpSiz);
     
    float *ptrTheta(ptrPhi + siz);
    const float pi(3.14159265), twopi(2*pi), halfpi = pi/2;
    const float incPhi(twopi/erpSize[0]), incTheta(pi/erpSize[1]);
    int index, indexTheta, indexPhi;
 
    for (int i=rowRange[0]; i <= rowRange[1]; i++){
        for (int j=0; j < imSize[0]; j++){
            indexPhi = round((pi + *ptrPhi++)/incPhi);
            indexTheta = round((halfpi - *ptrTheta++)/incTheta);
            index = indexTheta*erpSize[0]+indexPhi;
 
            if (index < 0) index =0;
            if (index > siz) indexPhi = siz;//????
            *ptrCamImR++ = *(ptrErpR + index);
            *ptrCamImG++ = *(ptrErpG + index);
            *ptrCamImB++ = *(ptrErpB + index);
 
        }
    }
}
    

void CameraIm(Parameters *par, int i){
	RotatedSensor(par->ptrRotated, par->camSize, par->rowThreadEntrances[i], par->phiWindow, par->thetaWindow, par->phi0, par->theta0);
	Tran2SphericalCoordinates(par->ptrRotated, par->ptrSpherical, par->camSize, par->rowThreadEntrances[i]);
	CropImage(par->ptrCamera, par->camSize, par->rowThreadEntrances[i], par->ptrImage, par->erpSize, par->ptrSpherical);
}


void ChangeCamData(Parameters &par){
    std::string aRatio("XXX");
    int aH, aV;
    std::cout << std::endl << std::endl << std::endl << "Change camera input data" << std::endl;
    std::cout << std::endl << "Specify the camera window" << std::endl;
    std::cout << "\t" << "Vertical angle of the virtual camera (degrees)  :  "; std::cin >> par.theta0;
    std::cout << "\t" << "Horizontal angle of the virtual camera (degrees):  "; std::cin >> par.phi0;
    
	while (!(aRatio == "1:1" || aRatio == "4:3" || aRatio == "16:9")) {
		std::cout << std::endl << "Specify the aspect ratio cf 1:1, 4:3, 16:9      :  ";
		std::cin >> aRatio;
	}
    
	if (aRatio.length() == 3) {
        aH = int(aRatio[0]) - 48; aV = int(aRatio[2])-48;
	}
    else {
        aH = 10 + int(aRatio[1]) - 48; aV = int(aRatio[3])-48;
	}

    std::cout << "\t" << "Horizontal resolution in pixels               :  "; std::cin >> par.camSize[0];   
    std::cout << "\t" << "Horizontal angle of view (degrees)            :  "; std::cin >> par.phiWindow;
    par.camSize[1] = par.camSize[0]/aH*aV;
	    
	par.theta0 *=M_PI/180;
    par.phi0 *=M_PI/180;
    par.phiWindow *=M_PI/180;
	par.thetaWindow = par.phiWindow*aV/aH;
	
	par.rowThreadEntrances[0][0] = 0; par.rowThreadEntrances[0][1] = par.camSize[1]/THREAD_NUMBER + par.camSize[1]%THREAD_NUMBER -1;
	for(int i=1; i<THREAD_NUMBER; i++){  
		par.rowThreadEntrances[i][0] = par.rowThreadEntrances[i-1][0] + par.camSize[1]/THREAD_NUMBER; 
		par.rowThreadEntrances[i][1] = par.rowThreadEntrances[i-1][1] + par.camSize[1]/THREAD_NUMBER;
	}

}

	
