/**************************************************************************************************************************
            DECLARATIONS OF FUNCTIONS of ERP_LIBRARY Version 1.0
			
Copyright (C), December 2018, Jan de Nijs

***************************************************************************************************************************/

#ifndef __ERPLIBv10__
#define __ERPLIBv10__

#include <CImg.h>
#ifdef Success
#undef Success
#endif

#include <thread>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <chrono>
#include <ctime>


using Eigen::Matrix3f;
using Eigen::Vector3f;
using namespace cimg_library;

#define THREAD_NUMBER   4
#define M_2PI			2*M_PI

struct Parameters{
		float       *ptrSensor;
		float       *ptrRotated;
		float       *ptrSpherical;
		unsigned char   *ptrImage;
		unsigned char	*ptrCamera;
		int         rowThreadEntrances[THREAD_NUMBER][2];
		int         erpSize[2];
		int			camSize[2];
		float		theta0;
		float		thetaWindow;
		float		phi0;
		float		phiWindow;
	};   

std::string GetFile(void);
//Collects file name and folder of the ERP image

void RotatedSensor(float *ptrRX, int const imSize[2], int const rowRange[2], float const phiWindow, float const thetaWindow, float const phi0, float const theta0);//RotatedSensor generates a grid/sensor located in the Z=1 plane and center at the (X,Y) = (0,0) and rotates the sensor over the angles phi0(azimuth), theta0(inclination).

void Tran2SphericalCoordinates(float *ptrX, float *ptrPHI, int const imSize[2], int const rowRange[2]);
//Tran2SphericalCoordinates calculates for each pixel of the rotated sensor (as defined in Cartesian coordinates XYZ) the Spherical angles phi and theta.  

void CropImage(unsigned char *ptrCamImR, int const imSize[2], int const rowRange[2], unsigned char *ptrERP_R, int const erpSize[2], float *ptrPHI);
//CropImage crops the ERP to the window/pixels defined by the sensor of the virtual camera given in spherical coordinates phi and theta.

void CameraIm(Parameters *par, int i);
//Aggregation of function called by threads 

void ChangeCamData(Parameters &par);
#endif
