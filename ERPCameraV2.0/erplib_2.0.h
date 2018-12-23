/******************************************************************************************************

			DECLARATIONS OF FUNCTIONS of ERP LIBRARY Version 2.0
			
Copyright (C), December 2018, Jan de Nijs

*******************************************************************************************************/
#ifndef __ERPLIBv20__
#define __ERPLIBv20__

#include <CImg.h>
#ifdef Success
#undef Success
#endif

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <chrono>
#include <ctime>
#include </CL/cl.hpp>

using Eigen::Matrix3f;
using Eigen::Vector3f;
using namespace cimg_library;

#define THREAD_NUMBER   3
#define HORIZONTAL_SEGMENT_OPTIONS  5  //Must be an odd number and 5 or larger!!!
#define __CL_ENABLE_EXCPTIONS


struct SensorPose {
	Eigen::Vector3f		origin;
	Eigen::Vector3f		leftTop;
	Eigen::Vector3f		rightBottom;
};

struct WindowSegmentation {
	int numberRows;
	struct  {
		int start[HORIZONTAL_SEGMENT_OPTIONS];
		int	end[HORIZONTAL_SEGMENT_OPTIONS];
	} index;

	struct  {
		float start[HORIZONTAL_SEGMENT_OPTIONS];
		float	end[HORIZONTAL_SEGMENT_OPTIONS];
	} angle;

};


void Demux(unsigned char *ptrMuxedImage, int const imageSize[2], unsigned char *ptrImR);
//demuxes the RGBA vectors into RGB channels of the CImg image format

std::string GetFile(void);
//utility function to read folder and name ERP image

SensorPose RotatedSensor(float const phiWindow, float const thetaWindow, float const theta0, float const phi0);
//Calculation of vectors defining 3 corners of the sensor after rotation
//phi0(azimuth) and theta0 (inclination) specify the normal of the sensor.
//RotatedSensor.origin, RotatedSensor.leftTop and RotatedSensor.rightBottom specify the Euclidean vectors defining the sensor. 
#endif
