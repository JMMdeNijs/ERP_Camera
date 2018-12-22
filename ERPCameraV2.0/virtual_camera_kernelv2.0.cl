// Calculates the polar angle or inclination (theta) and azimuth (phi)of a pixel in a grid that is rotated over
// an angle (phi0, theta0). Input are the 3 vectors that defines the origin, and the incremental 
// width and hight along the left and bottom sides per pixel of the rotated grid, given in Euclidian coordinates.
// The reference grid is centered at (0,0,1) and parallel to the XY plane with:
//     origin  =	  (-Xo, -Yo, 1)
//     left top = 	  (-Xo,  Yo, 1)
//     right bottom = ( Xo, -Yo, 1)
// The spherical coordinates theta and phi are calculated and mapped on the index space of the ERP image. 
// The RGB values are copied from the ERP image to the image of the virtual camera.
// Jan de Nijs, December 2018
 
__kernel
void virtual_camera_kernel(
	__write_only	image2d_t cameraImage,
	__constant		float	*origin,
	__constant		float	*inc_width,
	__constant		float	*inc_height,
	__constant		float	*resolutionERP,
	__constant		int		*erpSize,
	__constant		int		*erpRowRange,
	__read_only		image2d_t erpImageR,
	__read_only		image2d_t erpImageG,
	__read_only		image2d_t erpImageB,
					sampler_t sampler		) {					
	
	int col	= get_global_id(0);
	int row	= get_global_id(1);
	int2 coordinate = {col, row};
	
	//Calculation of the pixel coordinate in Euclidean coordinates
	float3 pixelCoordinate;
	pixelCoordinate.x = origin[0] + col*inc_width[0] + row*inc_height[0];
	pixelCoordinate.y = origin[1] + col*inc_width[1] + row*inc_height[1];
	pixelCoordinate.z = origin[2] + col*inc_width[2] + row*inc_height[2];

	//Conversion from Euclidean coordinates to polar coordinates (azimuth and polar angle) of the ERP image
	float phi, theta;
	phi = atan2(pixelCoordinate.x , pixelCoordinate.z);  //azimuth
	theta = -atan(pixelCoordinate.y / sqrt(pixelCoordinate.x * pixelCoordinate.x + pixelCoordinate.z * pixelCoordinate.z)); //inclination
	
	//Find index of pixel corresponding to azimuth and inclination angles
	int2 indexPhiTheta;
	indexPhiTheta.x = round((M_PI + phi) / resolutionERP[0]); //* erpSize[0] / 2 / M_PI);
	if(indexPhiTheta.x >= erpSize[0]) indexPhiTheta.x = erpSize[0] -1;
    indexPhiTheta.y = round((M_PI_2 - theta) / resolutionERP[1]) - erpRowRange[0]; //* erpSize[1] / M_PI); 
	
	//Find RGB values of pixel azimuth phi and inclination theta 
	uint4 pixelColor;
	pixelColor.x = read_imageui(erpImageR, sampler, indexPhiTheta).x; 
	pixelColor.y = read_imageui(erpImageG, sampler, indexPhiTheta).x; 
	pixelColor.z = read_imageui(erpImageB, sampler, indexPhiTheta).x;
	pixelColor.w = 0;
	
	write_imageui(cameraImage, coordinate, pixelColor);
}
