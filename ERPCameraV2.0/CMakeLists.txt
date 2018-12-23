cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project("ERP Camera")

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package (Threads REQUIRED)

find_package(X11 REQUIRED)
message(STATUS "X11_FOUND = ${X11_FOUND}")
message(STATUS "X11_INCLUDE_DIR = ${X11_INCLUDE_DIR}")
message(STATUS "X11_LIBRARIES = ${X11_LIBRARIES}")

# Handle OpenCL
find_package(OpenCL REQUIRED)
include_directories(${OpenCL_INCLUDE_DIRS})
#set_target_properties (OpenCL PROPERTIES VERSION "1.2" SOVERSION "1")

link_directories(${OpenCL_LIBRARY})

include_directories(	~/Libraries/CImg
						~/Libraries/Eigen32)

add_executable(ERPCamera 
				ERPCamera_2.0.cpp
				erplib_2.0.cpp		)
			
target_link_libraries (	ERPCamera
						${OpenCL_LIBRARY}
						Threads::Threads
						${X11_LIBRARIES})

configure_file(virtual_camera_kernelv2.0.cl ${CMAKE_CURRENT_BINARY_DIR}/virtual_camera_kernelv2.0.cl COPYONLY)

