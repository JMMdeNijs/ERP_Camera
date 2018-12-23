cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project("ERP Camera")

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package (Threads REQUIRED)

find_package(X11 REQUIRED)

message(STATUS "X11_FOUND = ${X11_FOUND}")
message(STATUS "X11_INCLUDE_DIR = ${X11_INCLUDE_DIR}")
message(STATUS "X11_LIBRARIES = ${X11_LIBRARIES}")


include_directories(
			~/Libraries/CImg
			~/Libraries/Eigen32
			${X11_INCLUDE_DIR})

add_executable(	ERPCamera ERPCamera1.0.cpp
				erplib_1.0.cpp)
			
target_link_libraries(ERPCamera  Threads::Threads ${X11_LIBRARIES})
