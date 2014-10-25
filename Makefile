raytracer: raytrace.cpp
	g++ -o raytrace raytrace.cpp --std=c++1y -Wall -Wextra && ./raytrace && ppmtobmp image.ppm > image.bmp
