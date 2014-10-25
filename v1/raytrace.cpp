#include <array>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <complex>
#include "kmath.hpp"
#include "ktypes.hpp"
#include "kdebug.hpp"


// Declare the image save function as a template function
template <unsigned short W, unsigned short H>
int save_image(kgl::raster_t<W, H>&, const std::string& = "image.ppm");

int main() {
    // Image dimensions
    const int width  = 400,
              height = 400;
    // Can also set color bits:
    // kgl::RGB_BITS = 8; // default value (try 4 or 2)

    // Image bitmap
    // 2-d array that stores a color vector for each element
    static kgl::raster_t<width, height> raster;
    // Set a background color
    kgl::rgb background_color = kgl::black;

    // Camera setup
    kgl::pos3 c_pos{6, 0, 0};               // Position of the camera
    kgl::pos3 c_up{0, 0, 1.0};              // Up vector
    kgl::pos3 c_look{0.0, 0.0, 0.0};        // Where the camera is pointing
    c_look = (c_look - c_pos).normalized(); // Normalized vector in direction of c_look

    // Focal distance (c_pos to sensor)
    // Like in a real camera, increase the number to zoom in 
    // This focal distance is really half the length of a focal distance describing a real camera system
    kgl::space_t focal_dist = 0.04;  // In meters (20mm = 0.02m) 

    // Sensor description
    /*  w = width; h = height
     *  x = wsensor; y = hsensor;
     *  d = dsensor
     *  __x__   
     * |\    |   d^2 = x^2 + y^2
     * |  d  y   x/y = w/h      
     * |____\|   x = y*w/h
     * 
     * d^2 = y^2 * w^2 / h^2 + y^2
     * d^2 = y^2 * (w^2/h^2 + 1)
     * y^2 = d^2 / (w^2/h^2 + 1)  => y = d / sqrt(w^2/h^2 + 1)
     *
     * x^2 = d^2 - y^2  => x = sqrt(d^2 - y^2)
     */
    kgl::space_t dsensor = 0.035;                                                           // diagonal of sensor
    kgl::space_t hsensor = dsensor / sqrt(width*width/(kgl::space_t)(height*height) + 1);   // sensor height
    kgl::space_t wsensor = sqrt(dsensor*dsensor - hsensor*hsensor);                         // sensor width
    
    // Field of view
    // To zoom in, have a narower field of view
    // I opted to have fov be a computed value, since this raytracer is supposed to mimic a camera.
    // To have a wider field, increase the sensor size and/or decrease the focal length
    // To have a narrower fied, decrease the sensor size and/or increase the focal length
    // field of view: 2*atan(x/(2*f)) where x is diagonal of sensor
    kgl::space_t fov = atan(dsensor/(2*focal_dist));
    kgl::space_t fov2 = fov * 2.0;

    // center of image: (c_pos + focal_dist * c_look)
    // I think of this as walking from c_pos to c_look for a distance of focal_dist
    kgl::pos3 i_center = (c_pos + focal_dist * c_look);  

    // Left vector
    // If looking at i_center and up is c_up, then left is... left!
    // Thanks to the magic of cross products, this is easy!
    kgl::pos3 left =  i_center.cross(c_up).normalized();

    // top left: (||i_center x c_up|| + c_up) * raperature + i_center 
    kgl::pos3 top_left = (left * wsensor + c_up * hsensor) / 2.0 + i_center - c_pos;


    // Print useful info
    std::cout << "dsensor: " << dsensor << std::endl;
    std::cout << "wsensor: " << wsensor << std::endl;
    std::cout << "hsensor: " << hsensor << std::endl;
    std::cout << "i_center: " << printv(i_center) << std::endl;
    std::cout << "dist (center): " << (i_center - c_pos).norm() << std::endl;
    std::cout << "left: " << printv(left) << std::endl;
    std::cout << "top_left: " << printv(top_left) << std::endl;
    std::cout << "c->tl: " << printv(top_left - i_center) << std::endl;
    std::cout << "dist (c->tl): " << (top_left - i_center).norm() << std::endl;
    std::cout << "fov2: " << fov2 << "; " << fov2 * 180 / (2 * 3.14159) << std::endl;

    kgl::space_t dc = wsensor / width;
    kgl::space_t dr = hsensor / height;


    std::cout << "(" << dc << ", " << dr << ")" << std::endl;
    
    // Make sphere
    const kgl::space_t radius = 1.0;
    kgl::pos3 center{0.0, 0.0, 0.0};

    // light
    kgl::pos3 light{2.0, 0.0, 1};


    kgl::pos3 ray;
    kgl::pos3 L;
    for (int r = 0; r < height; r++) {
        for (int c = 0; c < width; c++) {
            // Shoot a ray from the camera
            ray = (top_left - left * c * dc - c_up * r * dr);
            ray.normalize();
            L = (center - c_pos);
            kgl::space_t tca = ray.dot(L);
            if (tca < 0) {
                raster[r][c] = background_color;
                continue; // No intersection
            }
            kgl::space_t d2 = L.dot(L) - tca * tca;
            //std::cout << "d2: " << d2 << "; r2: " << radius*radius << std::endl;
            if (d2 > radius * radius) {
                raster[r][c] = background_color;
                continue; // No intersection
            }
            kgl::space_t thc = sqrt(radius * radius - d2);
            kgl::space_t t0 = tca - thc;
            kgl::space_t t1 = tca + thc;
            //std::cout << "tca, thc: (" << tca << ", " << thc << ")" << std::endl;
            //std::cout << "t: (" << t0 << ", " << t1 << ")" << std::endl;
            kgl::space_t t = t0;
            if (t1 < t0) 
                t = t1;
            //kgl::pos3 normal = (t*ray - center).normalized();
            kgl::pos3 normal = ((c_pos + t*ray) - center).normalized();
            //std::cout << "normal: " << printv(normal) << std::endl;
            kgl::color_t diffuse = (normal.dot(light.normalized()));
            if (diffuse < 0)
                diffuse = 0;
            else if (diffuse > 1)
                diffuse = 1;
            //std::cout << "diffuse: " << diffuse << std::endl;
            kgl::rgb color{diffuse, 0, 0};
            //std::cout << "color: " << printv(color) << std::endl;
            raster[r][c] = color;
                
        

            // Check for intersections
        }
    }

    // Save image
    std::cout << "Saving image..." << std::endl;
    return save_image<width, height>(raster);
}


template <unsigned short width, unsigned short height>
int save_image(kgl::raster_t<width, height> &raster, const std::string &filename) {
    std::ofstream f(filename);
    if (f.is_open()) {
        // Calculated maximum RGB value, based on number of bits per value to store in file
        // Must calculate seperately from kgl::RGB_MAX, to allow for 0 -> 1 representations for storage in raster
        // (e.g. raster can be 0 -> 1 (floating point) while the max value to store is 255 (8 bit)
        unsigned int max_rgb = pow(2, kgl::RGB_BITS) - 1;
        f << "P3\n"  <<                           // P3
             width   << " " << height << "\n" <<  // 800 600
             max_rgb << "\n";                     // 255

        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                raster[r][c] = raster[r][c] * max_rgb / kgl::RGB_MAX;
                f << (unsigned int)raster[r][c].x << " " << 
                     (unsigned int)raster[r][c].y << " " << 
                     (unsigned int)raster[r][c].z << "\n"; 
            }
        }

        f.close();
    } else {
        std::cerr << "Unable to open file." << std::endl;
        return -1;
    }
    return 0;
}
