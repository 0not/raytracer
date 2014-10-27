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

namespace kgl {
    class camera {
    public:
        pos3 pos;
        pos3 up;
        pos3 look;

        camera() : pos(), up(), look() {
            
        }
    };

};
int main() {
    // Image dimensions
    const int width  = 600,
              height = 400;
    // Can also set color bits:
    // kgl::RGB_BITS = 8; // default value (try 4 or 2)

    // Image bitmap
    // 2-d array that stores a color vector for each element
    static kgl::raster_t<width, height> raster;
    // Set a background color
    kgl::rgb background_color = kgl::black;

    // Camera setup
    kgl::camera cam;
    cam.pos.set(6, 0, 0);               // Position of the camera
    cam.up.set(0, 0, 1.0);              // Up vector
    cam.look.set(0.0, 0.0, 0.0);        // Where the camera is pointing
    cam.look = (cam.look - cam.pos).normalized(); // Normalized vector in direction of cam.look


    // Focal distance (cam.pos to sensor)
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
    kgl::space_t fov = atan(dsensor/(2*focal_dist));    // Half angle field of view (e.g. left to center)
    kgl::space_t fov2 = fov * 2.0;                      // Full angle field of view (e.g. left to right, total field)

    // center of image: (cam.pos + focal_dist * cam.look)
    // I think of this as walking from cam.pos to cam.look for a distance of focal_dist
    kgl::pos3 i_center = (cam.pos + focal_dist * cam.look);  

    // Left vector
    // If looking at i_center and up is cam.up, then left is... left!
    // Thanks to the magic of cross products, this is easy!
    kgl::pos3 left =  i_center.cross(cam.up).normalized();

    // top left: (||i_center x cam.up|| + cam.up) * raperature + i_center 
    kgl::pos3 top_left = (left * wsensor + cam.up * hsensor) / 2.0 + i_center - cam.pos;

    // Differential of each pixel
    // Size of a pixel along each dimension on the "sensor"
    // Used to construct the ray from the eye to the appropriate sensor location
    kgl::space_t dc = wsensor / width;
    kgl::space_t dr = hsensor / height;

    // Make sphere
    const kgl::space_t radius  = 1.0;
    const kgl::space_t radius2 = radius*radius;
    kgl::pos3 center{0.0, 0.0, 0.0};

    // light
    kgl::pos3 light{2.0, 0.0, 1};


    kgl::pos3 ray;       // Direction of ray
    kgl::pos3 to_obj;    // From camera position to center of sphere
    for (int r = 0; r < height; r++) {
        for (int c = 0; c < width; c++) {
            // Shoot a ray from the camera
            ray = (top_left - left * c * dc - cam.up * r * dr).normalized();
            
            // Check for intersections
            to_obj = (center - cam.pos);

            // Amount of "to_obj" in direction or "ray"
            // This is the parameter in the vector equation: ray = ray_orig + ray_dir * tca
            // It means how for down the ray you have to walk
            kgl::space_t tca = ray.dot(to_obj);

            // If less than 0, ray is behind us
            if (tca < 0) {
                raster[r][c] = background_color;
                continue; // No intersection
            }

            // d^2 + tca^2 = ||to_obj||^2
            // d2 is the closest approach of the ray to the obj
            // We don't solve for d because that would mean expensive square root operations that aren't necessary
            // (since we can just compare to the radius squared)
            kgl::space_t d2 = to_obj.norm2() - tca * tca;

            // If the ray never gets closer than the radius, we know there is no intersection
            if (d2 > radius2) {
                raster[r][c] = background_color;
                continue; // No intersection
            }
            
            // The ray intersects the sphere, going through two points (think entrance and exit wound)
            // tca is how far from the ray orign to the closest approach of the center.
            // Add thc to go to the exit, or subtract thc to go to the entrance.
            kgl::space_t thc = sqrt(radius2 - d2),
                         t0  = tca - thc,           // Closest intersection
                         t1  = tca + thc;           // Farthest intersection

            // Just grab the closest one
            kgl::space_t t = t0;

            // This is the surface normal at the intersection point
            // It is used to calculate the diffuse shading
            kgl::pos3 normal = ((cam.pos + t*ray) - center).normalized();

            // This is a value that indicates how much of the normal is pointing in the direction of the light
            // The more that's in the lights direction, the brighter that spot on the sphere.
            kgl::color_t diffuse = (normal.dot(light.normalized()));
            
            // We need to "diffuse" to be from 0 to 1, so that it is a valid color value
            if (diffuse < 0)
                diffuse = 0;
            else if (diffuse > 1)
                diffuse = 1;
            
            // Set the RGB color
            kgl::rgb color{1, 0, 0};

            // Multiply color by the diffuse value to get the correct diffuse shading
            color *= diffuse;

            // Set the raster value
            raster[r][c] = color;

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

        // Output the PPM header
        //    P3
        //    width height
        //    max_color_value
        f << "P3\n"  <<                           // P3
             width   << " " << height << "\n" <<  // 800 600
             max_rgb << "\n";                     // 255

        // For each pixel, output the color values
        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                // Set the value of raster to match the type needed in the PPM file, with correct bit value
                raster[r][c] = raster[r][c] * max_rgb / kgl::RGB_MAX;

                // Each value must be separated with whitespace
                // I chose to output each pixel on its own line. 
                // This is a red pixel:
                //     255 0 0
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
