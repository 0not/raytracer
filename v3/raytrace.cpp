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
    pos3 random_in_unit_sphere();

    class camera {
    public:
        pos3 pos;
        pos3 up;
        pos3 look;

        camera() : pos(), up(), look() {
            
        }
    };

    class ray {
        public:
            ray() {}
            ray(const pos3 &o, const pos3 &d)    { origin = o; direction = d; }
            pos3 point_at(space_t t) const       { return origin + t*direction; }
            space_t length() const               { return (origin - direction).length(); }

            pos3 origin;
            pos3 direction;
    };

    // hittable.h
    class material;

    struct hit_record {
        space_t t;
        pos3 point;
        pos3 normal;
        material *mat_ptr;
    };

    class hittable {
        public:
            virtual bool hit(const ray &r, space_t t_min, space_t t_max, hit_record &rec) const = 0;
    };

    // hittablelist.h
    class hittable_list : public hittable {
        public:
            hittable_list() {}
            hittable_list(hittable **l, int n) { list = l; list_size = n; }
            virtual bool hit(const ray &r, space_t t_min, space_t t_max, hit_record &rec) const;

            hittable **list;
            int list_size;
    };

    bool hittable_list::hit(const ray &r, space_t t_min, space_t t_max, hit_record &rec) const {
        hit_record temp_rec;
        bool hit_anything = false;
        space_t closest_so_far = t_max;
        for (int i = 0; i < list_size; i++) {
            if (list[i]->hit(r, t_min, closest_so_far, temp_rec)) {
                hit_anything   = true;
                closest_so_far = temp_rec.t;
                rec            = temp_rec;
            }
        }
        
        return hit_anything;
    }

    // sphere.h
    class sphere : public hittable {
        public:
            sphere() {}
            sphere(pos3 cen, space_t r, material *mat) : center(cen), radius(r), mat_ptr(mat) {};
            virtual bool hit(const ray &r, space_t t_min, space_t t_max, hit_record &rec) const;
            
            pos3      center;
            space_t   radius;
            material *mat_ptr;
    };

    class material {
        public:
            virtual bool scatter(const ray &r_in, const hit_record &rec, rgb &attenuation, ray &scattered) const = 0;
    };

    class lambertian : public material {
        public:
            lambertian(const rgb & a) : albedo(a) {}
            virtual bool scatter(const ray &r_in, const hit_record &rec, rgb &attenuation, ray &scattered) const {
                pos3 target = rec.normal + random_in_unit_sphere();
                scattered   = ray(rec.point, target);
                attenuation = albedo;
                return true;
            }
            rgb albedo;
    };

    class metal : public material {
        public:
            metal(const rgb & a) : albedo(a) { fuzz = 0; }
            metal(const rgb & a, float f) : albedo(a) { 
                if (f < 1) 
                    fuzz = f; 
                else fuzz = 1;
            }
            virtual bool scatter(const ray &r_in, const hit_record &rec, rgb &attenuation, ray &scattered) const {
                pos3 reflected = reflect(r_in.direction.normalized(), rec.normal);
                scattered      = ray(rec.point, reflected + fuzz*random_in_unit_sphere());
                attenuation    = albedo;
                return (scattered.direction.dot(rec.normal) > 0);
            }

            rgb albedo;
            float fuzz;
    };

    class dielectric : public material {
        public:
            dielectric(float ri) : ref_idx(ri) {}
            virtual bool scatter(const ray &r_in, const hit_record &rec, rgb &attenuation, ray &scattered) const {
                pos3 outward_normal;
                pos3 reflected = reflect(r_in.direction, rec.normal);
                float ni_over_nt;
                attenuation = rgb(1.0, 1.0, 1.0);
                pos3 refracted;
                float reflect_prob;
                float cosine;
                if (r_in.direction.dot(rec.normal) > 0) {
                    outward_normal = -rec.normal;
                    ni_over_nt     = ref_idx;
                    cosine = ref_idx * r_in.direction.dot(rec.normal) / r_in.direction.length();
                } else {
                    outward_normal = rec.normal;
                    ni_over_nt     = 1.0 / ref_idx;
                    cosine = -r_in.direction.dot(rec.normal) / r_in.direction.length();
                }
                if (refract(r_in.direction, outward_normal, ni_over_nt, refracted)) {
                    reflect_prob = schlick(cosine, ref_idx);
                } else {
                    reflect_prob = 1.0;
                }
                if (drand48() < reflect_prob) {
                    scattered = ray(rec.point, reflected);
                } else {
                    scattered = ray(rec.point, refracted);
                }
                return true;
            }

            float ref_idx;
    };

    //space_t hit_sphere(const pos3 &center, space_t radius, const ray &r);

    rgb get_color(const ray &r, hittable *world, int depth);
};

// Set a background color
kgl::rgb background_color = kgl::rgb(0.0, 1.0, 1.0);
// light
kgl::pos3 light{20.0, 0.0, 5.0};

int main() {
    // Image dimensions
    const int width  = 600*3,
              height = 400*3;

    const int num_alias_rays = 1000;
    // Can also set color bits:
    //kgl::RGB_BITS = 2; // default value is 8 (try 4 or 2)

    // Image bitmap
    // 2-d array that stores a color vector for each element
    static kgl::raster_t<width, height> raster;

    // Camera setup
    kgl::camera cam;
    cam.pos.set(20.0, 0, 5.0);           // Position of the camera
    cam.up.set(-5.0, 0.0, 20.0);              // Up vector
    cam.up.normalize();
    cam.look.set(0.0, 0.0, 0.0);        // Where the camera is pointing
    cam.look = (cam.look - cam.pos).normalized(); // Normalized vector in direction of cam.look


    // Focal distance (cam.pos to sensor)
    // Like in a real camera, increase the number to zoom in 
    // This focal distance is really half the length of a focal distance describing a real camera system
    kgl::space_t focal_dist = 0.04;  // In meters (20 mm = 0.02 m) 

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
    //kgl::space_t fov = atan(dsensor/(2*focal_dist));    // Half angle field of view (e.g. left to center)
    //kgl::space_t fov2 = fov * 2.0;                      // Full angle field of view (e.g. left to right, total field)

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
    //const kgl::space_t radius  = 1.0;
    //const kgl::space_t radius2 = radius*radius;
    //kgl::pos3 center{0.0, 1.0, 0.0};
    //kgl::sphere sp(center, radius);

    int num_hittables = 10;
    kgl::hittable *list[num_hittables];
    list[0] = new kgl::sphere(kgl::pos3(2, 4, 0.5), 0.5, new kgl::lambertian(kgl::rgb(0.8, 0.3, 0.3)));
    list[1] = new kgl::sphere(kgl::pos3(2, 0, 1), 1.0, new kgl::lambertian(kgl::rgb(1.0, 0.1, 0.1)));
    list[2] = new kgl::sphere(kgl::pos3(0, 0, -1000.0), 1000.0, new kgl::lambertian(kgl::rgb(0.0, 0.2, 0.8)));
    list[3] = new kgl::sphere(kgl::pos3(2.0, 2.0, 0.8), 0.8, new kgl::metal(kgl::rgb(0.8, 0.8, 0.8)));
    list[4] = new kgl::sphere(kgl::pos3(2.0, 2.0, 3.0), 0.8, new kgl::metal(kgl::rgb(0.9, 0.9, 0.7), 0.7));
    list[5] = new kgl::sphere(kgl::pos3(3.5, 2.5, 0.2), 0.2, new kgl::lambertian(kgl::rgb(0.9, 0.8, 0.5)));
    list[6] = new kgl::sphere(kgl::pos3(3.5, 1.5, 0.2), 0.2, new kgl::lambertian(kgl::rgb(0.4, 0.8, 0.8)));
    list[7] = new kgl::sphere(kgl::pos3(-6.0, -3.0, 4.0), 4.0, new kgl::metal(kgl::rgb(0.8, 0.8, 0.8)));
    list[8] = new kgl::sphere(kgl::pos3(9, 3, 3.5), 1.0, new kgl::dielectric(2.5));
    list[9] = new kgl::sphere(kgl::pos3(40, -15, 20), 10.0, new kgl::lambertian(kgl::rgb(0.9, 0.9, 0.1)));
    kgl::hittable *world = new kgl::hittable_list(list, num_hittables);


    //kgl::pos3 main_ray;       // Direction of ray
    //kgl::pos3 to_obj;         // From camera position to center of sphere
    kgl::ray main_ray;
    main_ray.origin = cam.pos;
    for (int r = 0; r < height; r++) {
        for (int c = 0; c < width; c++) {
            kgl::rgb avg_color(0, 0, 0);
            int n;
            for (n = 0; n < num_alias_rays; n ++) {
                // Shoot a ray from the camera
                if (num_alias_rays > 1) 
                    main_ray.direction = (top_left - left * (c + drand48() - 0.5) * dc - cam.up * (r + drand48() - 0.5) * dr).normalized();
                else
                    main_ray.direction = (top_left - left * c * dc - cam.up * r * dr).normalized();
                
                avg_color += kgl::get_color(main_ray, world, 0);
            }
            avg_color /= n;
            raster[r][c] = avg_color;
        }
    }

    // Save image
    std::cout << "Saving image..." << std::endl;
    return save_image<width, height>(raster);
}

kgl::rgb kgl::get_color(const kgl::ray &r, kgl::hittable *world, int depth) {
    // Stores information about the intersection
    kgl::hit_record rec;

    // Set the RGB color
    kgl::rgb color{0, 0, 0};

    // Check for hit
    if (world->hit(r, 0.0001, 1000000, rec)) {
        kgl::ray scattered;
        kgl::rgb attenuation;
        if (depth < 5 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
            // TODO Fix this bug. The sky color limits the color visible in the world. Cyan sky meant no red.
            color = attenuation * kgl::get_color(scattered, world, depth + 1);             
        } else {
            color = kgl::rgb(0, 0, 0);
        }
        
        /*
        // This is a value that indicates how much of the normal is pointing in the direction of the light.
        // The more that's in the lights direction, the brighter that spot on the sphere.
        kgl::color_t diffuse = 1;
        diffuse = (rec.normal.dot((light - rec.point).normalized()));
        
        // We need "diffuse" to be from 0 to 1, so that it is a valid color value
        if (diffuse < 0)
            diffuse = 0;
        else if (diffuse > 1)
            diffuse = 1;
        */
       
        //kgl::ray to_light = ray(rec.point, light.normalized());
        //kgl::hit_record shadow_rec;

        /*
        if ((rec.point - kgl::pos3(2, 0, 1)).length() <= 1){
            //color = rgb(0, 1, 0);
            //std::cout << "Attenuation: "  << printv(attenuation)  << std::endl;
            //std::cout << "Color: "  << printv(color)  << std::endl;
        }
        */

        /*
        if (rec.point.z > 0.0){
            color = rgb(0, 1, 0);
        }
        */

        /*
        // Cast shadow rays
        if (world->hit(to_light, 0.0001, to_light.length(), shadow_rec)) {
            diffuse = 0;
            //std::cout << to_light.length() << std::endl;
            //std::cout << shadow_rec.t << std::endl;
            //std::cout << kgl::printv(shadow_rec.point) << std::endl;
        }
        */
        // Multiply color by the diffuse value to get the correct diffuse shading
        //color *= diffuse;

        // Cast another ray for ambient lighting
        /*
        kgl::pos3 target       = rec.normal + kgl::random_in_unit_sphere();
        kgl::rgb ambient_color = kgl::get_color(ray(rec.point, target), world, depth + 1);
        if (ambient_color != background_color) {
            color                  = 0.4 * ambient_color + 0.6 * color;
        }
        */

        // Gamma correction?
        //color = kgl::rgb(sqrt(color.x), sqrt(color.y), sqrt(color.z));

        return color;
    } else {
        kgl::color_t t = 0.5*(r.direction.normalized().z + 1.0);
        return (1.0-t)*kgl::rgb(1.0, 1.0, 1.0) + t*kgl::rgb(0.5, 0.7, 1.0);
        //return background_color;
    }

    
}

//kgl::space_t kgl::hit_sphere(const kgl::pos3 &center, kgl::space_t radius, const kgl::ray &r) {
bool kgl::sphere::hit(const kgl::ray &r, kgl::space_t t_min, kgl::space_t t_max, kgl::hit_record &rec) const {
    // Check for intersections
    kgl::pos3 to_obj = (center - r.origin);

    // Amount of "to_obj" in direction or "ray"
    // This is the parameter in the vector equation: ray = ray_orig + ray_dir * tca
    // It means how far down the ray you have to walk
    kgl::space_t tca = r.direction.dot(to_obj);

    
    // If less than 0, ray is behind us
    if (tca < 0) {
        //return -1.0;
        return false;
    }

    // d^2 + tca^2 = ||to_obj||^2
    // d2 is the closest approach of the ray to the obj
    // We don't solve for d because that would mean expensive square root operations that aren't necessary
    // (since we can just compare to the radius squared)
    kgl::space_t d2 = to_obj.norm2() - tca * tca;

    // If the ray never gets closer than the radius, we know there is no intersection
    if (d2 > radius*radius) {
        //return -1.0;
        return false;
    }
    
    // The ray intersects the sphere, going through two points (think entrance and exit wound)
    // tca is how far from the ray orign to the closest approach of the center.
    // Add thc to go to the exit, or subtract thc to go to the entrance.
    kgl::space_t thc = sqrt(radius*radius - d2),
                 t0  = tca - thc;//,           // Closest intersection
                //t1  = tca + thc;           // Farthest intersection

    // Just grab the closest one
    kgl::space_t t = t0;

    if (t < t_min || t > t_max) {
        return false;
    }

    rec.t       = t;
    rec.point   = r.point_at(rec.t);
    rec.normal  = (rec.point - center).normalized();
    rec.mat_ptr = mat_ptr;

    /*
    if (r.origin.length() < 0.01) {
        std::cout << "Lights: " << std::endl;
        std::cout << "ray_d: "  << printv(r.direction) << std::endl;
        std::cout << "to_obj: " << printv(to_obj) << std::endl;
        std::cout << "point: "  << printv(rec.point) << std::endl;
        std::cout << "tca: "    << tca << std::endl;
        std::cout << "thc: "    << thc << std::endl;
        std::cout << "t: "      << t   << std::endl;
    }
    */
    
    //return t;
    return true;
}

kgl::pos3 kgl::random_in_unit_sphere() {
    kgl::pos3 p;
    do {
        p = 2.0*kgl::pos3(drand48(), drand48(), drand48()) - kgl::pos3(1, 1, 1);
    } while (p.length2() >= 1.0);
    return p;
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
