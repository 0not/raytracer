#ifndef KTYPES_H
#define KTYPES_H
#include "kmath.hpp"

namespace kgl {
    typedef float color_t;
    color_t RGB_MAX = 1.0;
    unsigned short RGB_BITS = 8;
    typedef vec3<color_t> rgb;
    
    template <unsigned short W, unsigned short H>
    using raster_t = std::array<std::array<rgb, W>, H>;

    // Some default colors
    rgb red{1.0, 0.0, 0.0};
    rgb yellow{1.0, 1.0, 0.0};
    rgb green{0.0, 1.0, 0.0};
    rgb cyan{0.0, 1.0, 1.0};
    rgb blue{0.0, 0.0, 1.0};
    rgb purple{1.0, 0.0, 1.0};
    rgb white{1.0, 1.0, 1.0};
    rgb black{0.0, 0.0, 0.0};
    
    // Spatial component type
    typedef double space_t; // 1.0 = 1.0m (meters)
    typedef vec3<space_t> pos3;

};

#endif
