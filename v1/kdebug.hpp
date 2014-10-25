#ifndef KDEBUG_H
#define KDEBUG_H
#include "kmath.hpp"

namespace kgl {
    std::string printv (kgl::vec3<kgl::color_t> v) {
        std::ostringstream ss;
        ss << "(" << v.x << ", " << v.y << ", " << v.z << ")";
        return ss.str();
    }
    std::string printv (kgl::pos3 v) {
        std::ostringstream ss;
        ss << "(" << v.x << ", " << v.y << ", " << v.z << ")";
        return ss.str();
    }
};
#endif
