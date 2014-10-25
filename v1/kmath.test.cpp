#include <iostream>
#include <string>
#include <sstream>
#include "kmath.hpp"

std::string printv (kgl::vec3<double> v) {
    std::ostringstream ss;
    ss << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return ss.str();
}

int main() {
    kgl::vec3<double> dbl;
    std::cout << "dbl: " << printv(dbl) << std::endl;
    dbl.set(1, 2, 3);
    std::cout << "dbl: " << printv(dbl) << std::endl;

    kgl::vec3<double> dbl2(1.0, 2.2, 100.123456800987654321);

    std::cout << "dbl2: " << printv(dbl2) << std::endl;

    kgl::vec3<double> dbl3(3, 5, 6);
    std::cout << "dbl3: " << printv(dbl3) << std::endl;

    double dot = dbl3.dot(dbl);
    std::cout << "dot = dbl3.dbl: " << dot << std::endl;

    dot = dbl.dot(dbl3);
    std::cout << "dot = dbl.dbl3: " << dot << std::endl;


    double norm = dbl.norm();
    std::cout << "norm: " << norm << std::endl;


    kgl::vec3<double> x{1, 0, 0};
    kgl::vec3<double> y(0, 1, 0);
    kgl::vec3<double> z = x.cross(y);
    std::cout << "z: " << printv(z) << std::endl;


    kgl::vec3<double> add = x + y + z;
    std::cout << "add: " << printv(add) << std::endl;

    kgl::vec3<double> sub = x - y + z;
    std::cout << "sub: " << printv(sub) << std::endl;

    kgl::vec3<double> mult = x * 11;
    std::cout << "mult: " << printv(mult) << std::endl;

    kgl::vec3<double> mult2 = 11 * y;
    std::cout << "mult2: " << printv(mult2) << std::endl;

    kgl::vec3<double> div = x / 11;
    std::cout << "div: " << printv(div) << std::endl;

    kgl::vec3<double> div2 = 11 / y;
    std::cout << "div2: " << printv(div2) << std::endl;

    kgl::vec3<double> eq = y;
    std::cout << "eq: " << printv(eq) << std::endl;

    kgl::vec3<double> pleq = x;
    pleq += y;
    std::cout << "pleq: " << printv(pleq) << std::endl;

    return 0;
}
