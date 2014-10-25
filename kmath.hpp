#ifndef KMATH_H
#define KMATH_H
#include <cmath>

namespace kgl {
    // 3-d vector
    // el is element type
    template <class T> class vec3 {
    public:
        T x, y, z;

        void set(const T &xt, const T &yt, const T &zt) {
            x = xt;
            y = yt;
            z = zt;
        }

        void set(const vec3<T> &v) {
            set(v.x, v.y, v.z);
        }

        vec3() : vec3((T)0, (T)0, (T)0) {}
        
        vec3(T x, T y, T z) {
            set(x, y, z);
        }

        vec3(const vec3<T> &v) {
            set(v.x, v.y, v.z);
        }

        static T dot(const vec3<T> &v1, const vec3<T> v2) {
            return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;    
        }

        T dot(const vec3<T> &v2) {
            return dot(*this, v2);  
        }

        T norm() {
            return sqrt(dot(*this));
        }
       
        // Norm squared
        T norm2() {
            return dot(*this);
        }

        T length()  { return norm(); }
        T length2() { return norm2(); }
        
        vec3<T> normalized() {
            T len = norm();
            if (len == 0) len = 1.0;
            return *this/len;
        }

        vec3<T> normalize() {
            T len = norm();
            if (len > 0)
                this->set(*this/len);
            return *this;
        }

        static vec3<T> cross(const vec3<T> &v1, const vec3<T> &v2) {
            return vec3(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
        }
        
        vec3<T> cross(const vec3<T> &v2) {
            return cross(*this, v2);
        }


        constexpr vec3<T> operator +(const vec3<T> &v2)  {
            return vec3(x + v2.x, y + v2.y, z + v2.z);
        }

        constexpr vec3<T> operator -(const vec3<T> &v2) const {
            return vec3(x - v2.x, y - v2.y, z - v2.z);
        }

        constexpr vec3<T> operator *(const T a) const {
            return vec3(a * x, a * y, a * z);
        }

        constexpr friend vec3<T> operator *(const T a, const vec3<T> &v1) {
            return v1 * a;
        }

        constexpr vec3<T> operator /(const T a) const {
            return vec3(1/a * x, 1/a * y, 1/a * z);
        }

        constexpr friend vec3<T> operator /(const T a, const vec3<T> &v1) {
            return v1 / a;
        }

        constexpr bool operator ==(const vec3<T> &v2) {
            return (x == v2.x && y == v2.y && z == v2.z);
        }

        constexpr bool operator !=(const vec3<T> &v2) {
            return !(x == v2.x && y == v2.y && z == v2.z);
        }

        vec3<T> & operator +=(const vec3<T> &v2) {
            x += v2.x;
            y += v2.y;
            z += v2.z;
            return *this;
        }
        
        vec3<T> & operator -=(const vec3<T> &v2) {
            x -= v2.x;
            y -= v2.y;
            z -= v2.z;
            return *this;
        }

        vec3<T> & operator *=(const vec3<T> &v2) {
            x *= v2.x;
            y *= v2.y;
            z *= v2.z;
            return *this;
        }

        vec3<T> & operator /=(const vec3<T> &v2) {
            x /= v2.x;
            y /= v2.y;
            z /= v2.z;
            return *this;
        }
    };
};

#endif
