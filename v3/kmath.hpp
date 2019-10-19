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

        constexpr vec3<T>& operator =(const vec3<T> &v) = default;

        T dot(const vec3<T> &v1, const vec3<T> v2)  const {
            return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;    
        }

        T dot(const vec3<T> &v2) const {
            return dot(*this, v2);  
        }

        T norm() const { return sqrt(dot(*this)); }
       
        // Norm squared
        T norm2() const { return dot(*this); }

        T length()  const { return norm();  }
        T length2() const { return norm2(); }
        
        vec3<T> normalized() const {
            T len = norm();
            if (len == 0) len = 1;
            return *this/len;
        }

        void normalize() {
            T len = norm();
            if (len > 0)
                this->set(*this/len);
        }

        static vec3<T> cross(const vec3<T> &v1, const vec3<T> &v2) {
            return vec3(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
        }
        
        vec3<T> cross(const vec3<T> &v2) {
            return cross(*this, v2);
        }


        constexpr vec3<T> operator +(const vec3<T> &v2) const {
            return vec3(x + v2.x, y + v2.y, z + v2.z);
        }

        constexpr vec3<T> operator -(const vec3<T> &v2) const {
            return vec3(x - v2.x, y - v2.y, z - v2.z);
        }

        constexpr vec3<T> operator -() const {
            return -1.0 * *this; 
        }

        constexpr vec3<T> operator *(const T a) const {
            return vec3(a * x, a * y, a * z);
        }

        constexpr friend vec3<T> operator *(const T a, const vec3<T> &v1) {
            return v1 * a;
        }

        constexpr vec3<T> operator *(const vec3<T> &v2) const {
            return vec3(x*v2.x, y*v2.y, z*v2.z);
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

        vec3<T> & operator *=(const T &a) {
            x *= a;
            y *= a;
            z *= a;
            return *this;
        }

        vec3<T> & operator *=(const vec3<T> &v2) {
            *this = *this * v2;
            /*
            x *= v2.x;
            y *= v2.y;
            z *= v2.z;
            */
            return *this;
        }

        vec3<T> & operator /=(const T &a) {
            x /= a;
            y /= a;
            z /= a;
            return *this;
        }
    };

    template <class T>
    vec3<T> reflect(const vec3<T> &v, const vec3<T> &n) {
        return v - 2*v.dot(n)*n;
    }

    template <class T>
    bool refract(const vec3<T> &v, const vec3<T> &n, float ni_over_nt, vec3<T> &refracted) {
        vec3<T> uv = v.normalized();
        T dt = uv.dot(n);
        T discriminant = 1.0 - ni_over_nt*ni_over_nt*(1 - dt*dt);
        if (discriminant > 0) {
            refracted = ni_over_nt*(uv - n*dt) - n*sqrt(discriminant);
            return true;
        } else
        {
            return false;
        }
        
    }

    float schlick(float cosine, float ref_idx) {
        float r0 = (1 - ref_idx) / (1 + ref_idx);
        r0 = r0 * r0;
        return r0 + (1 - r0)*pow((1 - cosine), 5);
    }
};

// From:
//    https://gist.github.com/mortennobel/8665258

#define RAND48_SEED_0 (0x330e)
#define RAND48_SEED_1 (0xabcd)
#define RAND48_SEED_2 (0x1234)
#define RAND48_MULT_0 (0xe66d)
#define RAND48_MULT_1 (0xdeec)
#define RAND48_MULT_2 (0x0005)
#define RAND48_ADD (0x000b)

unsigned short _rand48_seed[3] = {
        RAND48_SEED_0,
         RAND48_SEED_1,
         RAND48_SEED_2
};
unsigned short _rand48_mult[3] = {
         RAND48_MULT_0,
         RAND48_MULT_1,
         RAND48_MULT_2
 };
unsigned short _rand48_add = RAND48_ADD;

void
 _dorand48(unsigned short xseed[3])
 {
	         unsigned long accu;
	         unsigned short temp[2];
	
	         accu = (unsigned long)_rand48_mult[0] * (unsigned long)xseed[0] +
	          (unsigned long)_rand48_add;
	         temp[0] = (unsigned short)accu;        /* lower 16 bits */
	         accu >>= sizeof(unsigned short)* 8;
	         accu += (unsigned long)_rand48_mult[0] * (unsigned long)xseed[1] +
	          (unsigned long)_rand48_mult[1] * (unsigned long)xseed[0];
	         temp[1] = (unsigned short)accu;        /* middle 16 bits */
	         accu >>= sizeof(unsigned short)* 8;
	         accu += _rand48_mult[0] * xseed[2] + _rand48_mult[1] * xseed[1] + _rand48_mult[2] * xseed[0];
	         xseed[0] = temp[0];
	         xseed[1] = temp[1];
	         xseed[2] = (unsigned short)accu;
}

double erand48(unsigned short xseed[3])
{
         _dorand48(xseed);
         return ldexp((double) xseed[0], -48) +
                ldexp((double) xseed[1], -32) +
                ldexp((double) xseed[2], -16);
}

double drand48(){
	return erand48(_rand48_seed);
}

void srand48(long seed){
	_rand48_seed[0] = RAND48_SEED_0;
	_rand48_seed[1] = (unsigned short)seed;
	_rand48_seed[2] = (unsigned short)(seed >> 16);
	_rand48_mult[0] = RAND48_MULT_0;
	_rand48_mult[1] = RAND48_MULT_1;
	_rand48_mult[2] = RAND48_MULT_2;
	_rand48_add = RAND48_ADD;
}

#endif
