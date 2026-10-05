// vec3.hpp -- minimal 3-vector type and constants
#pragma once
#include <cmath>

namespace bem {

constexpr double PI = 3.14159265358979323846;
constexpr double EIGHT_PI = 8.0 * PI;

//=======================//    
// Vector operations.    //
//=======================// 
struct V3 {
    double x = 0, y = 0, z = 0;
    double& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
    double operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
};
inline V3 operator+(V3 a, V3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline V3 operator-(V3 a, V3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline V3 operator*(double s, V3 a) { return {s * a.x, s * a.y, s * a.z}; }
inline double dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline V3 cross(V3 a, V3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(V3 a) { return std::sqrt(dot(a, a)); }
inline V3 unit(V3 a) { return (1.0 / norm(a)) * a; }

}  // namespace bem
