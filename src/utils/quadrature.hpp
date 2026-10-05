// quadrature.hpp -- triangle rules (Dunavant 7-point and its subdivisions),
//
// weights and points are declared
// Gauss-Legendre on [0,1]
#pragma once
#include <array>
#include <cmath>
#include <vector>
#include "vec3.hpp"

namespace bem {

struct Rule {
    std::vector<std::array<double, 3>> L;  // barycentric coordinates
    std::vector<double> W;                 // weights, summing to 1
};

//============================//    
// The points and weights    //
//===========================// 
// 7-point, degree-5 rule
inline Rule dunavant7() {
    const double a1 = 0.059715871789770, b1 = 0.470142064105115, w1 = 0.132394152788506;
    const double a2 = 0.797426985353087, b2 = 0.101286507323456, w2 = 0.125939180544827;
    Rule r;
    r.L = {{1. / 3, 1. / 3, 1. / 3}, {a1, b1, b1}, {b1, a1, b1}, {b1, b1, a1}, //coordinates and its three rotations
           {a2, b2, b2}, {b2, a2, b2}, {b2, b2, a2}};
    r.W = {0.225, w1, w1, w1, w2, w2, w2};
    return r;
}


//=========================================//    
// Subdividing into 4 sub-triangle         //
//=========================================// 
// Dunavant-7 applied on the 4^level sub-triangles of the reference triangle
inline Rule subdivided(int level) {
    using T3 = std::array<std::array<double, 3>, 3>;
    std::vector<T3> tris = {T3{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}};
    auto mid = [](const std::array<double, 3>& p, const std::array<double, 3>& q) {
        return std::array<double, 3>{(p[0] + q[0]) / 2, (p[1] + q[1]) / 2, (p[2] + q[2]) / 2};
    };
    for (int l = 0; l < level; ++l) { //depending on level each big triangle is divided into 4
        std::vector<T3> nw;
        for (auto& t : tris) {
            auto ab = mid(t[0], t[1]), bc = mid(t[1], t[2]), ca = mid(t[2], t[0]); // middle points of a triangle sides are calculated
            nw.push_back({t[0], ab, ca});
            nw.push_back({ab, t[1], bc}); //triangle split into four, 3 on each corner and one on center
            nw.push_back({ca, bc, t[2]});
            nw.push_back({ab, bc, ca});
        }
        tris = nw;
    }
    Rule d = dunavant7(), r;
    for (auto& t : tris)
        for (size_t q = 0; q < d.W.size(); ++q) {
            std::array<double, 3> l{0, 0, 0};
            for (int k = 0; k < 3; ++k)
                for (int c = 0; c < 3; ++c) l[c] += d.L[q][k] * t[k][c];
            r.L.push_back(l);
            r.W.push_back(d.W[q] / tris.size()); // weight renormalizing according to number of subdivided triangles
        }
    return r;
}



// rules()[0] = 7 points, [1] = 28, [2] = 112, [3] = 448
inline const std::vector<Rule>& rules() {
    static const std::vector<Rule> R = {dunavant7(), subdivided(1), subdivided(2), subdivided(3)};
    return R;
}


//=========================================//    
// Gauss-Legendre techniques               //
//=========================================// 
// n-point Gauss-Legendre nodes/weights mapped to [0, 1]
inline void gauss_legendre01(int n, std::vector<double>& x, std::vector<double>& w) {
    x.resize(n);
    w.resize(n);
    for (int i = 0; i < n; ++i) {
        double z = std::cos(PI * (i + 0.75) / (n + 0.5)), pp = 0; //guess value: it begins with a guess based on a cosine formula.
        for (int it = 0; it < 100; ++it) { // evaluating legendre polunomial
            double p1 = 1, p2 = 0;
            for (int j = 1; j <= n; ++j) {
                double p3 = p2;
                p2 = p1;  // degree n-1 Legendre p.
                p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;// degree n Legendre p.
            }
            pp = n * (z * p1 - p2) / (z * z - 1); //computes the polynomial's derivative
            double dz = p1 / pp;
            z -= dz;
            if (std::fabs(dz) < 1e-15) break;
        }
        x[i] = 0.5 * (1 - z); // move from [−1, 1] to [0, 1]
        w[i] = 1.0 / ((1 - z * z) * pp * pp); // compute the weights
    }
}

}  // namespace bem
