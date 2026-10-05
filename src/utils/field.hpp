// field.hpp -- velocity at points inside the fluid, from the boundary solution
//
//   u_j(x0) = 1/(8 pi mu) sum_p SL_p,ij f_p,i  -  1/(8 pi) sum_p DL_p,ij u_p,i
//
// (boundary integral representation, normals out of the fluid, f = sigma.m).
// Valid for any x0 strictly inside the fluid; accuracy degrades within about
// one panel size of a boundary, where adaptive quadrature is only partly
// able to resolve the nearly singular integrals.


//
#pragma once
#include <vector>
#include "kernels.hpp"

namespace bem {

inline std::vector<V3> velocity_at(const Mesh& M, const std::vector<V3>& u,
                                   const std::vector<V3>& f, double mu,
                                   const std::vector<V3>& points) {
    std::vector<V3> out(points.size());
    double SL[9], DL[9];
    for (size_t n = 0; n < points.size(); ++n) {
        V3 s{};
        for (int p = 0; p < M.P(); ++p) {
            panel_integrals(M, p, points[n], false, SL, DL);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    s[j] += SL[3 * i + j] * f[p][i] / (EIGHT_PI * mu) - DL[3 * i + j] * u[p][i] / EIGHT_PI;
        }
        out[n] = s;
    }
    return out;
}

}  // namespace bem
