// kernels.hpp -- free-space Green's functions of Stokes flow integrated over panels
//
//   Stokeslet   G_ij(r)  = d_ij / |r| + r_i r_j / |r|^3
//   Stresslet   T_ijk(r) = -6 r_i r_j r_k / |r|^5          r = x - x0
//
// panel_integrals() returns, for panel p and target point x0,
//   SL_ij = \int_p G_ij dS          DL_ij = \int_p T_ijk m_k dS
// (row-major 3x3 arrays). Both are symmetric in (i, j).
//
// Quadrature choice by d = |centroid_p - x0| / h_p :
//   d >= 3     : 7 points          1.5 <= d < 3 : 28 points
//   0.75<=d<1.5: 112 points        d < 0.75     : 448 points
//   x0 inside p (self panel): Duffy transform for the 1/r singularity of G;
//   the stresslet term vanishes on a flat panel (r . m = 0).
#pragma once
#include <vector>
#include "mesh.hpp"
#include "quadrature.hpp"
#include "vec3.hpp"



namespace bem {

//=============================================//    
// defining stokeslet and stresslet kernels    //
//=============================================// 
// accumulate w*G(r) into SL and w*T(r).m into DL
inline void accum_kernels(const V3& r, double w, const V3& m, double* SL, double* DL) {
    double r2 = dot(r, r), rn = std::sqrt(r2), ir = 1 / rn, ir3 = ir / r2;
    double c[3] = {r.x, r.y, r.z};
    double d = -6 * w * dot(r, m) * ir3 / r2, s = w * ir3;
    for (int i = 0; i < 3; ++i) {
        SL[4 * i] += w * ir;
        for (int j = 0; j < 3; ++j) {
            SL[3 * i + j] += s * c[i] * c[j];
            DL[3 * i + j] += d * c[i] * c[j];
        }
    }
}


//===========================//    
// doing the integrations    //
//===========================// 
// \int_T G dS over a flat triangle that contains the singular point x0.
// The triangle is split into three sub-triangles (x0, v_a, v_b); each is mapped
// from the unit square by the Duffy transform x = x0 + u[(1-v) e1 + v e2], whose
// Jacobian (~u) cancels the 1/r singularity.
inline void self_single_layer(const std::array<V3, 3>& V, const V3& x0, double* SL) {
    static std::vector<double> g, gw;
    if (g.empty()) gauss_legendre01(10, g, gw);
    for (int k = 0; k < 9; ++k) SL[k] = 0;
    double dummy[9];
    for (int e = 0; e < 3; ++e) {
        V3 e1 = V[e] - x0, e2 = V[(e + 1) % 3] - x0;
        double J = norm(cross(e1, e2));
        for (size_t a = 0; a < g.size(); ++a)
            for (size_t b = 0; b < g.size(); ++b) {
                double u = g[a], v = g[b];
                V3 r = u * ((1 - v) * e1 + v * e2);
                accum_kernels(r, gw[a] * gw[b] * u * J, V3{0, 0, 0}, SL, dummy);
            }
    }
}

inline void panel_integrals(const Mesh& M, int p, const V3& x0, bool self, double* SL,
                            double* DL) {
    for (int k = 0; k < 9; ++k) SL[k] = DL[k] = 0;
    if (self) {
        self_single_layer(M.V[p], x0, SL);
        return;
    }
    double d = norm(M.centroid[p] - x0) / M.h[p];
    int lev = d < 0.75 ? 3 : (d < 1.5 ? 2 : (d < 3.0 ? 1 : 0));
    const Rule& R = rules()[lev];
    const auto& Vp = M.V[p];
    for (size_t q = 0; q < R.W.size(); ++q) {
        V3 x = R.L[q][0] * Vp[0] + R.L[q][1] * Vp[1] + R.L[q][2] * Vp[2];
        accum_kernels(x - x0, R.W[q] * M.area[p], M.normal[p], SL, DL);
    }
}

// far panel
//    ↓
//   7 points

// nearer
//    ↓
//  28 points

// close
//    ↓
// 112 points

// very close
//    ↓
// 448 points


//===========================//    
// squirmer surface slip     //
//===========================// 
// Tangential squirming velocity of the Blake / Ishikawa-Pedley squirmer:
//   u_s = (B1 + B2 cos th) sin th e_th,   cos th = rhat . e
// using  sin th e_th = (cos th) rhat - e.
inline V3 squirmer_slip(V3 x, V3 center, V3 e, double B1, double B2) {
    V3 rh = unit(x - center);
    double ct = dot(rh, e);
    return (B1 + B2 * ct) * (ct * rh - e);
}

}  // namespace bem
//
    //          panel p
    //             │
    //             ▼
    //    distance / panel size
    //             │
    //             ▼
    //     ┌───────────────┐
    //     │ choose rule   │
    //     └───────┬───────┘
    //             │
    //    ┌────────┼─────────┐
    //    │        │         │
    //  far      near      self
    //    │        │         │
    //   7        28/112     Duffy
    // points     /448        │
    //    │        │          │
    //    └────────┼──────────┘
    //             ▼
    //       evaluate G,T
    //             │
    //             ▼
    //       integrate over
    //          triangle
    //             │
    //             ▼
    //         SL, DL