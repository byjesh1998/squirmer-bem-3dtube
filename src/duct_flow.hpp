// duct_flow.hpp -- Stage 1: pressure-driven Stokes flow in a finite rigid tube
//
// Domain: tube 0 < z < L, radius R. Only a pressure difference is imposed.
//   wall   : u = 0                                   (no slip)
//   inlet  : u_x = u_y = 0,  f_z = sigma.m|_z = +p_in   (m = -e_z)
//   outlet : u_x = u_y = 0,  f_z = -p_out               (m = +e_z)
// Unknowns: all wall tractions; on the caps u_z, f_x, f_y.
// Exact answer (Hagen-Poiseuille):
//   u_z(r) = G (R^2 - r^2)/(4 mu),  G = (p_in - p_out)/L,  Q = pi R^4 G/(8 mu)
//
// The system matrix depends only on the geometry, so it is LU-factorised once;
// every pressure drop is then just a new right-hand side (Stokes flow is linear).
#pragma once
#include <vector>
#include "utils/field.hpp"
#include "utils/kernels.hpp"
#include "utils/linalg.hpp"
#include "utils/mesh.hpp"

namespace bem {

class DuctFlow {
public:
    Mesh mesh;
    double R, L, mu;

    struct Solution {
        std::vector<V3> u, f;
        double Q = 0;  // volume flow rate through the outlet
    };

    DuctFlow(double R_, double L_, int n_theta, int n_z, int n_r, double mu_ = 1.0)
        : R(R_), L(L_), mu(mu_) {
        mesh = tube_mesh(R, uniform_z(0, L, n_z), n_theta, n_r);
        assemble();
    }

    // true if the velocity component (panel p, component b) is prescribed
    bool u_known(int p, int b) const { return !(mesh.tags[p] != WALL && b == 2); }

    Solution solve(double p_in, double p_out) const {
        int n = 3 * mesh.P();
        Mat rhs(n, 1), fk(Fk_.m, 1);
        for (size_t k = 0; k < cap_.size(); ++k)
            fk(k, 0) = mesh.tags[cap_[k]] == INLET ? p_in : -p_out;
        gemm(-1, Fk_, fk, 0, rhs);  // known u are all zero -> no u contribution
        lu_.solve(rhs);
        Solution s;
        s.u.assign(mesh.P(), V3{});
        s.f.assign(mesh.P(), V3{});
        for (int p = 0; p < mesh.P(); ++p)
            for (int b = 0; b < 3; ++b) {
                double x = rhs(3 * p + b, 0);
                if (u_known(p, b)) s.f[p][b] = x;
                else s.u[p][b] = x;
            }
        for (size_t k = 0; k < cap_.size(); ++k) {
            int p = cap_[k];
            s.f[p].z = fk(k, 0);
            if (mesh.tags[p] == OUTLET) s.Q += s.u[p].z * mesh.area[p];
        }
        return s;
    }

    std::vector<V3> velocity(const Solution& s, const std::vector<V3>& pts) const {
        return velocity_at(mesh, s.u, s.f, mu, pts);
    }

    // exact Hagen-Poiseuille solution
    double exact_uz(double r, double dP) const { return dP / L * (R * R - r * r) / (4 * mu); }
    double exact_Q(double dP) const { return PI * std::pow(R, 4) * dP / L / (8 * mu); }

private:
    LUFactor lu_;
    Mat Fk_;                 // columns of A_f belonging to the known cap tractions f_z
    std::vector<int> cap_;   // cap panel indices (order of the columns of Fk_)

    // Completed boundary integral equation at every panel centroid:
    //   A_f f + A_u u = 0,  A_f = SL/(8 pi mu),  A_u = -[DL - diag(sum DL)]/(8 pi)
    // Each unknown scalar takes its column from A_f (traction unknown) or A_u
    // (velocity unknown).
    void assemble() {
        int P = mesh.P(), n = 3 * P;
        for (int p = 0; p < P; ++p)
            if (mesh.tags[p] != WALL) cap_.push_back(p);
        std::vector<int> col_of_cap(P, -1);
        for (size_t k = 0; k < cap_.size(); ++k) col_of_cap[cap_[k]] = (int)k;
        Mat M(n, n);
        Fk_ = Mat(n, (int)cap_.size());
        double SL[9], DL[9];
        for (int i = 0; i < P; ++i) {
            double DLsum[9] = {0};
            for (int p = 0; p < P; ++p) {
                panel_integrals(mesh, p, mesh.centroid[i], i == p, SL, DL);
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b) {
                        double af = SL[3 * a + b] / (EIGHT_PI * mu), au = -DL[3 * a + b] / EIGHT_PI;
                        M(3 * i + a, 3 * p + b) = u_known(p, b) ? af : au;
                        if (b == 2 && col_of_cap[p] >= 0) Fk_(3 * i + a, col_of_cap[p]) = af;
                    }
                for (int k = 0; k < 9; ++k) DLsum[k] += DL[k];
            }
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    if (!u_known(i, b)) M(3 * i + a, 3 * i + b) += DLsum[3 * a + b] / EIGHT_PI;
        }
        lu_.factor(std::move(M));
    }
};

}  // namespace bem
