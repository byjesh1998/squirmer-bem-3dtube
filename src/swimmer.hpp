// swimmer.hpp -- Stage 2: a spherical squirmer in a closed rigid tube (or in
// unbounded fluid)
//
// Boundary conditions
//   tube wall and caps : u = 0      (closed tube, as in Zhu, Lauga & Brandt 2013)
//   sphere surface     : u = U + Omega x (x - xc) + u_s(x)
//   unknowns           : traction everywhere, U, Omega
//   closure            : force-free  \int f dS = 0,  torque-free \int (x-xc) x f dS = 0
//

// Linear system (t = tube, s = sphere), with T the fixed tube-tube block:
//   [ T  B  Ct ] [f_t ]   [r_t]
//   [ D  E  Cs ] [f_s ] = [r_s]
//   [ 0  F  0  ] [UOm ]   [ 0 ]
// T is LU-factorised once (class Tube); the swimmer is eliminated through the
// Schur complement  G = [E | Cs] - D T^{-1} [B | Ct], of size 3*Ns + 6.
//
// Null modes: on each closed surface with velocity prescribed everywhere, a
// uniform normal traction f = c m produces no flow (incompressibility), so the
// single-layer block is singular. It is removed by Wielandt deflation (adding
// a rank-one term); this changes no velocity, force, torque or power.
#pragma once
#include <vector>
#include "utils/field.hpp"
#include "utils/io.hpp"
#include "utils/kernels.hpp"
#include "utils/linalg.hpp"
#include "utils/mesh.hpp"

namespace bem {

// ---------------------------------------------------------------- tube
struct Tube {
    Mesh mesh;
    int Nt = 0;
    double R = 1, L = 4, mu = 1;
    LUFactor lu;

    Tube() = default;
    Tube(double R_, double L_, double a, int n_theta, int n_r = 4, double mu_ = 1)
        : R(R_), L(L_), mu(mu_) {
        mesh = swimmer_tube_mesh(R, L, a, n_theta, n_r);
        Nt = mesh.P();
        factor();
    }

private:
    void factor() {
        int n = 3 * Nt;
        Mat T(n, n);
        double SL[9], DL[9];
        for (int i = 0; i < Nt; ++i)
            for (int p = 0; p < Nt; ++p) {
                panel_integrals(mesh, p, mesh.centroid[i], i == p, SL, DL);
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b) T(3 * i + a, 3 * p + b) = SL[3 * a + b] / (EIGHT_PI * mu);
            }
        deflate(mesh, 0, Nt, T, 0, mu);
        lu.factor(std::move(T));
    }

public:
    // rank-one deflation of the null mode f = c m on panels [p0, p0+np) of M
    static void deflate(const Mesh& M, int p0, int np, Mat& A, int off, double mu) {
        double At = 0;
        for (int p = p0; p < p0 + np; ++p) At += M.area[p];
        double sc = 1.0 / (EIGHT_PI * mu * std::sqrt(At));
        for (int i = 0; i < np; ++i)
            for (int p = 0; p < np; ++p)
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b)
                        A(off + 3 * i + a, off + 3 * p + b) +=
                            sc * M.normal[p0 + i][a] * M.normal[p0 + p][b] * M.area[p0 + p];
    }
};

// ---------------------------------------------------------------- swimmer
struct SwimmerSolution {
    V3 U, Om;
    std::vector<V3> f_sphere, f_tube;  // tractions (f = sigma.m, m out of the fluid)
};

class SwimmerSolver {
public:
    const Tube* tube;  // nullptr -> unbounded fluid
    Mesh S;            // sphere mesh
    V3 xc;
    double mu;
    int Ns, Nt;
    double t_assembly = 0, t_schur = 0;

    SwimmerSolver(const Tube* t, const Mesh& sphere, V3 center, double mu_ = 1.0)
        : tube(t), S(sphere), xc(center), mu(mu_) {
        Ns = S.P();
        Nt = tube ? tube->Nt : 0;
        build();
    }

    // Squirming slip on every sphere panel for orientation e and modes B1, B2
    std::vector<V3> slip(V3 e, double B1, double B2) const {
        std::vector<V3> u(Ns);
        for (int p = 0; p < Ns; ++p) u[p] = squirmer_slip(S.centroid[p], xc, e, B1, B2);
        return u;
    }

    // Force- and torque-free swimming, one solution per slip field
    std::vector<SwimmerSolution> solve_free(const std::vector<std::vector<V3>>& slips) {
        int ns = 3 * Ns, n = ns + 6, nr = (int)slips.size();
        Mat M(n, n), Rh(n, nr);
        for (int j = 0; j < ns; ++j)
            for (int i = 0; i < ns; ++i) M(i, j) = Gs_(i, j);
        for (int j = 0; j < 6; ++j)
            for (int i = 0; i < ns; ++i) M(i, ns + j) = Gc_(i, j);
        double amax = 0, Amax = 0;
        for (double v : Gs_.a) amax = std::max(amax, std::fabs(v));
        for (double A : S.area) Amax = std::max(Amax, 5.0 * A);
        double sc = amax / Amax;  // row scaling of the constraint equations
        for (int p = 0; p < Ns; ++p) {
            V3 r = S.centroid[p] - xc;
            double A = S.area[p] * sc;
            for (int a = 0; a < 3; ++a) M(ns + a, 3 * p + a) = A;           // sum f dA = 0
            M(ns + 3, 3 * p + 1) = -A * r.z; M(ns + 3, 3 * p + 2) = A * r.y;  // sum r x f dA = 0
            M(ns + 4, 3 * p + 0) = A * r.z;  M(ns + 4, 3 * p + 2) = -A * r.x;
            M(ns + 5, 3 * p + 0) = -A * r.y; M(ns + 5, 3 * p + 1) = A * r.x;
        }
        std::vector<Mat> rt(nr);
        for (int k = 0; k < nr; ++k) {
            auto rs = reduced_rhs(slips[k], rt[k]);
            for (int i = 0; i < ns; ++i) Rh(i, k) = rs[i];
        }
        solve_dense(M, Rh);
        std::vector<SwimmerSolution> out(nr);
        for (int k = 0; k < nr; ++k) {
            Mat x(ns + 6, 1);
            for (int i = 0; i < ns + 6; ++i) x(i, 0) = Rh(i, k);
            out[k] = unpack(x, rt[k]);
        }
        return out;
    }

    // Rigid sphere with prescribed (U, Omega) and no slip ("towed" problem)
    SwimmerSolution solve_prescribed(V3 U, V3 Om) {
        int ns = 3 * Ns;
        Mat w(6, 1), rhs(ns, 1);
        for (int a = 0; a < 3; ++a) { w(a, 0) = U[a]; w(3 + a, 0) = Om[a]; }
        gemm(-1, Gc_, w, 0, rhs);
        solve_dense(Gs_, rhs);
        Mat x(ns + 6, 1);
        for (int i = 0; i < ns; ++i) x(i, 0) = rhs(i, 0);
        for (int i = 0; i < 6; ++i) x(ns + i, 0) = w(i, 0);
        Mat rt(3 * Nt, 1);
        return unpack(x, rt);
    }

    // Hydrodynamic force ON the sphere: F = -sum f dA (since m points into it)
    V3 force(const SwimmerSolution& s) const {
        V3 F{};
        for (int p = 0; p < Ns; ++p) F = F - S.area[p] * s.f_sphere[p];
        return F;
    }

    // Velocity at points in the fluid (lab frame)
    std::vector<V3> velocity(const SwimmerSolution& s, const std::vector<V3>& slip_field,
                             const std::vector<V3>& pts) const {
        Mesh all = tube ? merge(tube->mesh, S) : S;
        std::vector<V3> u(all.P()), f(all.P());
        for (int p = 0; p < Nt; ++p) f[p] = s.f_tube[p];
        for (int p = 0; p < Ns; ++p) {
            u[Nt + p] = s.U + cross(s.Om, S.centroid[p] - xc) + slip_field[p];
            f[Nt + p] = s.f_sphere[p];
        }
        return velocity_at(all, u, f, mu, pts);
    }

private:
    Mat Gs_, Gc_;       // Schur operator: sphere tractions, rigid-body columns
    Mat AuBt_, AuBs_;   // A_u restricted to sphere columns (tube rows, sphere rows)
    Mat D_, Y_;         // D (sphere rows <- tube tractions), Y = T^{-1}[B | Ct]

    void build() {
        Timer clock;
        int ns = 3 * Ns, nt = 3 * Nt;
        Mat B(nt, ns), E(ns, ns);
        AuBt_ = Mat(nt, ns);
        AuBs_ = Mat(ns, ns);
        D_ = Mat(ns, nt);
        double SL[9], DL[9];
        for (int i = 0; i < Nt; ++i) {  // tube rows <- sphere sources
            V3 x0 = tube->mesh.centroid[i];
            for (int p = 0; p < Ns; ++p) {
                panel_integrals(S, p, x0, false, SL, DL);
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b) {
                        B(3 * i + a, 3 * p + b) = SL[3 * a + b] / (EIGHT_PI * mu);
                        AuBt_(3 * i + a, 3 * p + b) = -DL[3 * a + b] / EIGHT_PI;
                    }
            }
        }
        for (int i = 0; i < Ns; ++i) {  // sphere rows <- tube and sphere sources
            V3 x0 = S.centroid[i];
            double DLsum[9] = {0};
            for (int p = 0; p < Nt; ++p) {
                panel_integrals(tube->mesh, p, x0, false, SL, DL);
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b) D_(3 * i + a, 3 * p + b) = SL[3 * a + b] / (EIGHT_PI * mu);
                for (int k = 0; k < 9; ++k) DLsum[k] += DL[k];
            }
            for (int p = 0; p < Ns; ++p) {
                panel_integrals(S, p, x0, p == i, SL, DL);
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b) {
                        E(3 * i + a, 3 * p + b) = SL[3 * a + b] / (EIGHT_PI * mu);
                        AuBs_(3 * i + a, 3 * p + b) = -DL[3 * a + b] / EIGHT_PI;
                    }
                for (int k = 0; k < 9; ++k) DLsum[k] += DL[k];
            }
            // completed double layer; exterior (unbounded) problem gets an extra -u(x0)
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    AuBs_(3 * i + a, 3 * i + b) += DLsum[3 * a + b] / EIGHT_PI - (tube ? 0.0 : (a == b));
        }
        Tube::deflate(S, 0, Ns, E, 0, mu);
        // rigid-body map K: u = U + Omega x r = U - [r]x Omega
        Mat K(ns, 6);
        for (int p = 0; p < Ns; ++p) {
            V3 r = S.centroid[p] - xc;
            for (int a = 0; a < 3; ++a) K(3 * p + a, a) = 1;
            K(3 * p + 0, 4) = r.z;  K(3 * p + 0, 5) = -r.y;
            K(3 * p + 1, 3) = -r.z; K(3 * p + 1, 5) = r.x;
            K(3 * p + 2, 3) = r.y;  K(3 * p + 2, 4) = -r.x;
        }
        Mat Ct(nt, 6), Cs(ns, 6);
        gemm(1, AuBt_, K, 0, Ct);
        gemm(1, AuBs_, K, 0, Cs);
        t_assembly = clock.seconds();
        Gs_ = E;
        Gc_ = Cs;
        if (tube) {
            Y_ = Mat(nt, ns + 6);
            std::copy(B.a.begin(), B.a.end(), Y_.a.begin());
            std::copy(Ct.a.begin(), Ct.a.end(), Y_.a.begin() + (size_t)nt * ns);
            tube->lu.solve(Y_);
            Mat DY(ns, ns + 6);
            gemm(1, D_, Y_, 0, DY);
            for (int j = 0; j < ns; ++j)
                for (int i = 0; i < ns; ++i) Gs_(i, j) -= DY(i, j);
            for (int j = 0; j < 6; ++j)
                for (int i = 0; i < ns; ++i) Gc_(i, j) -= DY(i, ns + j);
        }
        t_schur = clock.seconds() - t_assembly;
    }



    // reduced right-hand side r_s - D T^{-1} r_t; also returns T^{-1} r_t
    std::vector<double> reduced_rhs(const std::vector<V3>& us, Mat& rt) {
        int ns = 3 * Ns, nt = 3 * Nt;
        Mat u(ns, 1), rs(ns, 1);
        for (int p = 0; p < Ns; ++p)
            for (int a = 0; a < 3; ++a) u(3 * p + a, 0) = us[p][a];
        gemm(-1, AuBs_, u, 0, rs);
        rt = Mat(nt, 1);
        if (tube) {
            gemm(-1, AuBt_, u, 0, rt);
            tube->lu.solve(rt);
            gemm(-1, D_, rt, 1, rs);
        }
        return rs.a;
    }


    
    // x = [f_s; U; Omega] -> full solution, recovering f_t = T^{-1} r_t - Y x
    SwimmerSolution unpack(const Mat& x, const Mat& rt) const {
        int ns = 3 * Ns;
        SwimmerSolution s;
        s.U = {x(ns, 0), x(ns + 1, 0), x(ns + 2, 0)};
        s.Om = {x(ns + 3, 0), x(ns + 4, 0), x(ns + 5, 0)};
        s.f_sphere.resize(Ns);
        for (int p = 0; p < Ns; ++p) s.f_sphere[p] = {x(3 * p, 0), x(3 * p + 1, 0), x(3 * p + 2, 0)};
        if (tube) {
            Mat ft = rt;
            gemm(-1, Y_, x, 1, ft);
            s.f_tube.resize(Nt);
            for (int p = 0; p < Nt; ++p) s.f_tube[p] = {ft(3 * p, 0), ft(3 * p + 1, 0), ft(3 * p + 2, 0)};
        }
        return s;
    }
};

}  // namespace bem
