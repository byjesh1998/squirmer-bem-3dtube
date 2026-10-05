// main.cpp -- command-line driver
//
//   squirmer_bem <input-file> [key=value ...]
//
// The input file selects a mode:
//   duct        Stage 1: pressure-driven flow in a tube (no swimmer)
//   kinematics  Stage 2: swimming velocity/rotation (+ towed drag) for lists of
//               confinements a/R and eccentricities beta
//   field       Stage 2: velocity field around a squirmer on a grid in the x-z plane
//   trajectory  Stage 3: time-integrated swimmer path
// See inputs/*.in for documented examples and README.md for all keys.

#include <cstdio>
#include <memory>
#include <string>
#include "duct_flow.hpp"
#include "swimmer.hpp"
#include "trajectory.hpp"
#include "utils/config.hpp"
#include "utils/io.hpp"

using namespace bem;

//============================//    
// duct flow                  //
//===========================// 
static void run_duct(const Config& c, const std::string& out) {
    double R = c.num("R", 1.0), L = c.num("L", 4.0), mu = c.num("mu", 1.0);
    int nth = c.integer("n_theta", 24), nz = c.integer("n_z", 24), nr = c.integer("n_r", 6);
    double p_out = c.num("p_out", 0.0), dP_detail = c.num("dP_detail", 1.0);
    auto dPs = c.list("dP_list", {0.5, 1, 2, 4});
    Timer clock;
    DuctFlow duct(R, L, nth, nz, nr, mu);
    std::printf("duct flow: %d panels, %d unknowns, assembled + factorised in %.1f s\n",
                duct.mesh.P(), 3 * duct.mesh.P(), clock.seconds());
    Csv q(out + "/flow_rate.csv", "dP,Q,Q_exact,rel_err");
    for (double dP : dPs) {
        auto s = duct.solve(p_out + dP, p_out);
        double Qe = duct.exact_Q(dP);
        q.row({dP, s.Q, Qe, std::fabs(s.Q - Qe) / Qe});
        std::printf("  dP = %-5g  Q = %.6f   exact %.6f   rel. error %.2e\n", dP, s.Q, Qe,
                    std::fabs(s.Q - Qe) / Qe);
    }
    // detailed solution for one pressure drop
    auto s = duct.solve(p_out + dP_detail, p_out);
    write_panels(out + "/panels.csv", duct.mesh, s.u, s.f);
    int nprof = c.integer("n_profile", 20);
    std::vector<V3> pts;
    for (double ang : {0.0, PI / 5, PI / 2})
        for (int k = 0; k < nprof; ++k) {
            double r = 0.95 * R * k / (nprof - 1);
            pts.push_back({r * std::cos(ang), r * std::sin(ang), L / 2});
        }
    auto v = duct.velocity(s, pts);
    Csv pr(out + "/profile.csv", "x,y,z,r,ux,uy,uz,uz_exact",
           {"interior velocity at z = L/2, dP = " + std::to_string(dP_detail)});
    for (size_t k = 0; k < pts.size(); ++k) {
        double r = std::hypot(pts[k].x, pts[k].y);
        pr.row({pts[k].x, pts[k].y, pts[k].z, r, v[k].x, v[k].y, v[k].z, duct.exact_uz(r, dP_detail)});
    }
    std::printf("wrote %s/{flow_rate,panels,profile}.csv  (%.1f s)\n", out.c_str(), clock.seconds());
}



//============================//    
// kinematics.               //
//===========================// 
static void run_kinematics(const Config& c, const std::string& out) {
    const double a = 1.0, U0 = 2.0 / 3;
    auto lams = c.list("a_over_R", {0.3});
    auto betas = c.list("beta", {0.0});
    int nth = c.integer("n_theta", 24), lev = c.integer("sphere_level", 3);
    double LoR = c.num("L_over_R", 4.0);
    bool towed = c.integer("towed_drag", 1) != 0;
    Csv k(out + "/kinematics.csv",
          "a_over_R,beta,panels,UB1x,UB1y,UB1z,OB1x,OB1y,OB1z,UB2x,UB2y,UB2z,OB2x,OB2y,OB2z,drag_K,seconds",
          {"U, Omega for B1=1 (neutral mode) and B2=1 (dipole mode); a=1, mu=1, e=z. "
           "Any squirmer: U = U_B1 + alpha*U_B2. drag_K = axial drag/(6 pi mu a U), sphere on axis only. "
           "a_over_R=0 means unbounded fluid."});
    const V3 e{0, 0, 1};
    for (double lam : lams) {
        std::unique_ptr<Tube> tube;
        if (lam > 0) tube = std::make_unique<Tube>(a / lam, LoR * a / lam, a, nth);
        for (double beta : betas) {
            Timer clock;
            double b = lam > 0 ? beta * (a / lam - a) : 0.0;
            SwimmerSolver s(tube.get(), icosphere(a, lev, {b, 0, 0}), {b, 0, 0});
            auto m = s.solve_free({s.slip(e, 1, 0), s.slip(e, 0, 1)});
            double K = NAN;
            if (towed && std::fabs(b) < 1e-12) K = -s.force(s.solve_prescribed({0, 0, 1}, {})).z / (6 * PI);
            k.row({lam, beta, (double)(s.Ns + s.Nt), m[0].U.x, m[0].U.y, m[0].U.z, m[0].Om.x, m[0].Om.y,
                   m[0].Om.z, m[1].U.x, m[1].U.y, m[1].U.z, m[1].Om.x, m[1].Om.y, m[1].Om.z, K,
                   clock.seconds()});
            std::printf("a/R=%.2f beta=%.2f: U_B1/U0=(%+.4f,%+.4f,%.4f) Om_B1=(%+.4f,%+.4f,%+.4f) "
                        "U_B2/U0=(%+.4f,%+.4f,%+.4f) K=%.4f [%.1fs]\n",
                        lam, beta, m[0].U.x / U0, m[0].U.y / U0, m[0].U.z / U0, m[0].Om.x, m[0].Om.y,
                        m[0].Om.z, m[1].U.x / U0, m[1].U.y / U0, m[1].U.z / U0, K, clock.seconds());
        }
    }
}



//============================//    
// fields                    //
//===========================// 

static void run_field(const Config& c, const std::string& out) {
    const double a = 1.0;
    double lam = c.num("a_over_R", 0.3), beta = c.num("beta", 0.0);
    int nth = c.integer("n_theta", 24), lev = c.integer("sphere_level", 3);
    int nx = c.integer("nx", 25), nz = c.integer("nz", 57);
    double zext = c.num("z_extent", 7.0);  // grid spans |z| <= z_extent (units of a)
    Timer clock;
    std::unique_ptr<Tube> tube;
    double R = lam > 0 ? a / lam : c.num("x_extent", 4.0) / 0.97;
    if (lam > 0) tube = std::make_unique<Tube>(R, c.num("L_over_R", 4.0) * R, a, nth);
    double b = lam > 0 ? beta * (R - a) : 0.0;
    const V3 e{0, 0, 1};
    SwimmerSolver s(tube.get(), icosphere(a, lev, {b, 0, 0}), {b, 0, 0});
    auto sl1 = s.slip(e, 1, 0), sl2 = s.slip(e, 0, 1);
    auto m = s.solve_free({sl1, sl2});
    std::vector<V3> pts;
    std::vector<int> inside;
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < nz; ++j) {
            V3 p{-0.97 * R + 1.94 * R * i / (nx - 1), 0, -zext + 2 * zext * j / (nz - 1)};
            inside.push_back(norm(p - V3{b, 0, 0}) < 1.03 * a ? 1 : 0);
            pts.push_back(p);
        }
    std::vector<V3> eval;
    for (size_t k = 0; k < pts.size(); ++k)
        if (!inside[k]) eval.push_back(pts[k]);
    for (int mode = 0; mode < 2; ++mode) {
        auto v = s.velocity(m[mode], mode == 0 ? sl1 : sl2, eval);
        Csv f(out + (mode == 0 ? "/field_B1.csv" : "/field_B2.csv"), "x,z,ux,uy,uz",
              {"lab-frame velocity in the plane y=0 for B" + std::to_string(mode + 1) +
                   "=1; points inside the sphere are NaN. a_over_R=" + std::to_string(lam) +
                   " beta=" + std::to_string(beta),
               "swimmer U=(" + std::to_string(m[mode].U.x) + "," + std::to_string(m[mode].U.y) + "," +
                   std::to_string(m[mode].U.z) + ")"});
        size_t q = 0;
        for (size_t k = 0; k < pts.size(); ++k) {
            if (inside[k]) f.row({pts[k].x, pts[k].z, NAN, NAN, NAN});
            else { f.row({pts[k].x, pts[k].z, v[q].x, v[q].y, v[q].z}); ++q; }
        }
    }
    Csv k(out + "/field_info.csv", "a_over_R,beta,R,UB1z,UB2z,OB1y");
    k.row({lam, beta, lam > 0 ? R : 0.0, m[0].U.z, m[1].U.z, m[0].Om.y});
    std::printf("field: %zu points x 2 modes, U_B1 = %.5f  (%.1f s)\n", eval.size(), m[0].U.z,
                clock.seconds());
}


//============================//    
// trajectory                 //
//===========================// 
static void run_traj(const Config& c, const std::string& out) {
    TrajectoryParams P;
    P.alpha = c.num("alpha", 0);
    P.a_over_R = c.num("a_over_R", 0.3);
    P.beta0 = c.num("beta0", 0.5);
    P.pitch0_deg = c.num("pitch0_deg", 0);
    P.yaw0_deg = c.num("yaw0_deg", 0);
    P.dt = c.num("dt", 0.5);
    P.t_max = c.num("t_max", 100);
    P.L_over_R = c.num("L_over_R", 4);
    P.n_theta = c.integer("n_theta", 30);
    P.sphere_level = c.integer("sphere_level", 2);
    P.beta_stop = c.num("beta_stop", 0.95);
    P.budget_s = c.num("budget_s", 0);
    ensure_dir(out);
    run_trajectory(P, out + "/" + c.str("name", "trajectory"));
}




int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <input-file> [key=value ...]\n", argv[0]);
        return 1;
    }
    try {
        Config c = Config::load(argv[1]);
        for (int i = 2; i < argc; ++i) c.parse_line(argv[i]);  // command-line overrides
        std::string mode = c.str("mode");
        std::string out = c.str("output_dir", "outputs/" + mode);
        ensure_dir(out);
        if (mode == "duct") run_duct(c, out);
        else if (mode == "kinematics") run_kinematics(c, out);
        else if (mode == "field") run_field(c, out);
        else if (mode == "trajectory") run_traj(c, out);
        else throw std::runtime_error("unknown mode '" + mode + "' (duct|kinematics|field|trajectory)");
    } catch (const std::exception& ex) {
        std::fprintf(stderr, "error: %s\n", ex.what());
        return 1;
    }
    return 0;
}
