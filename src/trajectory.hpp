// trajectory.hpp -- Stage 3: time integration of a squirmer's path in a tube
//
//   dX/dt = U(X, e),     de/dt = Omega(X, e) x e
//
// The tube translates axially with the swimmer (it is always centred at z = 0 in
// the computational frame), so the tube block is factorised only once; each step
// rebuilds only the swimmer couplings (see swimmer.hpp).
//
// Linearity in the squirming modes: U = U[B1=1] + alpha * U[B2=1] (and likewise
// Omega), so both modes are solved together as two right-hand sides.
//
// Time stepping: 4th-order Adams-Bashforth (Zhu, Lauga & Brandt 2013), started
// with three RK4 steps. Checkpoint/restart allows long runs to be split.
#pragma once
#include <array>
#include <fstream>
#include <string>
#include "swimmer.hpp"
#include "utils/io.hpp"

namespace bem {

struct TrajectoryParams {
    double alpha = 0;         // B2/B1 (pusher < 0 < puller)
    double a_over_R = 0.3;    // confinement
    double beta0 = 0.5;       // initial offset from the axis, in units of (R - a)
    double pitch0_deg = 0;    // initial tilt toward +x (radial direction)
    double yaw0_deg = 0;      // initial tilt toward +y (azimuthal direction)
    double dt = 0.5;          // time step, units a/B1
    double t_max = 100;
    double L_over_R = 4;      // truncated tube length
    int n_theta = 30;         // tube panels around the circumference
    int sphere_level = 2;     // icosphere refinement (2 -> 320 panels)
    double beta_stop = 0.95;  // stop when the swimmer gets this close to the wall
    double budget_s = 0;      // wall-clock budget per call (0 = unlimited)
};



// Runs (or resumes) a trajectory. Writes <out_prefix>.csv; keeps
// <out_prefix>.ckpt for restart. Returns a status string.
inline std::string run_trajectory(const TrajectoryParams& P, const std::string& out_prefix,
                                  bool verbose = true) {
    using State = std::array<double, 6>;  // x, y, z, ex, ey, ez
    Timer clock;
    const double a = 1, R = a / P.a_over_R;
    Tube tube(R, P.L_over_R * R, a, P.n_theta);
    Mesh sph0 = icosphere(a, P.sphere_level);
    int nevals = 0;

    struct Eval { State d; V3 U, Om; };
    auto evaluate = [&](const State& y) {
        V3 c{y[0], y[1], 0}, e = unit(V3{y[3], y[4], y[5]});
        Mesh S = sph0;
        S.translate(c);
        SwimmerSolver s(&tube, S, c);
        auto k = s.solve_free({s.slip(e, 1, 0), s.slip(e, 0, 1)});
        Eval r;
        r.U = k[0].U + P.alpha * k[1].U;
        r.Om = k[0].Om + P.alpha * k[1].Om;
        V3 de = cross(r.Om, e);
        r.d = {r.U.x, r.U.y, r.U.z, de.x, de.y, de.z};
        ++nevals;
        return r;
    };
    auto beta_of = [&](const State& y) { return std::hypot(y[0], y[1]) / (R - a); };


    // ---- initialise or restart
    std::string ck = out_prefix + ".ckpt", csvp = out_prefix + ".csv";
    State y;
    double t = 0;
    int step = 0;
    std::vector<State> hist;  // up to 3 previous derivatives, oldest first
    std::ifstream in(ck);
    bool restart = (bool)in;
    if (restart) {
        size_t nh;
        in >> t >> step >> nh;
        for (auto& v : y) in >> v;
        hist.resize(nh);
        for (auto& h : hist)
            for (auto& v : h) in >> v;
        if (verbose) std::printf("restart at t = %.3f (step %d)\n", t, step);
    } else {
        double b = P.beta0 * (R - a), xi = P.pitch0_deg * PI / 180, chi = P.yaw0_deg * PI / 180;
        y = {-b, 0, 0, std::cos(chi) * std::sin(xi), std::sin(chi), std::cos(chi) * std::cos(xi)};
    }
    char info[256];
    std::snprintf(info, sizeof info,
                  "alpha=%g a_over_R=%g beta0=%g pitch0_deg=%g yaw0_deg=%g dt=%g n_theta=%d sphere_level=%d",
                  P.alpha, P.a_over_R, P.beta0, P.pitch0_deg, P.yaw0_deg, P.dt, P.n_theta, P.sphere_level);
    Csv out(csvp, "t,x,y,z,ex,ey,ez,Ux,Uy,Uz,Omx,Omy,Omz,beta", {info}, restart);
    auto save = [&]() {
        std::ofstream o(ck);
        o.precision(17);
        o << t << " " << step << " " << hist.size() << "\n";
        for (double v : y) o << v << " ";
        o << "\n";
        for (auto& h : hist) {
            for (double v : h) o << v << " ";
            o << "\n";
        }
    };

    

    // ---- time loop
    std::string status = "finished";
    while (t < P.t_max - 1e-9) {
        if (beta_of(y) > P.beta_stop) { status = "wall contact"; break; }
        double el = clock.seconds();
        if (P.budget_s > 0 && nevals > 0 &&
            el + 1.3 * el / nevals * (hist.size() < 3 ? 4 : 1) > P.budget_s) {
            status = "paused";
            break;
        }
        Eval f = evaluate(y);
        out.row({t, y[0], y[1], y[2], y[3], y[4], y[5], f.U.x, f.U.y, f.U.z, f.Om.x, f.Om.y,
                 f.Om.z, beta_of(y)});
        if (hist.size() < 3) {  // RK4 start-up
            auto axpy = [](const State& s, double h, const State& d) {
                State o;
                for (int i = 0; i < 6; ++i) o[i] = s[i] + h * d[i];
                return o;
            };
            Eval k2 = evaluate(axpy(y, P.dt / 2, f.d));
            Eval k3 = evaluate(axpy(y, P.dt / 2, k2.d));
            Eval k4 = evaluate(axpy(y, P.dt, k3.d));
            for (int i = 0; i < 6; ++i)
                y[i] += P.dt / 6 * (f.d[i] + 2 * k2.d[i] + 2 * k3.d[i] + k4.d[i]);
        } else {  // Adams-Bashforth 4
            const State &f1 = hist[2], &f2 = hist[1], &f3 = hist[0];
            for (int i = 0; i < 6; ++i)
                y[i] += P.dt / 24 * (55 * f.d[i] - 59 * f1[i] + 37 * f2[i] - 9 * f3[i]);
        }
        double n = std::sqrt(y[3] * y[3] + y[4] * y[4] + y[5] * y[5]);
        for (int i = 3; i < 6; ++i) y[i] /= n;  // keep |e| = 1
        hist.push_back(f.d);
        if (hist.size() > 3) hist.erase(hist.begin());
        t += P.dt;
        ++step;
        save();
    }
    if (verbose)
        std::printf("%s: t = %.2f, beta = %.3f, z = %.2f  (%d solves, %.1f s, %.2f s/solve)\n",
                    status.c_str(), t, beta_of(y), y[2], nevals, clock.seconds(),
                    nevals ? clock.seconds() / nevals : 0.0);
    if (status != "paused") std::remove(ck.c_str());
    return status;
}

}  // namespace bem
