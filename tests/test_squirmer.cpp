// test_squirmer.cpp -- Stage 2 test: squirmer velocity and flow field
//
//   ./build/test_squirmer          quick (coarse meshes, ~30 s)
//   ./build/test_squirmer --full   production meshes (~3 min), used for the README figure
//
// Checks
//   1. unbounded fluid: Stokes drag 6 pi mu a U, swimming speed 2 B1/3, B2 mode
//      gives no motion, flow field = Blake (1971) analytical solution
//   2. sphere on the tube axis: towed drag vs Haberman & Sayre (1958); squirmer
//      speed decreases with confinement; B2 mode gives no motion
//   3. Lorentz reciprocal theorem reproduces the directly computed speed
//   4. off the axis: B1 mode rotates away from the wall and does not drift
//      sideways; B2 (puller) mode drifts away from the wall and does not rotate
//   5. flow field on a grid around the squirmer (for the plot)
// Writes outputs/tests/squirmer/*.csv   ->   python3 tests/plot_squirmer.py
#include <vector>
#include "swimmer.hpp"
#include "test_common.hpp"
#include "utils/io.hpp"

using namespace bem;
using namespace testing;

static double haberman_sayre(double l) {
    return (1 - 0.75857 * std::pow(l, 5)) /
           (1 - 2.1050 * l + 2.0865 * std::pow(l, 3) - 1.7068 * std::pow(l, 5) + 0.72603 * std::pow(l, 6));
}

// Blake's squirmer flow (lab frame), swimming along e = z, radius 1
static V3 blake(V3 x, double B1, double B2) {
    double r = norm(x), c = x.z / r, s = std::sqrt(std::max(0.0, 1 - c * c));
    double ur = 2. / 3 * B1 * c / std::pow(r, 3) + (std::pow(r, -4) - std::pow(r, -2)) * B2 * (3 * c * c - 1) / 2;
    double ut = 1. / 3 * B1 * s / std::pow(r, 3) + B2 * s * c / std::pow(r, 4);
    V3 rh = (1 / r) * x, eth = (1 / std::max(s, 1e-12)) * (c * rh - V3{0, 0, 1});
    return ur * rh + ut * eth;
}

int main(int argc, char** argv) {
    const bool full = has_flag(argc, argv, "--full");
    const int lev = full ? 3 : 2, nth_axis = 24, nth_off = full ? 36 : 24;
    const double a = 1, U0 = 2. / 3;
    const V3 e{0, 0, 1};
    const std::string out = "outputs/tests/squirmer";
    Timer clock;
    std::printf("== Stage 2: squirmer (%s mode: sphere level %d) ==\n", full ? "full" : "quick", lev);

    // ---------------------------------------------------------------- 1. free space
    std::printf("-- 1. unbounded fluid\n");
    {
        SwimmerSolver s(nullptr, icosphere(a, lev), {0, 0, 0});
        double K = -s.force(s.solve_prescribed({0, 0, 1}, {})).z / (6 * PI);
        auto m = s.solve_free({s.slip(e, 1, 0), s.slip(e, 0, 1)});
        check_close("drag / (6 pi mu a U)", K, 1.0, 5e-3);
        check_close("swimming speed U / (2 B1 / 3)", m[0].U.z / U0, 1.0, 6e-3);
        check_close("B2 mode: |U| + |Omega|", norm(m[1].U) + norm(m[1].Om), 0.0, 1e-10);
        // Blake field on shells r = 1.1 .. 4, pusher alpha = -3
        double alpha = -3, maxerr[3] = {0, 0, 0};
        std::vector<V3> pts;
        for (double r : {1.1, 1.5, 3.0})
            for (int i = 0; i < 24; ++i) {
                double th = PI * (i + 0.5) / 24;
                pts.push_back({r * std::sin(th), 0, r * std::cos(th)});
            }
        auto sl1 = s.slip(e, 1, 0), sl2 = s.slip(e, 0, 1);
        auto v1 = s.velocity(m[0], sl1, pts), v2 = s.velocity(m[1], sl2, pts);
        Csv fb(out + "/blake.csv", "r,theta,ur_bem,ut_bem,ur_blake,ut_blake",
               {"unbounded pusher alpha=-3, lab frame"});
        for (size_t k = 0; k < pts.size(); ++k) {
            V3 ub = v1[k] + alpha * v2[k], ue = blake(pts[k], 1, alpha);
            double r = norm(pts[k]);
            V3 rh = (1 / r) * pts[k], eth = unit(dot(rh, e) * rh - e);
            fb.row({r, std::acos(pts[k].z / r), dot(ub, rh), dot(ub, eth), dot(ue, rh), dot(ue, eth)});
            int sh = k / 24;
            maxerr[sh] = std::max(maxerr[sh], norm(ub - ue) / U0);
        }
        // Points closer to the surface than about one panel size are resolved only on
        // the fine mesh (quick mode: panels ~0.2a, so r = 1.1a is checked in --full only).
        if (full) check_close("Blake field error / U0 at r = 1.1a", maxerr[0], 0.0, 2e-2);
        check_close("Blake field error / U0 at r = 1.5a", maxerr[1], 0.0, full ? 1e-2 : 2e-2);
        check_close("Blake field error / U0 at r = 3.0a", maxerr[2], 0.0, full ? 5e-3 : 1e-2);
    }

    // ---------------------------------------------------------------- 2. on the axis
    std::printf("-- 2. sphere on the tube axis\n");
    Csv fa(out + "/axis.csv", "a_over_R,drag_K,drag_K_HS,U_B1,U_B2");
    double Uprev = U0;
    bool slower = true;
    double b2max = 0;
    for (double lam : {0.1, 0.2, 0.3, 0.4, 0.5, 0.6}) {
        Tube tube(a / lam, 4 * a / lam, a, nth_axis);
        SwimmerSolver s(&tube, icosphere(a, lev), {0, 0, 0});
        double K = -s.force(s.solve_prescribed({0, 0, 1}, {})).z / (6 * PI);
        auto m = s.solve_free({s.slip(e, 1, 0), s.slip(e, 0, 1)});
        fa.row({lam, K, haberman_sayre(lam), m[0].U.z, m[1].U.z});
        if (lam <= 0.4 + 1e-9)
            check_close("towed drag / Haberman-Sayre, a/R = " + std::to_string(lam).substr(0, 3), K / haberman_sayre(lam), 1.0, 1e-2);
        slower &= m[0].U.z < Uprev;
        Uprev = m[0].U.z;
        b2max = std::max(b2max, norm(m[1].U) + norm(m[1].Om));
    }
    check_true("squirmer speed decreases monotonically with a/R", slower);
    check_close("B2 mode on axis: max |U| + |Omega|", b2max, 0.0, 1e-6);

    // ---------------------------------------------------------------- 3. reciprocal theorem
    std::printf("-- 3. Lorentz reciprocal theorem (a/R = 0.3)\n");
    const double lam = 0.3, R = a / lam;
    {
        Tube tube(R, 4 * R, a, nth_axis);
        SwimmerSolver s(&tube, icosphere(a, lev), {0, 0, 0});
        auto sl = s.slip(e, 1, 0);
        auto m = s.solve_free({sl});
        auto tow = s.solve_prescribed({0, 0, 1}, {});
        double Fz = s.force(tow).z, num = 0;
        for (int p = 0; p < s.Ns; ++p) num += dot(sl[p], tow.f_sphere[p]) * s.S.area[p];
        double Urec = num / Fz;  // U.F_hat = -\int u_s . t_hat dS,  t_hat = -f
        Csv fr(out + "/reciprocal.csv", "U_direct,U_reciprocal");
        fr.row({m[0].U.z, Urec});
        check_close("U(reciprocal theorem) / U(direct)", Urec / m[0].U.z, 1.0, 1e-2);

        // ------------------------------------------------------------ 5. flow field
        std::printf("-- 5. flow field around the squirmer (a/R = 0.3, on axis)\n");
        const int nx = 25, nz = 57;
        std::vector<V3> pts;
        std::vector<int> in;
        for (int i = 0; i < nx; ++i)
            for (int j = 0; j < nz; ++j) {
                V3 p{-0.97 * R + 1.94 * R * i / (nx - 1), 0, -2.2 * R + 4.4 * R * j / (nz - 1)};
                in.push_back(norm(p) < 1.03 * a);
                if (!in.back()) pts.push_back(p);
            }
        auto m2 = s.solve_free({s.slip(e, 0, 1)});
        auto v1 = s.velocity(m[0], sl, pts), v2 = s.velocity(m2[0], s.slip(e, 0, 1), pts);
        Csv ff(out + "/field.csv", "x,z,ux_B1,uz_B1,ux_B2,uz_B2",
               {"a_over_R=0.3 on axis, lab frame, R=" + std::to_string(R)});
        size_t q = 0;
        double umax = 0;
        for (int i = 0, k = 0; i < nx; ++i)
            for (int j = 0; j < nz; ++j, ++k) {
                double x = -0.97 * R + 1.94 * R * i / (nx - 1), z = -2.2 * R + 4.4 * R * j / (nz - 1);
                if (in[k]) { ff.row({x, z, NAN, NAN, NAN, NAN}); continue; }
                ff.row({x, z, v1[q].x, v1[q].z, v2[q].x, v2[q].z});
                if (std::fabs(z) > 2 * R) umax = std::max(umax, norm(v1[q]));
                ++q;
            }
        check_close("flow decays along the tube: max |u|/U0 at |z| > 2R", umax / U0, 0.0, 2e-2);
    }

    // ---------------------------------------------------------------- 4. off the axis
    std::printf("-- 4. off-axis kinematics (a/R = 0.3, wall on the +x side)\n");
    {
        Tube tube(R, 4 * R, a, nth_off);
        Csv fo(out + "/offaxis.csv", "beta,UB1x,UB1z,OB1y,UB2x,UB2z,OB2y");
        bool rot_away = true, puller_away = true;
        double side_B1 = 0, rot_B2 = 0, Uz_prev = 1e9;
        bool uz_dec = true;
        for (double beta : {0.0, 0.2, 0.4, 0.6, 0.75, 0.85}) {
            double b = beta * (R - a);
            SwimmerSolver s(&tube, icosphere(a, lev, {b, 0, 0}), {b, 0, 0});
            auto m = s.solve_free({s.slip(e, 1, 0), s.slip(e, 0, 1)});
            fo.row({beta, m[0].U.x, m[0].U.z, m[0].Om.y, m[1].U.x, m[1].U.z, m[1].Om.y});
            if (beta > 0) {
                rot_away &= m[0].Om.y < 0;      // e turns toward -x, away from the wall
                puller_away &= m[1].U.x < 0;    // B2 > 0 (puller) moves toward -x
            }
            side_B1 = std::max(side_B1, std::fabs(m[0].U.x));
            rot_B2 = std::max(rot_B2, std::fabs(m[1].Om.y));
            uz_dec &= m[0].U.z < Uz_prev + 1e-9;
            Uz_prev = m[0].U.z;
        }
        check_true("B1 mode rotates away from the nearest wall", rot_away);
        check_true("B2 > 0 (puller) drifts away from the wall", puller_away);
        check_close("B1 mode: max transverse velocity", side_B1, 0.0, 1e-8);
        check_close("B2 mode: max rotation rate", rot_B2, 0.0, 1e-8);
        check_true("axial speed decreases with eccentricity", uz_dec);
    }
    std::printf("total %.1f s; data in %s\n", clock.seconds(), out.c_str());
    return summary("test_squirmer");
}
