// test_duct_flow.cpp -- Stage 1 test: Poiseuille flow recovered from a pressure drop
//
// Checks
//   * flow rate Q vs Hagen-Poiseuille, and exact linearity Q ~ dP
//   * interior velocity profile at mid-length vs the parabola
//   * wall pressure (linear in z) and wall shear stress (-G R / 2)
//   * error decreases under mesh refinement
// Writes outputs/tests/duct_flow/*.csv   ->   python3 tests/plot_duct_flow.py
#include <vector>
#include "duct_flow.hpp"
#include "test_common.hpp"
#include "utils/io.hpp"

using namespace bem;
using namespace testing;

int main() {
    const std::string out = "outputs/tests/duct_flow";
    const double R = 1, L = 4, mu = 1;
    std::printf("== Stage 1: pressure-driven flow in a tube ==\n");
    Timer clock;
    DuctFlow duct(R, L, 24, 24, 6, mu);
    std::printf("mesh: %d panels (%.1f s)\n", duct.mesh.P(), clock.seconds());

    // --- flow rate and linearity
    Csv fq(out + "/flow_rate.csv", "dP,Q,Q_exact");
    double ratio0 = 0, max_dev = 0;
    for (double dP : {0.5, 1.0, 2.0, 4.0}) {
        auto s = duct.solve(dP, 0.0);
        fq.row({dP, s.Q, duct.exact_Q(dP)});
        double ratio = s.Q / dP;
        if (ratio0 == 0) ratio0 = ratio;
        max_dev = std::max(max_dev, std::fabs(ratio / ratio0 - 1));
    }
    auto s = duct.solve(1.0, 0.0);
    check_close("flow rate Q / Q_exact  (dP = 1)", s.Q / duct.exact_Q(1.0), 1.0, 1e-2);
    check_close("linearity: max deviation of Q/dP", max_dev, 0.0, 1e-10);

    // --- interior profile at mid-length along three directions
    std::vector<V3> pts;
    for (double ang : {0.0, PI / 5, PI / 2})
        for (int k = 0; k < 20; ++k) {
            double r = 0.95 * R * k / 19;
            pts.push_back({r * std::cos(ang), r * std::sin(ang), L / 2});
        }
    auto v = duct.velocity(s, pts);
    Csv fp(out + "/profile.csv", "r,uz,uz_exact,ur");
    double err = 0, err_r = 0, wmax = duct.exact_uz(0, 1.0);
    for (size_t k = 0; k < pts.size(); ++k) {
        double r = std::hypot(pts[k].x, pts[k].y);
        fp.row({r, v[k].z, duct.exact_uz(r, 1.0), std::hypot(v[k].x, v[k].y)});
        err = std::max(err, std::fabs(v[k].z - duct.exact_uz(r, 1.0)));
        err_r = std::max(err_r, std::hypot(v[k].x, v[k].y));
    }
    check_close("interior u_z: max error / u_max", err / wmax, 0.0, 2e-2);
    check_close("interior transverse velocity / u_max", err_r / wmax, 0.0, 1e-3);

    // --- wall traction: pressure p(z) = dP (1 - z/L), shear -G R/2
    write_panels(out + "/panels.csv", duct.mesh, s.u, s.f);
    double G = 1.0 / L, shear = 0, perr = 0;
    int n = 0;
    for (int p = 0; p < duct.mesh.P(); ++p) {
        V3 c = duct.mesh.centroid[p];
        if (duct.mesh.tags[p] != WALL || c.z < 0.5 || c.z > L - 0.5) continue;
        V3 er = unit(V3{c.x, c.y, 0});
        perr = std::max(perr, std::fabs(-dot(s.f[p], er) - (1.0 - G * c.z)));
        shear += s.f[p].z;
        ++n;
    }
    check_close("mean wall shear / (-G R / 2)", shear / n / (-G * R / 2), 1.0, 2e-2);
    check_close("wall pressure: max error / dP", perr, 0.0, 1e-2);

    // --- mesh convergence
    Csv fc(out + "/convergence.csv", "panels,h,rel_err_Q");
    std::vector<double> errs;
    for (auto m : std::vector<std::array<int, 3>>{{12, 12, 3}, {16, 16, 4}, {20, 20, 5}, {24, 24, 6}}) {
        DuctFlow d(R, L, m[0], m[1], m[2], mu);
        double e = std::fabs(d.solve(1.0, 0.0).Q / d.exact_Q(1.0) - 1);
        fc.row({(double)d.mesh.P(), std::sqrt(d.mesh.total_area() / d.mesh.P()), e});
        errs.push_back(e);
    }
    bool mono = true;
    for (size_t i = 1; i < errs.size(); ++i) mono &= errs[i] < errs[i - 1];
    check_true("error in Q decreases under refinement", mono);
    std::printf("total %.1f s; data in %s\n", clock.seconds(), out.c_str());
    return summary("test_duct_flow");
}
