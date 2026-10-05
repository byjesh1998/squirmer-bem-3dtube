// test_trajectory.cpp -- Stage 3 test (optional extension): swimmer trajectories
//
//   ./build/test_trajectory          quick (~1.5 min): neutral half-wave + pusher crash
//   ./build/test_trajectory --full   the four production runs of the README
//                                    (~30 min on one core; resumable)
//
// Checks (quick)
//   * neutral squirmer (alpha = 0) started at beta = 0.7 parallel to the axis
//     crosses the axis and turns at the SAME distance on the other side
//     (amplitude A = 2 b_I, Zhu, Lauga & Brandt 2013)
//   * pusher (alpha = -3) started at beta = 0.3 reaches the wall
// Writes outputs/tests/trajectories/*.csv  ->  python3 tests/plot_trajectories.py
#include <fstream>
#include <sstream>
#include <vector>
#include "test_common.hpp"
#include "trajectory.hpp"

using namespace bem;
using namespace testing;

// read the numeric rows of a trajectory CSV (skips '#' comments and header)
static std::vector<std::vector<double>> read_csv(const std::string& path) {
    std::ifstream in(path);
    std::vector<std::vector<double>> rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 't') continue;
        std::vector<double> r;
        std::stringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, ',')) r.push_back(std::stod(tok));
        rows.push_back(r);
    }
    return rows;
}

int main(int argc, char** argv) {
    const std::string out = "outputs/tests/trajectories";
    ensure_dir(out);
    Timer clock;
    if (has_flag(argc, argv, "--full")) {
        std::printf("== Stage 3 (full): production trajectories, a/R = 0.3 ==\n");
        struct Case { const char* name; double alpha, beta0, dt, tmax; };
        for (Case c : {Case{"neutral", 0, 0.7, 0.5, 150}, Case{"pusher", -3, 0.3, 0.5, 200},
                       Case{"puller_a3", 3, 0.7, 1.0, 124}, Case{"puller_a5", 5, 0.3, 1.0, 250}}) {
            TrajectoryParams P;
            P.alpha = c.alpha; P.beta0 = c.beta0; P.dt = c.dt; P.t_max = c.tmax;
            std::printf("-- %s\n", c.name);
            run_trajectory(P, out + "/" + c.name);
        }
        std::printf("total %.1f s; data in %s\n", clock.seconds(), out.c_str());
        return 0;
    }

    std::printf("== Stage 3 (quick): trajectories, a/R = 0.3 ==\n");
    const double R = 1 / 0.3, bmax = R - 1;

    std::printf("-- neutral squirmer, one half-wave\n");
    TrajectoryParams P;
    P.alpha = 0; P.beta0 = 0.7; P.dt = 1.0; P.t_max = 44; P.n_theta = 24;
    std::remove((out + "/neutral_quick.ckpt").c_str());
    run_trajectory(P, out + "/neutral_quick");
    auto d = read_csv(out + "/neutral_quick.csv");
    double xmax = -1e9, xmin = 1e9;
    for (auto& r : d) { xmax = std::max(xmax, r[1]); xmin = std::min(xmin, r[1]); }
    double b0 = P.beta0 * bmax;
    check_true("neutral: crosses the axis", xmin < 0 && xmax > 0);
    check_close("neutral: turning distance / initial distance", xmax / b0, 1.0, 1e-2);
    check_true("neutral: stays clear of the wall", xmax < 0.8 * bmax);

    std::printf("-- pusher (alpha = -3)\n");
    TrajectoryParams Q;
    Q.alpha = -3; Q.beta0 = 0.3; Q.dt = 1.0; Q.t_max = 60; Q.n_theta = 24;
    std::remove((out + "/pusher_quick.ckpt").c_str());
    std::string st = run_trajectory(Q, out + "/pusher_quick");
    check_true("pusher: reaches the wall (beta > 0.95) before t = 60", st == "wall contact", st);

    std::printf("total %.1f s; data in %s\n", clock.seconds(), out.c_str());
    return summary("test_trajectory");
}
