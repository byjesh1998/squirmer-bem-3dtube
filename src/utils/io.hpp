// io.hpp -- CSV output, directories, timing
#pragma once
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <string>
#include <vector>
#include "mesh.hpp"

namespace bem {

inline void ensure_dir(const std::string& d) { std::filesystem::create_directories(d); }

// CSV file with optional '# ...' comment lines followed by one header line.
class Csv {
public:
    Csv(const std::string& path, const std::string& header,
        const std::vector<std::string>& comments = {}, bool append = false) {
        std::filesystem::path p(path);
        if (p.has_parent_path()) ensure_dir(p.parent_path().string());
        bool exists = std::filesystem::exists(p);
        out_.open(path, append ? std::ios::app : std::ios::trunc);
        out_.precision(10);
        if (!append || !exists) {
            for (auto& c : comments) out_ << "# " << c << "\n";
            out_ << header << "\n";
        }
    }
    void row(std::initializer_list<double> v) { row(std::vector<double>(v)); }
    void row(const std::vector<double>& v) {
        for (size_t i = 0; i < v.size(); ++i) out_ << (i ? "," : "") << v[i];
        out_ << "\n";
        out_.flush();
    }

private:
    std::ofstream out_;
};

// one row per panel: geometry, tag, velocity and traction
inline void write_panels(const std::string& path, const Mesh& M, const std::vector<V3>& u,
                         const std::vector<V3>& f) {
    Csv c(path, "cx,cy,cz,nx,ny,nz,area,tag,ux,uy,uz,fx,fy,fz");
    for (int p = 0; p < M.P(); ++p)
        c.row({M.centroid[p].x, M.centroid[p].y, M.centroid[p].z, M.normal[p].x, M.normal[p].y,
               M.normal[p].z, M.area[p], (double)M.tags[p], u[p].x, u[p].y, u[p].z, f[p].x,
               f[p].y, f[p].z});
}

struct Timer {
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    double seconds() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
};

}  // namespace bem
