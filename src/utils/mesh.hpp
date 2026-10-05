////////////////////////////////////////////////////////////////////////////////
// // mesh.hpp -- flat-triangle surface meshes: closed tube and icosphere
//
// Convention used everywhere in the code: panel normals point OUT OF THE FLUID.
//   tube wall: radially outward      tube caps: -e_z (inlet), +e_z (outlet)
//   sphere   : toward the sphere centre
//
// Tags: 0 = tube wall, 1 = inlet cap (z = z_min), 2 = outlet cap (z = z_max),
//       10 = sphere (swimmer)
////////////////////////////////////////////////////////////////////////////////

//                  MESH GENERATION
//                        │
//           ┌────────────┴────────────┐
//           │                         │
//       Tube mesh                Icosphere
//           │                         │
//     ┌─────┼─────┐              subdivision
//     │     │     │                   │
//   wall  inlet outlet                │
//     │     │     │                   │
//     └─────┼─────┘                   │
//           │                         │
//        merge                    volume scale
//           │                         │
//           └──────────┬──────────────┘
//                      │
//                   geometry
//                      │
//           ┌──────────┼──────────┐
//           │          │          │
//         area       normal    centroid


#pragma once
#include <algorithm>
#include <array>
#include <functional>
#include <map>
#include <vector>
#include "vec3.hpp"

namespace bem {

enum Tag { WALL = 0, INLET = 1, OUTLET = 2, SPHERE = 10 };

//============================//    
// defining meshes            //
//===========================// 
struct Mesh {
    std::vector<V3> nodes;
    std::vector<std::array<int, 3>> tris;
    std::vector<int> tags;
    // derived geometry (filled by geometry())
    std::vector<std::array<V3, 3>> V;  // vertex coordinates per panel
    std::vector<V3> normal, centroid;
    std::vector<double> area, h;       // h = sqrt(2*area), a panel length scale

    int P() const { return (int)tris.size(); }

    void geometry() {
        int n = P();
        V.resize(n); normal.resize(n); centroid.resize(n); area.resize(n); h.resize(n);
        for (int p = 0; p < n; ++p) {
            for (int k = 0; k < 3; ++k) V[p][k] = nodes[tris[p][k]];
            V3 c = cross(V[p][1] - V[p][0], V[p][2] - V[p][0]); //area by cross product
            area[p] = 0.5 * norm(c);
            normal[p] = (1.0 / (2 * area[p])) * c;
            centroid[p] = (1.0 / 3) * (V[p][0] + V[p][1] + V[p][2]);
            h[p] = std::sqrt(2 * area[p]);
        }
    }


    // flip triangles so that normal . outward(centroid) > 0
    void orient(const std::function<V3(V3)>& outward) {
        geometry();
        for (int p = 0; p < P(); ++p)
            if (dot(normal[p], outward(centroid[p])) < 0) std::swap(tris[p][1], tris[p][2]);
        geometry();//recomputes all geometric quantities after the orientation change.
    }
    void translate(V3 d) {
        for (auto& x : nodes) x = x + d;
        geometry();
    }
    double total_area() const {
        double s = 0;
        for (double A : area) s += A; //calculates total area
        return s;
    }
};



inline Mesh merge(const Mesh& a, const Mesh& b) { //This combines two independent meshes into one.
    Mesh m = a;
    int off = (int)a.nodes.size();
    m.nodes.insert(m.nodes.end(), b.nodes.begin(), b.nodes.end());
    for (auto t : b.tris) m.tris.push_back({t[0] + off, t[1] + off, t[2] + off});
    m.tags.insert(m.tags.end(), b.tags.begin(), b.tags.end());
    m.geometry();
    return m;
}




//============================//    
// axial nodes.               //
//===========================// 
// creates equally spaced points along the tube axis.
inline std::vector<double> uniform_z(double z0, double z1, int nz) {
    std::vector<double> z(nz + 1);
    for (int k = 0; k <= nz; ++k) z[k] = z0 + (z1 - z0) * k / nz;
    return z; // produce axial points like (-5,-2.5,0,2.5,5)
}


//creates equally spaced points along the tube axis.
// Nodes on [-L/2, L/2]: spacing dzmin for |z| < core, then growing geometrically
// (factor `growth`) up to dzmax toward the ends.
inline std::vector<double> graded_z(double L, double dzmin, double dzmax, double core,
                                    double growth = 1.25) {
    std::vector<double> half = {0.0};
    double dz = dzmin;
    while (half.back() < L / 2 - 1e-12) {
        if (half.back() >= core) dz = std::min(dz * growth, dzmax);
        half.push_back(std::min(half.back() + dz, L / 2));
    }
    if (half.size() > 2 && L / 2 - half[half.size() - 2] < 0.3 * dz)
        half.erase(half.end() - 2);  // avoid a sliver panel at the end
    std::vector<double> z;
    for (int i = (int)half.size() - 1; i >= 0; --i) z.push_back(-half[i]);
    for (size_t i = 1; i < half.size(); ++i) z.push_back(half[i]);
    return z;
}




//============================//    
// Tube                      //
//===========================// 
        //      outlet
        //   +-----------+
        //  /             \
        // |               |
        // |      fluid    |
        // |               |
        //  \             /
        //   +-----------+
        //      inlet
//It consists of: cylindrical wall, inlet disk, outlet disk
// Closed tube of radius R, axis z, axial nodes z[], n_theta panels around, n_r
// rings on each cap. The polygon vertices are placed on a slightly larger radius
// so that the faceted cross-section has the same AREA as the circle (removes the
// leading O(h^2) geometric error of flat panels; see README).
inline Mesh tube_mesh(double R, const std::vector<double>& z, int nth, int nr) {
    double al = 2 * PI / nth; // the number of angular divisions
// The actual polygon vertices are placed at radius `Rp`, not `R`.
// A flat-sided polygon inscribed in a circle of radius `R` has an area
// smaller than the true circular area:
//Rp = R * sqrt(dtheta / sin(dtheta))
    double Rp = R * std::sqrt(al / std::sin(al));
    int nz = (int)z.size() - 1;
    Mesh wall;
    for (int i = 0; i < nth; ++i)
        for (int k = 0; k <= nz; ++k)
            wall.nodes.push_back({Rp * std::cos(i * al), Rp * std::sin(i * al), z[k]});//(x,y,z)=(R_p costheta,R_p sintheta,zk)
    auto id = [&](int i, int k) { return (i % nth) * (nz + 1) + k; };
    for (int i = 0; i < nth; ++i)
        for (int k = 0; k < nz; ++k) {
            int a = id(i, k), b = id(i + 1, k), c = id(i + 1, k + 1), d = id(i, k + 1);
            wall.tris.push_back({a, b, c}); //each rectangular surface cell becomes two flat triangles.
            wall.tris.push_back({a, c, d});
        }
    wall.tags.assign(wall.tris.size(), WALL);
    wall.orient([](V3 x) { return V3{x.x, x.y, 0}; });

    auto cap = [&](double z0, int tag, double sg) {//z0 → z-coordinate of the cap, tag → INLET or OUTLET, sg → desired normal direction, -1 or +1
        Mesh m;
        m.nodes.push_back({0, 0, z0});
        for (int k = 1; k <= nr; ++k)
            for (int i = 0; i < nth; ++i)
                m.nodes.push_back({Rp * k / nr * std::cos(i * al), Rp * k / nr * std::sin(i * al), z0});
        auto ring = [&](int k, int i) { return 1 + k * nth + (i % nth); };
        for (int i = 0; i < nth; ++i) m.tris.push_back({0, ring(0, i), ring(0, i + 1)});
        for (int k = 0; k < nr - 1; ++k)
            for (int i = 0; i < nth; ++i) {
                int a = ring(k, i), b = ring(k, i + 1), c = ring(k + 1, i + 1), d = ring(k + 1, i);
                m.tris.push_back({a, d, c});
                m.tris.push_back({a, c, b});
            }
        m.tags.assign(m.tris.size(), tag);
        m.orient([sg](V3) { return V3{0, 0, sg}; });
        return m;
    };
    return merge(merge(wall, cap(z.front(), INLET, -1.0)), cap(z.back(), OUTLET, 1.0));
}

// Tube for swimmer problems: centred on z = 0, length L, axially graded
// (fine within |z| < a + R/2, where the swimmer sits).
inline Mesh swimmer_tube_mesh(double R, double L, double a, int nth, int nr = 4) { //near the swimmer → fine axial resolution, farther away → progressively coarser resolution
    double h = 2 * PI * R / nth;//The fine region extends to approximately: |z| < a+\frac{R}{2}.
    return tube_mesh(R, graded_z(L, 0.6 * h, 1.4 * h, a + 0.5 * R), nth, nr);//So if the swimmer has characteristic length a, the tube gets finer around it.
}




//============================//    
// Sphere                     //
//===========================// 
// Subdivided icosahedron (20 * 4^level panels), scaled so the polyhedron has the
// same VOLUME as the sphere. Normals point into the sphere (out of the fluid).
// An icosahedron produces much more uniform triangular elements than ordinary spherical latitude/longitude grids.
inline Mesh icosphere(double a, int level, V3 center = {}, int tag = SPHERE) {
    double t = (1 + std::sqrt(5.0)) / 2; // golden ratio
    std::vector<V3> v = {{-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0}, {0, -1, t}, {0, 1, t}, //12 standard vertices
                         {0, -1, -t}, {0, 1, -t}, {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1}};
    for (auto& x : v) x = unit(x);
    std::vector<std::array<int, 3>> F = {
        {0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4},
        {11, 10, 2}, {10, 7, 6}, {7, 1, 8}, {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8},
        {3, 8, 9}, {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};
    for (int l = 0; l < level; ++l) {
        std::map<std::pair<int, int>, int> cache;
        auto mid = [&](int i, int j) {
            auto k = std::make_pair(std::min(i, j), std::max(i, j));
            auto it = cache.find(k);
            if (it != cache.end()) return it->second;
            v.push_back(unit(v[i] + v[j]));
            return cache[k] = (int)v.size() - 1;
        };
        std::vector<std::array<int, 3>> nf;// subdivision
        for (auto f : F) {
            int ab = mid(f[0], f[1]), bc = mid(f[1], f[2]), ca = mid(f[2], f[0]);
            nf.push_back({f[0], ab, ca});
            nf.push_back({f[1], bc, ab});
            nf.push_back({f[2], ca, bc});
            nf.push_back({ab, bc, ca});
        }
        F = nf;
    }
    double vol = 0;
    for (auto f : F) vol += dot(v[f[0]], cross(v[f[1]], v[f[2]])); // projecting the midpoint onto spheres
    double s = std::cbrt(4 * PI / 3 / (std::fabs(vol) / 6));
    Mesh m;
    for (auto x : v) m.nodes.push_back(center + (a * s) * x);
    m.tris = F;
    m.tags.assign(F.size(), tag);
    m.orient([center](V3 x) { return center - x; });
    return m;
}

}  // namespace bem
