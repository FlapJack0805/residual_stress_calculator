#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cmath>
#include <array>
#include "mesh_handler.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <fstream>
#include <sstream>
#include <CGAL/boost/graph/IO/STL.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Surface_mesh;

struct TestNode {
    int id;
    double x, y, z;
};

struct TestStress {
    int id;
    double Sxx, Syy, Szz, Sxy, Sxz, Syz;
};

struct NodeStress {
    int id;
    std::array<double, 3> pos;     // x, y, z
    std::array<double, 6> stress;  // Sxx, Syy, Szz, Sxy, Sxz, Syz
};

// === CSV Loaders ===
std::vector<TestNode> load_csv_nodes(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<TestNode> nodes;
    std::string line;
    std::getline(file, line); // skip header
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        TestNode node;
        std::string val;
        std::getline(ss, val, ','); node.id = std::stoi(val);
        std::getline(ss, val, ','); node.x = std::stod(val);
        std::getline(ss, val, ','); node.y = std::stod(val);
        std::getline(ss, val, ','); node.z = std::stod(val);
        nodes.push_back(node);
    }
    return nodes;
}

std::vector<TestStress> load_csv_stresses(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<TestStress> stresses;
    std::string line;
    std::getline(file, line); // skip header
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        TestStress s;
        std::string val;
        std::getline(ss, val, ','); s.id = std::stoi(val);
        std::getline(ss, val, ','); s.Sxx = std::stod(val);
        std::getline(ss, val, ','); s.Syy = std::stod(val);
        std::getline(ss, val, ','); s.Szz = std::stod(val);
        std::getline(ss, val, ','); s.Sxy = std::stod(val);
        std::getline(ss, val, ','); s.Sxz = std::stod(val);
        std::getline(ss, val, ','); s.Syz = std::stod(val);
        stresses.push_back(s);
    }
    return stresses;
}

// === Merge ===
std::vector<NodeStress> merge_node_stresses(
    const std::vector<TestNode>& nodes,
    const std::vector<TestStress>& stresses)
{
    std::unordered_map<int, TestStress> stress_map;
    for (const auto& s : stresses)
        stress_map[s.id] = s;

    std::vector<NodeStress> merged;
    merged.reserve(nodes.size());
    for (const auto& n : nodes) {
        if (!stress_map.count(n.id)) continue;
        const auto& s = stress_map[n.id];
        NodeStress ns;
        ns.id = n.id;
        ns.pos = {n.x, n.y, n.z};
        ns.stress = {s.Sxx, s.Syy, s.Szz, s.Sxy, s.Sxz, s.Syz};
        merged.push_back(ns);
    }
    return merged;
}

// === 1D Chunk Averaging ===
std::vector<NodeStress> chunk_average(
    const std::vector<NodeStress>& data,
    char axis, int bins)
{
    double min_c = 1e9, max_c = -1e9;
    for (const auto& ns : data) {
        double c = (axis == 'x') ? ns.pos[0] :
                   (axis == 'y') ? ns.pos[1] : ns.pos[2];
        min_c = std::min(min_c, c);
        max_c = std::max(max_c, c);
    }

    double bin_width = (max_c - min_c) / bins;
    struct Bin {
        int count = 0;
        std::array<double, 6> stress_sum{};
        double center;
    };
    std::vector<Bin> binned(bins);
    for (int i = 0; i < bins; ++i)
        binned[i].center = min_c + (i + 0.5) * bin_width;

    for (const auto& ns : data) {
        double c = (axis == 'x') ? ns.pos[0] :
                   (axis == 'y') ? ns.pos[1] : ns.pos[2];
        int idx = std::min(int((c - min_c) / bin_width), bins - 1);
        binned[idx].count++;
        for (int j = 0; j < 6; ++j)
            binned[idx].stress_sum[j] += ns.stress[j];
    }

    std::vector<NodeStress> averaged;
    for (const auto& bin : binned) {
        if (bin.count == 0) continue;
        NodeStress ns;
        ns.pos = {0.0, 0.0, 0.0};
        if (axis == 'x') ns.pos[0] = bin.center;
        else if (axis == 'y') ns.pos[1] = bin.center;
        else ns.pos[2] = bin.center;

        for (int j = 0; j < 6; ++j)
            ns.stress[j] = bin.stress_sum[j] / bin.count;
        averaged.push_back(ns);
    }
    return averaged;
}

// === RMS Error ===
double rms_error(const std::vector<NodeStress>& pred,
                 const std::vector<NodeStress>& truth)
{
    double sum_sq = 0;
    int count = std::min(pred.size(), truth.size());

    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < 6; ++j) {
            double diff = pred[i].stress[j] - truth[i].stress[j];
            sum_sq += diff * diff;
        }
    }
    return std::sqrt(sum_sq / (count * 6.0));
}



std::vector<Point> load_points(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<Point> pts;
    std::string line;
    std::getline(file, line); // skip header
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string val;
        int id;
        double x, y, z;
        std::getline(ss, val, ','); id = std::stoi(val);
        std::getline(ss, val, ','); x = std::stod(val);
        std::getline(ss, val, ','); y = std::stod(val);
        std::getline(ss, val, ','); z = std::stod(val);
        pts.emplace_back(x, y, z);
    }
    return pts;
}

int main() {
    std::vector<Point> points = load_points("abaqus_before_cut_nodes.csv");

    Surface_mesh mesh;
    CGAL::convex_hull_3(points.begin(), points.end(), mesh);

    std::cout << "Convex hull has " 
              << mesh.number_of_vertices() << " vertices and "
              << mesh.number_of_faces() << " faces.\n";

    std::ofstream before_cut_mesh("before_cut_mesh.stl");
    std::ofstream boundary_condition("boundary_condition.stl");

    MeshHandler test_mesh(before_cut_mesh, boundary_condition);
}

