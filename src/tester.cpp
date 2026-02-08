#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cmath>
#include <array>
#include "mesh_handler.hpp"
#include "fea_solver.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <fstream>
#include <sstream>
#include <CGAL/boost/graph/IO/STL.h>
#include <CLI/CLI.hpp>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/IO/polygon_mesh_io.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Surface_mesh;

struct BBox {
    double min_x, max_x;
    double min_y, max_y;
    double min_z, max_z;
};

struct TestNode {
    int id;
    double x, y, z;
};

struct TestStress {
    int id;
    double Sxx, Syy, Szz;
};

struct StressPoint {
    double x, y, z;
    double Sxx, Syy, Szz;
};

struct NodeStress {
    int id;
    std::array<double, 3> pos;     // x, y, z
    std::array<double, 3> stress;  // Sxx, Syy, Szz
};


// Key for 3D voxel index
struct VoxelKey {
    int i, j, k;

    bool operator==(const VoxelKey& other) const {
        return i == other.i && j == other.j && k == other.k;
    }
};

// Hash function for unordered_map
struct VoxelKeyHash {
    std::size_t operator()(const VoxelKey& v) const {
        return ((std::hash<int>()(v.i) ^
               (std::hash<int>()(v.j) << 1)) >> 1) ^
               (std::hash<int>()(v.k) << 1);
    }
};

struct AccumulatedStress {
    double pred[3] = {0,0,0};
    double truth[3] = {0,0,0};
    int count = 0;
};

// === CSV Loaders ===
std::vector<TestNode> load_csv_nodes(const std::filesystem::path& filename) {
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

std::vector<StressPoint> load_csv_stresses(const std::filesystem::path& filename, bool from_abaqus) {
    std::ifstream file(filename);
    std::vector<StressPoint> stresses;
    std::string line;
    std::getline(file, line); // skip header
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        StressPoint s;
        std::string val;
        std::getline(ss, val, ','); // instance -> element_id
        if (from_abaqus)
        {
            std::getline(ss, val, ','); // instance -> garbage
            std::getline(ss, val, ','); // instance -> MISES
        }
        std::getline(ss, val, ','); s.x = std::stod(val);
        std::getline(ss, val, ','); s.y = std::stod(val);
        std::getline(ss, val, ','); s.z = std::stod(val);
        std::getline(ss, val, ','); s.Sxx = std::stod(val);
        std::getline(ss, val, ','); s.Syy = std::stod(val);
        std::getline(ss, val, ','); s.Szz = std::stod(val);
        stresses.push_back(s);
    }
    return stresses;
}

/*
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
        ns.stress = {s.Sxx, s.Syy, s.Szz};
        merged.push_back(ns);
    }
    return merged;
}
*/

// === 3D Chunk Averaging ===
// === 3D Chunk Averaging ===
std::vector<NodeStress> chunk_average(
    const std::vector<StressPoint>& data,
    const BBox& bb,
    int nx, int ny, int nz)
{
    std::vector<NodeStress> averaged;
    if (data.empty() || nx <= 0 || ny <= 0 || nz <= 0) {
        return averaged;
    }

    double min_x = bb.min_x, max_x = bb.max_x;
    double min_y = bb.min_y, max_y = bb.max_y;
    double min_z = bb.min_z, max_z = bb.max_z;

    // Guard against degenerate bbox in any dimension
    if (max_x == min_x) max_x = min_x + 1e-9;
    if (max_y == min_y) max_y = min_y + 1e-9;
    if (max_z == min_z) max_z = min_z + 1e-9;

    double dx = (max_x - min_x) / nx;
    double dy = (max_y - min_y) / ny;
    double dz = (max_z - min_z) / nz;

    struct Voxel {
        int count = 0;
        double sum_Sxx = 0.0;
        double sum_Syy = 0.0;
        double sum_Szz = 0.0;
    };

    std::vector<Voxel> grid(nx * ny * nz);

    auto index = [&](int ix, int iy, int iz) {
        return ix + nx * (iy + ny * iz);
    };

    // 1) Assign points to voxels and accumulate stress
    for (const auto& sp : data) {
        int ix = std::min(int((sp.x - min_x) / dx), nx - 1);
        int iy = std::min(int((sp.y - min_y) / dy), ny - 1);
        int iz = std::min(int((sp.z - min_z) / dz), nz - 1);

        Voxel& v = grid[index(ix, iy, iz)];
        v.count++;
        v.sum_Sxx += sp.Sxx;
        v.sum_Syy += sp.Syy;
        v.sum_Szz += sp.Szz;
    }

    // 2) Build averaged nodes at voxel centers
    averaged.reserve(nx * ny * nz);
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int iz = 0; iz < nz; ++iz) {
                const Voxel& v = grid[index(ix, iy, iz)];
                if (v.count == 0) continue;

                NodeStress ns;
                ns.id = index(ix, iy, iz);

                ns.pos = {
                    min_x + (ix + 0.5) * dx,
                    min_y + (iy + 0.5) * dy,
                    min_z + (iz + 0.5) * dz
                };

                ns.stress[0] = v.sum_Sxx / v.count; // Sxx
                ns.stress[1] = v.sum_Syy / v.count; // Syy
                ns.stress[2] = v.sum_Szz / v.count; // Szz

                averaged.push_back(ns);
            }
        }
    }

    return averaged;
}

// === RMS Error ===
double rms_error(const std::vector<NodeStress>& pred,
                 const std::vector<NodeStress>& truth)
{
    double sum_sq = 0;
    int count = std::min(pred.size(), truth.size());

    std::ofstream out("rms_chunks.csv");
    out << "id,x,y,z,Sxx,Syy,Szz\n";


    size_t id = 1;
    for (int i = 0; i < count; ++i) {
        double diff[3];
        for (int j = 0; j < 3; ++j) {
             diff[j] = pred[i].stress[j] - truth[i].stress[j];
            sum_sq += diff[j] * diff[j];
        }

        out << id << ","
            << pred[i].pos[0] << ","
            << pred[i].pos[1] << ","
            << pred[i].pos[2] << ","
            << diff[0] << ","
            << diff[1] << ","
            << diff[2] << "\n";
    }

    out.close();

    return std::sqrt(sum_sq / (count * 3.0));
}



std::vector<Point> load_points(const std::filesystem::path& filename, bool includes_instance=true) {
    std::ifstream file(filename);
    std::vector<Point> pts;
    std::string line;
    std::getline(file, line); // skip header
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string val;
        int id;
        double x, y, z;

        if (includes_instance)
        {
            std::getline(ss, val, ','); //instance -> garbage
        }

        std::getline(ss, val, ','); id = std::stoi(val);
        std::getline(ss, val, ','); x = std::stod(val);
        std::getline(ss, val, ','); y = std::stod(val);
        std::getline(ss, val, ','); z = std::stod(val);
        pts.emplace_back(x, y, z);
    }
    return pts;
}


/**
 * Compare predicted vs truth only near the toolpath mesh.
 *
 * @param truth         Ground-truth stress points
 * @param pred          Predicted stress points
 * @param toolpath      Surface_mesh representing removed material boundary
 * @param radius        Distance (in mm) from toolpath to consider “near”
 *
 * @return RMS error over only near-cut stresses
 */
double rms_error_near_cut(
    const std::vector<NodeStress>& truth,
    const std::vector<NodeStress>& pred,
    const Surface_mesh& toolpath,
    double radius)
{
    // Build AABB tree for nearest point queries
    Tree tree(faces(toolpath).begin(), faces(toolpath).end(), toolpath);
    tree.accelerate_distance_queries();

    double sum_sq = 0.0;
    int count = 0;

    // Loop through ground truth nodes
    for (size_t i = 0; i < truth.size(); ++i)
    {
        const auto& t = truth[i];
        const auto& p = pred[i];

        Kernel::Point_3 q(t.pos[0], t.pos[1], t.pos[2]);

        // Compute squared distance to toolpath
        double dist2 = tree.squared_distance(q);

        if (dist2 <= radius * radius)
        {
            // Compare this node only
            for (int j = 0; j < 3; ++j)
            {
                double diff = p.stress[j] - t.stress[j];
                sum_sq += diff * diff;
            }
            count += 3;
        }
    }

    if (count == 0)
        return 0.0; // no near-cut nodes found

    return std::sqrt(sum_sq / count);
}

BBox compute_bbox(const std::vector<StressPoint>& a,
                  const std::vector<StressPoint>& b)
{
    BBox bb;
    bb.min_x = bb.min_y = bb.min_z =  1e18;
    bb.max_x = bb.max_y = bb.max_z = -1e18;

    auto update = [&](const StressPoint& p){
        bb.min_x = std::min(bb.min_x, p.x);
        bb.min_y = std::min(bb.min_y, p.y);
        bb.min_z = std::min(bb.min_z, p.z);
        bb.max_x = std::max(bb.max_x, p.x);
        bb.max_y = std::max(bb.max_y, p.y);
        bb.max_z = std::max(bb.max_z, p.z);
    };

    for (auto& p : a) update(p);
    for (auto& p : b) update(p);

    return bb;
}

void write_bins_to_csv(const std::vector<NodeStress>& bins,
                       const std::string& filename)
{
    std::ofstream out(filename);
    out << "id,x,y,z,Sxx,Syy,Szz\n";

    int id = 1;
    for (const auto& b : bins) {
        out << id++ << ","
            << b.pos[0] << ","
            << b.pos[1] << ","
            << b.pos[2] << ","
            << b.stress[0] << ","
            << b.stress[1] << ","
            << b.stress[2] << "\n";
    }
}


int main(int argc, char **argv) 
{
    CLI::App app{"MeshHandler Tester"};

    std::filesystem::path before_cut_nodes_csv, boundary_nodes_csv, pre_deformation_nodes_csv, post_deformation_nodes_csv, stress_values_csv;

    app.add_option("--before_cut_nodes", before_cut_nodes_csv)->required();
    app.add_option("--boundary_nodes", boundary_nodes_csv)->required();
    app.add_option("--predeformation_nodes", pre_deformation_nodes_csv)->required();
    app.add_option("--post_deformation_nodes", post_deformation_nodes_csv)->required();
    auto opt = app.add_option("--stress_values", stress_values_csv);

    CLI11_PARSE(app, argc, argv);

    bool use_stresses = opt->count() > 0;

    std::cout << "1" << std::endl;
    std::vector<Point> before_cut_points = load_points(before_cut_nodes_csv);
    std::cout << "2" << std::endl;
    std::vector<Point> post_deformation_points = load_points(post_deformation_nodes_csv);
    std::cout << "3" << std::endl;
    std::vector<Point> pre_deformation_points = load_points(pre_deformation_nodes_csv);
    std::cout << "4" << std::endl;
    std::vector<Point> boundary_points = load_points(boundary_nodes_csv);
    std::cout << "5" << std::endl;

    Surface_mesh before_cut_mesh, pre_deformation_mesh, post_deformation_mesh, tool_path_mesh;
    CGAL::convex_hull_3(before_cut_points.begin(), before_cut_points.end(), before_cut_mesh);
    CGAL::convex_hull_3(pre_deformation_points.begin(), pre_deformation_points.end(), pre_deformation_mesh);
    CGAL::convex_hull_3(post_deformation_points.begin(), post_deformation_points.end(), post_deformation_mesh);


    std::filesystem::path tool_path_stl("tool_path.stl");
    std::filesystem::path before_cut_stl("pre_cut.stl");
    std::filesystem::path pre_deformation_stl("pre_deformation.stl");
    std::filesystem::path post_deformation_stl("post_deformation.stl");

    CGAL::Polygon_mesh_processing::corefine_and_compute_difference(before_cut_mesh, pre_deformation_mesh, tool_path_mesh);

    if (!CGAL::IO::write_polygon_mesh(before_cut_stl.string(), before_cut_mesh, CGAL::parameters::stream_precision(17))) {
        throw std::runtime_error("Failed to write surface mesh to " + before_cut_stl.string());
    }

    if (!CGAL::IO::write_polygon_mesh(pre_deformation_stl.string(), pre_deformation_mesh, CGAL::parameters::stream_precision(17))) {
        throw std::runtime_error("Failed to write surface mesh to " + pre_deformation_stl.string());
    }


    if (!CGAL::IO::write_polygon_mesh(tool_path_stl.string(), tool_path_mesh, CGAL::parameters::stream_precision(17))) {
        throw std::runtime_error("Failed to write surface mesh to " + tool_path_stl.string());
    }

    if (!CGAL::IO::write_polygon_mesh(post_deformation_stl.string(), post_deformation_mesh, CGAL::parameters::stream_precision(17))) {
        throw std::runtime_error("Failed to write surface mesh to " + post_deformation_stl.string());
    }

    MeshHandler test_mesh(pre_deformation_mesh, boundary_points);
    test_mesh.set_stress_in_cut(tool_path_stl, post_deformation_stl);

    if (!CGAL::IO::write_polygon_mesh(tool_path_stl.string(), tool_path_mesh, CGAL::parameters::stream_precision(17))) {
        throw std::runtime_error("Failed to write surface mesh to " + tool_path_stl.string());
    }


    if (use_stresses)
    {
        std::vector<StressPoint> ground_truth_stresses = load_csv_stresses(stress_values_csv, true);
        std::vector<StressPoint> calculated_stresses = load_csv_stresses("stress_field.csv", false);

        int num_bins_x = 10;
        int num_bins_y = 10;
        int num_bins_z = 10;

        BBox bbox = compute_bbox(ground_truth_stresses, calculated_stresses);

        std::vector<NodeStress> ground_truth_chunks = chunk_average(ground_truth_stresses, bbox, num_bins_x, num_bins_y, num_bins_z);
        std::vector<NodeStress> calculated_chunks = chunk_average(calculated_stresses, bbox, num_bins_x, num_bins_y, num_bins_z); 

        double full_material_rms_error = rms_error(calculated_chunks, ground_truth_chunks);
        double near_cut_rms_error= rms_error_near_cut(calculated_chunks, ground_truth_chunks, tool_path_mesh, 1);

        std::cout << "Full Material RMS Error: " << full_material_rms_error << std::endl;
        std::cout << "Near Cut RMS Error: " << near_cut_rms_error << std::endl;

        write_bins_to_csv(calculated_chunks, "calculated_chunks.csv");
        write_bins_to_csv(ground_truth_chunks, "ground_truth_chunks.csv");
    }
}
