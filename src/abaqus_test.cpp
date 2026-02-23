#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <sstream>
#include <tuple>
#include "mesh_handler.hpp"
#include <CLI/CLI.hpp>
#include "../../surfacic_toolpaths/include/toolpath.hxx"

#define LENGTH_OF_CUT 10 // TODO: Find out a reasonable size for each cut

namespace fs = std::filesystem;


struct Config 
{
    std::string abaqus_cmd = "abaqus";
    std::string abaqus_script = "abaqus_tester.py"; // python file
    fs::path case_dir = "case_000";
    int seed = 7;

    int max_steps = 50;
};

struct Position
{
    double x;
    double y;
    double z;

    Position operator+(const Position& other)
    {
        return {x + other.x, y + other.y , z + other.z};
    }

    Position operator-(const Position& other)
    {
        return {x - other.x, y - other.y , z - other.z};
    }
};

// -------------------------------
// Simple utilities
// -------------------------------

Point3D position_to_point3D(Position pos)
{
    return Point3D{pos.x, pos.y, pos.z};
}

Vec3D delta_to_vec3D(const Position& a, const Position& b) 
{
    return Vec3D{b.x - a.x, b.y - a.y, b.z - a.z};
}

double wrap_0_2pi(double a)
{
    const double TWO_PI = 2.0 * M_PI;
    while (a < 0) a += TWO_PI;
    while (a >= TWO_PI) a -= TWO_PI;
    return a;
}

std::pair<Position, double> arcMidpoint(const Position& center, const Position& start, const Position& end, bool ccw)
{
    double th1 = atan2(start.y - center.y, start.x - center.x);
    double th2 = atan2(end.y - center.y, end.x - center.x);

    double delta = ccw
        ? wrap_0_2pi(th2 - th1)
        : wrap_0_2pi(th1 - th2);

    double thm = ccw
        ? th1 + delta * 0.5
        : th1 - delta * 0.5;

    double R = hypot(start.x - center.x, start.y - center.y);

    double arc_length = R * delta;

    return {{
        center.x + R * cos(thm),
        center.y + R * sin(thm),
        center.z
    }, arc_length};
}

static void write_text_file(const fs::path& p, const std::string& s) 
{
    std::ofstream out(p);
    if (!out) throw std::runtime_error("Failed to open for write: " + p.string());
    out << s;
}

static std::string read_text_file(const fs::path& p) 
{
    std::ifstream in(p);
    if (!in) throw std::runtime_error("Failed to open for read: " + p.string());
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

static void ensure_dir(const fs::path& p) 
{
    fs::create_directories(p);
}

static void require_exists(const fs::path& p, const std::string& msg) 
{
    if (!fs::exists(p)) throw std::runtime_error(msg + ": " + p.string());
}

// -------------------------------
// Abaqus process launcher
// -------------------------------
static int run_cmd(const std::string& cmd, const fs::path& workdir) 
{
    // run command in working directory by prepending "cd"
    std::ostringstream full;
    full << "cd " << workdir.string() << " && " << cmd;
    return std::system(full.str().c_str());
}

static void abaqus_init_case(const Config& cfg) 
{
    ensure_dir(cfg.case_dir);

    // Run init from inside case_dir (keeps BASE.* there)
    std::ostringstream cmd;
    cmd << cfg.abaqus_cmd
        << " cae noGUI=" << cfg.abaqus_script
        << " -- init"
        << " --case ."
        << " --seed " << cfg.seed;

    int rc = run_cmd(cmd.str(), cfg.case_dir);
    if (rc != 0) throw std::runtime_error("Abaqus init failed rc=" + std::to_string(rc));

    require_exists(cfg.case_dir / "DONE_INIT.txt", "Missing DONE_INIT.txt after init");
    require_exists(cfg.case_dir / "state.json", "Missing state.json after init");
    require_exists(cfg.case_dir / "mesh_BASE.vtk", "Missing mesh_BASE.vtk after init");
    require_exists(cfg.case_dir / "disp_BASE.csv", "Missing disp_BASE.csv after init");
}

static void abaqus_run_cut_step(const Config& cfg,
                                int step_idx,
                                const fs::path& cut_labels_file,
                                const std::string& job_name) 
{
    // Run cut from inside case_dir (python script will create step_### and copy labels).
    std::ostringstream cmd;
    cmd << cfg.abaqus_cmd
        << " cae noGUI=" << cfg.abaqus_script
        << " -- cut"
        << " --case ."
        << " --cut_labels " << cut_labels_file.string()
        << " --job " << job_name;

    int rc = run_cmd(cmd.str(), cfg.case_dir);
    if (rc != 0) throw std::runtime_error("Abaqus cut failed rc=" + std::to_string(rc));

    // The python creates case_dir/step_###/DONE.txt
    std::ostringstream step_folder;
    step_folder << "step_" << std::setw(3) << std::setfill('0') << step_idx;
    fs::path step_dir = cfg.case_dir / step_folder.str();
    require_exists(step_dir / "DONE.txt", "Missing DONE.txt after cut step");
}

// -------------------------------
// VTK mesh reader (minimal)
// -------------------------------
// You need:
// - points (node coords)
// - cells (connectivity)
// - cell_data elem_label (Abaqus element labels)
// - optionally cell types (ignore for centroid)

// Keep it simple: parse ASCII VTK UnstructuredGrid
struct VtkMesh 
{
    struct Point { double x,y,z; };
    std::vector<Point> points;
    std::vector<std::vector<int>> cells; // indices into points
    std::vector<int> elem_label;         // per cell
};

static VtkMesh load_vtk_unstructured_with_elem_label(const fs::path& vtk_path) 
{
    std::ifstream in(vtk_path);
    if (!in) throw std::runtime_error("Failed to open VTK: " + vtk_path.string());

    VtkMesh m;
    std::string tok;

    // very basic token parsing
    while (in >> tok) 
    {
        if (tok == "POINTS") 
        {
            int n; in >> n;
            std::string dtype; in >> dtype;
            m.points.resize(n);
            for (int i=0;i<n;i++) in >> m.points[i].x >> m.points[i].y >> m.points[i].z;
        } 
        else if (tok == "CELLS") 
        {
            int ncells, nints; in >> ncells >> nints;
            m.cells.resize(ncells);
            for (int c=0;c<ncells;c++) 
            {
                int k; in >> k;
                m.cells[c].resize(k);
                for (int j=0;j<k;j++) in >> m.cells[c][j];
            }

        } 
        else if (tok == "SCALARS") 
        {
            std::string name; in >> name;
            if (name == "elem_label") 
            {
                std::string type; in >> type;
                int ncomp; in >> ncomp; (void)ncomp;
                // consume LOOKUP_TABLE default
                in >> tok; // LOOKUP_TABLE
                in >> tok; // default
                // then one int per cell
                m.elem_label.resize(m.cells.size());
                for (size_t i=0;i<m.cells.size();i++) 
                {
                    in >> m.elem_label[i];
                }
            } 
            else 
            {
                // skip this scalar block crudely (not robust)
                // can extend if needed.
            }
        }
    }

    if (m.points.empty() || m.cells.empty() || m.elem_label.size() != m.cells.size()) 
    {
        throw std::runtime_error("VTK missing required sections (POINTS/CELLS/elem_label).");
    }
    return m;
}

static VtkMesh::Point cell_centroid(const VtkMesh& m, size_t ci) 
{
    const auto& c = m.cells[ci];
    double x=0,y=0,z=0;
    for (int vi : c) {
        x += m.points[vi].x;
        y += m.points[vi].y;
        z += m.points[vi].z;
    }
    double inv = 1.0 / double(c.size());
    return {x*inv, y*inv, z*inv};
}

// -------------------------------
// Cut volume classification (stub)
// -------------------------------
// Replace this with your CGAL implementation.
//
// Input:
// - current active FE mesh (VTK) with element labels
// - cut volume mesh (watertight) OR tool surface
//
// Output:
// - element labels to deactivate
static std::vector<int> compute_elements_to_remove(
    const VtkMesh& active_fe_mesh,
    const fs::path& cut_volume_mesh_path,
    int step_idx
) 
{
    (void)cut_volume_mesh_path;
    (void)step_idx;

    std::vector<int> labels;

    // ---- TODO: implement using CGAL ----
    // Pseudocode:
    // 1) load cut volume Surface_mesh from file
    // 2) create CGAL::Side_of_triangle_mesh for inside test
    // 3) for each cell in active_fe_mesh:
    //     centroid = cell_centroid(...)
    //     if inside(centroid): labels.push_back(elem_label[cell])
    //
    // Optional: add tolerance band or intersection rule.

    return labels;
}

// -------------------------------
// Write cut labels file
// -------------------------------
static fs::path write_cut_labels_file(const fs::path& step_dir,
                                      const std::vector<int>& elem_labels) 
{
    ensure_dir(step_dir);
    fs::path out = step_dir / "cut_elems.txt";
    std::ofstream f(out);
    if (!f) throw std::runtime_error("Failed to write cut labels: " + out.string());

    for (int lab : elem_labels) f << lab << "\n";
    return out;
}

// -------------------------------
// Your estimator hook (stub)
// -------------------------------
// This should:
// - read truth deformation/mesh
// - update your residual stress estimate
// - propose next cut geometry (and write out cut volume mesh)
static fs::path propose_next_cut_volume_mesh(
    const Config& cfg,
    int step_idx,
    const fs::path& current_truth_mesh_vtk,
    const fs::path& current_truth_disp_csv
) 
{
    (void)cfg; (void)step_idx; (void)current_truth_mesh_vtk; (void)current_truth_disp_csv;

    // ---- TODO: call your pipeline here ----
    // Return a mesh path of the cutter swept-volume (watertight) for this step.
    // Example output path:
    fs::path cut_mesh = cfg.case_dir / "cuts" / ("cut_volume_" + std::to_string(step_idx) + ".stl");
    ensure_dir(cut_mesh.parent_path());

    // For now, error out so you don't silently run nonsense
    throw std::runtime_error("propose_next_cut_volume_mesh: implement this (write cut volume mesh).");

    return cut_mesh;
}

// -------------------------------
// Main closed-loop driver
// -------------------------------
int main(int argc, char** argv) 
{
    (void)argc; (void)argv;

    CLI::App app{"Automated Abaqus Tester"};

    fs::path starting_mesh_stl, g_code_file;

    app.add_option("--starting_mesh_stl", starting_mesh_stl)->required();
    app.add_option("--g_code", g_code_file)->required();

    CLI11_PARSE(app, argc, argv);

    CylindricalTool tool{ .radius = 3.0, .height = 10.0 };

    std::vector<Line> lines;
    std::vector<ArcOfCircle> arcs;
    std::vector<InterpolatedCurve> splines;
    std::vector<Circle> circles;
    std::vector<ToolPath> toolpaths;


    std::ifstream g_code(g_code_file);
    std::string line_string;
    Position current_pos = {0, 0, 0};

    double cut_length = 0;
    while (std::getline(g_code, line_string))
    {
        std::stringstream line(line_string);

        std::string command;
        std::getline(line, command, ' ');

        // straight line cut
        if (command == "G0" || command == "G1") // line
        {
            std::string x_pos;
            std::string y_pos;
            std::string z_pos;
            std::getline(line, x_pos, ' ');
            std::getline(line, y_pos, ' ');
            std::getline(line, z_pos, ' ');

            // Remove the letter from the position value
            x_pos.erase(0, 1);
            y_pos.erase(0, 1);
            z_pos.erase(0, 1);

            Position next_pos = {std::stod(x_pos), std::stod(y_pos), std::stod(z_pos)};
            double line_length = std::sqrt((current_pos.x - next_pos.x) * (current_pos.x - next_pos.x) + (current_pos.y - next_pos.y) * (current_pos.y - next_pos.y));

            lines.push_back(Line{position_to_point3D(current_pos), delta_to_vec3D(current_pos, next_pos)});

            cut_length += line_length;
            current_pos = next_pos;

        }

        // clock wise arc
        else if (command == "G2")
        {
            std::string end_pos_x, end_pos_y, x_offset, y_offset;
            std::getline(line, end_pos_x, ' ');
            std::getline(line, end_pos_y, ' ');
            std::getline(line, x_offset, ' ');
            std::getline(line, y_offset, ' ');

            end_pos_x.erase(0, 1);
            end_pos_y.erase(0, 1);
            y_offset.erase(0, 2);

            Position end_pos = {stod(end_pos_x), stod(end_pos_y), 0};

            Position arc_center = current_pos - Position({stod(x_offset), stod(y_offset), 0});

            std::pair<Position, double> arc_info = arcMidpoint(arc_center, current_pos, end_pos, false);

            Position mid_point = arc_info.first;
            double arc_length = arc_info.second;

            //ArcOfCircle::ArcOfCircle(const std::pair<Point3D, Point3D>& arc_endpoints, const Point3D& arc_interior_point)
            arcs.push_back(ArcOfCircle({position_to_point3D(current_pos), position_to_point3D(end_pos)}, 
                                       position_to_point3D(mid_point)));

            cut_length += arc_length;
            current_pos = end_pos;
        }

        // conuter clock wise arc
        else if (command == "G3")
        {
            std::string end_pos_x, end_pos_y, x_offset, y_offset;
            std::getline(line, end_pos_x, ' ');
            std::getline(line, end_pos_y, ' ');
            std::getline(line, x_offset, ' ');
            std::getline(line, y_offset, ' ');

            end_pos_x.erase(0, 1);
            end_pos_y.erase(0, 1);
            y_offset.erase(0, 2);

            Position end_pos = {stod(end_pos_x), stod(end_pos_y), 0};

            Position arc_center = current_pos - Position({stod(x_offset), stod(y_offset), 0});

            std::pair<Position, double> arc_info = arcMidpoint(arc_center, current_pos, end_pos, true);

            Position mid_point = arc_info.first;
            double arc_length = arc_info.second;

            //ArcOfCircle::ArcOfCircle(const std::pair<Point3D, Point3D>& arc_endpoints, const Point3D& arc_interior_point)
            arcs.push_back(ArcOfCircle({position_to_point3D(current_pos), position_to_point3D(end_pos)}, 
                                       position_to_point3D(mid_point)));

            cut_length += arc_length;
            current_pos = end_pos;
        }

        if (cut_length > LENGTH_OF_CUT)
        {
            //ToolPath::ToolPath(const std::tuple<std::vector<Line>, 
            //                        std::vector<ArcOfCircle>,
            //                        std::vector<InterpolatedCurve>,
            //                        std::vector<Circle>> compound,
            //       const CylindricalTool& profile,
            //       const bool display)

            auto compound = std::make_tuple(lines, arcs, splines, circles);

            ToolPath toolpath(compound, tool, false);

            toolpaths.push_back(toolpath);

            // triangulate it (controls triangle quality)
            toolpath.mesh_surface(/*angle=*/0.35, /*deflection=*/0.1);

            // write swept-volume STL (this is your “cut volume”)
            toolpath.shape_to_stl("cut_step_001", "/abs/path/cut_step_001.stl");

            cut_length = 0;
        }

    }


    Config cfg;

    // 1) Init truth baseline
    abaqus_init_case(cfg);

    // current truth after BASE
    fs::path current_mesh = cfg.case_dir / "mesh_BASE.vtk";
    fs::path current_disp = cfg.case_dir / "disp_BASE.csv";

    // 2) Iterative loop
    for (int step = 1; step <= cfg.max_steps; ++step) {
        std::cout << "\n=== STEP " << step << " ===\n";

        // a) Your code proposes next cut volume mesh (watertight)
        fs::path cut_volume_mesh = propose_next_cut_volume_mesh(cfg, step, current_mesh, current_disp);

        // b) Load current active FE mesh (truth mesh)
        VtkMesh fe = load_vtk_unstructured_with_elem_label(current_mesh);

        // c) Classify which elements to remove based on cut volume
        std::vector<int> remove_labels = compute_elements_to_remove(fe, cut_volume_mesh, step);
        if (remove_labels.empty()) 
        {
            std::cout << "No elements selected for removal; stopping.\n";
            break;
        }

        // d) Write cut labels file into step folder (for record)
        std::ostringstream step_name;
        step_name << "step_" << std::setw(3) << std::setfill('0') << step;
        fs::path step_dir = cfg.case_dir / step_name.str();

        fs::path cut_labels_file = write_cut_labels_file(step_dir, remove_labels);

        // e) Call Abaqus to run ONE truth cut job
        std::ostringstream job;
        job << "CUT_" << std::setw(3) << std::setfill('0') << step;
        abaqus_run_cut_step(cfg, step, cut_labels_file, job.str());

        // f) Update current truth file paths
        current_mesh = step_dir / ("mesh_" + job.str() + ".vtk");
        current_disp = step_dir / ("disp_" + job.str() + ".csv");

        require_exists(current_mesh, "Missing truth mesh after cut");
        require_exists(current_disp, "Missing truth disp after cut");

        // g) Stopping rule hooks (optional):
        // - compute remaining volume from mesh
        // - stop if volume < threshold or if close to target geometry
    }

    std::cout << "\nDone.\n";
    return 0;
}

