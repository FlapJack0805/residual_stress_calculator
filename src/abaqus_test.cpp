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
#include <cmath>
#include <iomanip>
#include <optional>
#include <cctype>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>

#define LENGTH_OF_CUT 10 // TODO: Find out a reasonable size for each cut

namespace fs = std::filesystem;

struct GCode
{
    std::optional<int> g;   // 0,1,2,3 if present
    std::optional<double> x;
    std::optional<double> y;
    std::optional<double> z;
    std::optional<double> i;
    std::optional<double> j;
};



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

    Position operator+(const Position& other) const
    {
        return {x + other.x, y + other.y , z + other.z};
    }

    Position operator-(const Position& other) const
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
// We need:
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

static Point cell_centroid(const VtkMesh& m, size_t ci) 
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
// Write cut labels file
// -------------------------------
static fs::path write_cut_labels_file(const fs::path& step_dir,
                                      const Surface_mesh &cut_mesh,
                                      VtkMesh &vtk_mesh) 
{
    ensure_dir(step_dir);
    fs::path out = step_dir / "cut_elems.txt";
    std::ofstream f(out);
    if (!f) throw std::runtime_error("Failed to write cut labels: " + out.string());

    Tree cut_tree(faces(cut_mesh).begin(), faces(cut_mesh).end(), cut_mesh);
    cut_tree.accelerate_distance_queries();

    CGAL::Side_of_triangle_mesh<Surface_mesh, Kernel> inside_cut_test(cut_tree);

    for (size_t i = 0; i < vtk_mesh.elem_label.size(); ++i) 
    {
        Point centroid = cell_centroid(vtk_mesh, i);
        if (inside_cut_test(centroid) == CGAL::ON_BOUNDED_SIDE)// all elements inside of the cut get added to the cut file
        {
            f << vtk_mesh.elem_label[i] << "\n";
        }
    }
    return out;
}



// -------------------------------
// Your estimator hook (stub)
// -------------------------------
// This should:
// - read truth deformation/mesh
// - update your residual stress estimate
// - propose next cut geometry (and write out cut volume mesh)
static Surface_mesh propose_next_cut_volume_mesh(
    const Config& cfg,
    int step_idx,
    const fs::path& current_truth_mesh_vtk,
    const fs::path& current_truth_disp_csv,
    MeshHandler &mesh_handler
) 
{
    (void)cfg; (void)step_idx; (void)current_truth_mesh_vtk; (void)current_truth_disp_csv;


    Surface_mesh cut_mesh = mesh_handler.make_best_cut(); // TODO: Make the make_cut() function and make sure it returns the cut it's making as a surface mesh

    return cut_mesh;
}

static std::string trim(const std::string& s)
{
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

static std::string strip_comments(const std::string& line)
{
    std::string out = line;

    // Remove semicolon comments
    auto semi = out.find(';');
    if (semi != std::string::npos) {
        out = out.substr(0, semi);
    }

    // Remove (...) comments
    std::string cleaned;
    bool in_paren = false;
    for (char c : out) {
        if (c == '(') {
            in_paren = true;
            continue;
        }
        if (c == ')') {
            in_paren = false;
            continue;
        }
        if (!in_paren) {
            cleaned.push_back(c);
        }
    }

    return trim(cleaned);
}

static bool parse_word_value(const std::string& token, char expected_prefix, double& value)
{
    if (token.size() < 2) return false;
    if (std::toupper(static_cast<unsigned char>(token[0])) != expected_prefix) return false;

    try {
        value = std::stod(token.substr(1));
        return true;
    } catch (...) {
        return false;
    }
}

static bool parse_word_int(const std::string& token, char expected_prefix, int& value)
{
    if (token.size() < 2) return false;
    if (std::toupper(static_cast<unsigned char>(token[0])) != expected_prefix) return false;

    try {
        value = std::stoi(token.substr(1));
        return true;
    } catch (...) {
        return false;
    }
}

static GCode parse_gcode_line(const std::string& raw_line)
{
    GCode move;
    std::string line = strip_comments(raw_line);
    if (line.empty()) return move;

    std::stringstream ss(line);
    std::string token;

    while (ss >> token) {
        if (token.empty()) continue;

        // normalize token prefix to uppercase
        token[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(token[0])));

        int gcode_int;
        double value;


        if (parse_word_int(token, 'G', gcode_int)) {
            move.g = gcode_int;
        } else if (parse_word_value(token, 'X', value)) {
            move.x = value;
        } else if (parse_word_value(token, 'Y', value)) {
            move.y = value;
        } else if (parse_word_value(token, 'Z', value)) {
            move.z = value;
        } else if (parse_word_value(token, 'I', value)) {
            move.i = value;
        } else if (parse_word_value(token, 'J', value)) {
            move.j = value;
        }
    }

    return move;
}

static void flush_toolpath_chunk(
    const Config& cfg,
    const CylindricalTool& tool,
    std::vector<Line>& lines,
    std::vector<ArcOfCircle>& arcs,
    std::vector<InterpolatedCurve>& splines,
    std::vector<Circle>& circles,
    std::vector<Surface_mesh>& cut_meshes)
{
    if (lines.empty() && arcs.empty() && splines.empty() && circles.empty()) {
        return;
    }

    auto compound = std::make_tuple(lines, arcs, splines, circles);
    ToolPath tp(compound, tool, false);

    tp.mesh_surface(0.35, 0.1);

    std::ostringstream name;
    name << "cut_step_" << std::setw(3) << std::setfill('0') << (cut_meshes.size() + 1);

    fs::path cut_stl = cfg.case_dir / "cuts" / (name.str() + ".stl");
    ensure_dir(cut_stl.parent_path());
    tp.shape_to_stl(name.str(), cut_stl.string());


    Surface_mesh cut_mesh;
    // load the cut stl into a surface mesh
    if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(cut_stl.string(), cut_mesh)) {
        throw std::runtime_error("Failed to load mesh: " + cut_stl.string());
    }

    cut_meshes.push_back(std::move(cut_mesh));

    lines.clear();
    arcs.clear();
    splines.clear();
    circles.clear();
}

static double xy_distance(const Position& a, const Position& b)
{
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    return std::sqrt(dx * dx + dy * dy);
}



void convert_stl_to_vtk_python(const std::string& stl_path)
{
    std::string output = "starting_mesh_vtk.vtk";

    std::string python = "/home/jstifter/code/research/residual_stress/residual_stress_estimator/src/venv/bin/python";

    std::string cmd = python + " /home/jstifter/code/research/residual_stress/residual_stress_estimator/src/convert_stl_to_vtk.py \"" +
                      stl_path + "\" \"" +
                      output + "\"";

    std::cout << "Running: " << cmd << std::endl;

    int rc = std::system(cmd.c_str());

    if (rc != 0)
    {
        throw std::runtime_error("STL → VTK conversion failed");
    }

    std::cout << "Converted STL to " << output << std::endl;
}

// -------------------------------
// Main closed-loop driver
// -------------------------------
int main(int argc, char** argv) 
{
    (void)argc; (void)argv;

    CLI::App app{"Automated Abaqus Tester"};

    fs::path starting_mesh_stl, boundary_mesh_stl, starting_mesh_vtk, g_code_file;

    app.add_option("--starting_mesh_stl", starting_mesh_stl)->required();
    app.add_option("--boundary_mesh_stl", boundary_mesh_stl)->required();
    app.add_option("--g_code", g_code_file)->required();

    CLI11_PARSE(app, argc, argv);

    convert_stl_to_vtk_python(starting_mesh_stl);

    Config cfg;

    CylindricalTool tool{ .radius = 1.0, .height = 10.0 };

    std::vector<Line> lines;
    std::vector<ArcOfCircle> arcs;
    std::vector<InterpolatedCurve> splines;
    std::vector<Circle> circles;
    std::vector<Surface_mesh> cut_meshes;

    std::ifstream g_code(g_code_file);
    if (!g_code) {
        throw std::runtime_error("Failed to open g-code file: " + g_code_file.string());
    }

    std::string line_string;
    Position current_pos{0.0, 0.0, 0.0};
    std::optional<int> modal_g;
    double cut_length = 0.0;

    while (std::getline(g_code, line_string))
    {
        std::cout << line_string << std::endl;
        GCode move = parse_gcode_line(line_string);

        // Support modal G-codes: if a line omits G, reuse last one
        if (move.g.has_value()) {
            modal_g = move.g;
        }

        if (!modal_g.has_value()) {
            continue;
        }

        const int g = *modal_g;

        // Only care about motion commands we support
        if (g != 0 && g != 1 && g != 2 && g != 3) {
            continue;
        }

        Position next_pos = current_pos;

        if (move.x.has_value()) next_pos.x = *move.x;
        if (move.y.has_value()) next_pos.y = *move.y;
        if (move.z.has_value()) next_pos.z = *move.z;

        // Ignore zero-length motion
        const bool changed =
            (next_pos.x != current_pos.x) ||
            (next_pos.y != current_pos.y) ||
            (next_pos.z != current_pos.z);

        if (!changed && g != 2 && g != 3) {
            continue;
        }

        // G0 = rapid move, usually NOT cutting
        if (g == 0) {
            current_pos = next_pos;
            continue;
        }

        // G1 = linear cutting move
        if (g == 1)
        {
            const bool xy_changed =
                (next_pos.x != current_pos.x) ||
                (next_pos.y != current_pos.y);

            // Skip pure-Z plunges / retracts for now
            if (!xy_changed) {
                current_pos = next_pos;
                continue;
            }

            const double line_length = xy_distance(current_pos, next_pos);

            if (line_length > 0.0) {
                lines.push_back(Line{
                    position_to_point3D(current_pos),
                    delta_to_vec3D(current_pos, next_pos)
                });

                cut_length += line_length;
            }

            current_pos = next_pos;
        }

        // G2/G3 = circular arcs in XY plane using I/J center offsets
        else if (g == 2 || g == 3)
        {
            if (!move.i.has_value() || !move.j.has_value()) {
                throw std::runtime_error("Arc move missing I/J offsets: " + line_string);
            }

            Position arc_center{
                current_pos.x + *move.i,
                current_pos.y + *move.j,
                current_pos.z
            };

            const bool ccw = (g == 3);
            auto [mid_point, arc_length] = arcMidpoint(arc_center, current_pos, next_pos, ccw);

            arcs.emplace_back(
                std::make_pair(position_to_point3D(current_pos),
                               position_to_point3D(next_pos)),
                position_to_point3D(mid_point)
            );

            cut_length += arc_length;
            current_pos = next_pos;
        }

        if (cut_length >= LENGTH_OF_CUT)
        {
            flush_toolpath_chunk(cfg, tool, lines, arcs, splines, circles, cut_meshes);
            cut_length = 0.0;
        }
    }

    // Flush final partial chunk if we have one
    // NOTE: With how we measure our best cuts at the moment this smaller cut will probably be seen as optimal simply because it is small.
    // I don't know if Biruk's code accounts for size when finding the best cut but it should be something to look into. I also don't think it's a big deal
    // though because it's only 1 cut, everything else is of uniform size.
    if (cut_length > 0)
    {
        flush_toolpath_chunk(cfg, tool, lines, arcs, splines, circles, cut_meshes);
    }

    // 1) Init truth baseline
    abaqus_init_case(cfg);
    MeshHandler mesh_handler(starting_mesh_stl, boundary_mesh_stl, cut_meshes);
    VtkMesh vtk_mesh = load_vtk_unstructured_with_elem_label(starting_mesh_vtk);

    // current truth after BASE
    fs::path current_mesh = cfg.case_dir / "mesh_BASE.vtk";
    fs::path current_disp = cfg.case_dir / "disp_BASE.csv";

    // 2) Iterative loop
    for (int step = 1; step <= cfg.max_steps; ++step) {
        std::cout << "\n=== STEP " << step << " ===\n";

        // Stop if we ran out of precomputed cut chunks
        if (step > static_cast<int>(cut_meshes.size())) {
            std::cout << "No more cut chunks available. Stopping.\n";
            break;
        }

        // a) Get next cut volume mesh
        Surface_mesh cut_mesh = cut_meshes[step - 1];

        // Optional debug export of the cut mesh
        fs::path cut_mesh_stl = cfg.case_dir / "cuts" / ("cut_volume_" + std::to_string(step) + ".stl");
        ensure_dir(cut_mesh_stl.parent_path());
        CGAL::IO::write_polygon_mesh(cut_mesh_stl.string(), cut_mesh);

        // b) Load current active FE mesh (truth mesh)
        VtkMesh vtk_mesh = load_vtk_unstructured_with_elem_label(current_mesh);

        // d) Write cut labels file into step folder (this also does c: classification)
        std::ostringstream step_name;
        step_name << "step_" << std::setw(3) << std::setfill('0') << step;
        fs::path step_dir = cfg.case_dir / step_name.str();

        fs::path cut_labels_file = write_cut_labels_file(step_dir, cut_mesh, vtk_mesh);

        // e) Call Abaqus to run ONE truth cut job
        std::ostringstream job;
        job << "CUT_" << std::setw(3) << std::setfill('0') << step;
        abaqus_run_cut_step(cfg, step, cut_labels_file, job.str());

        // f) Update current truth file paths
        current_mesh = step_dir / ("mesh_" + job.str() + ".vtk");
        current_disp = step_dir / ("disp_" + job.str() + ".csv");

        require_exists(current_mesh, "Missing truth mesh after cut");
        require_exists(current_disp, "Missing truth disp after cut");
    }

    std::cout << "\nDone.\n";
    return 0;
}

