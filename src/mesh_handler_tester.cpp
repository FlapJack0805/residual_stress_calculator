#include "mesh_handler.hpp"
#include <CLI/CLI.hpp>
#include <iostream>

int main(int argc, char** argv)
{
    CLI::App app{"MeshHandler Tester"};
    std::string mesh_stl, boundary_stl, tool_path_stl;

    app.add_option("--mesh", mesh_stl)->required();
    app.add_option("--boundary", boundary_stl)->required();
    app.add_option("--toolpath", tool_path_stl);

    CLI11_PARSE(app, argc, argv);

    MeshHandler handler(mesh_stl, boundary_stl);

    std::cout << "Before tool path applied\n";
    handler.print_summary();

    if (!tool_path_stl.empty()) {
        handler.set_tool_path(tool_path_stl);
        std::cout << "After tool path applied\n";
        handler.print_summary(); // shows changes after toolpath applied
    }
}
