cmake_minimum_required(VERSION 3.15)
project(residual_stress_estimator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# -----------------------------
# Dependencies
# -----------------------------

# CLI11 (via FetchContent)
include(FetchContent)
FetchContent_Declare(
    CLI11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.4.1
)
FetchContent_MakeAvailable(CLI11)

# CGAL (system-installed)
find_package(CGAL REQUIRED)

# NLOPT (system-installed)
find_library(NLOPT_LIB nlopt REQUIRED)

# -----------------------------
# Source files
# -----------------------------
set(SRC_FILES
    src/mesh_handler.cpp
    src/fea_solver.cpp
    src/optimizer.cpp
    src/stress_field.cpp
)

set(HDR_FILES
    src/mesh_handler.hpp
    src/fea_solver.hpp
    src/optimizer.hpp
    src/stress_field.hpp
)

# -----------------------------
# Main executable
# -----------------------------
add_executable(residual_stress_estimator
    ${SRC_FILES}
    src/main.cpp
)

target_compile_options(residual_stress_estimator PRIVATE -frounding-math)

target_link_libraries(residual_stress_estimator
    PRIVATE
        CLI11::CLI11
        CGAL::CGAL
        ${NLOPT_LIB}
        gmp
        mpfr
)

# -----------------------------
# Unit test executable
# -----------------------------
add_executable(mesh_handler_tester
    ${SRC_FILES}
    src/mesh_handler_tester.cpp
)

target_link_libraries(mesh_handler_tester
    PRIVATE
        CLI11::CLI11
        CGAL::CGAL
        ${NLOPT_LIB}
        gmp
        mpfr
)

# Optional: organize binaries in bin/
set_target_properties(residual_stress_estimator mesh_handler_tester
    PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
