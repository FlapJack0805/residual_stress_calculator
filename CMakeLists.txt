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
find_package(NLopt REQUIRED)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

# NLOPT (system-installed)
find_library(NLOPT_LIB nlopt REQUIRED)

# -----------------------------
# Source files
# -----------------------------
set(SRC_FILES
    src/mesh_handler.cpp
    src/fea_solver.cpp
    ~/code/research/residual_stress/surfacic_toolpaths/src/*.cpp #Change this to your path
)

set(HDR_FILES
    src/mesh_handler.hpp
    src/fea_solver.hpp
    ~/code/research/residual_stress/surfacic_toolpaths/include/toolpath.hxx #Change this to your path
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
# Unit test executables
# -----------------------------
add_executable(mesh_handler_tester
    ${SRC_FILES}
    src/mesh_handler_tester.cpp
)

add_executable(set_stress_test
    ${SRC_FILES}
    src/set_stress_test.cpp
)

add_executable(absolute_truth_tester
    ${SRC_FILES}
    src/tester.cpp
)


target_link_libraries(set_stress_test
    PRIVATE
        CLI11::CLI11
        CGAL::CGAL
        ${NLOPT_LIB}
        gmp
        mpfr
)

target_link_libraries(mesh_handler_tester
    PRIVATE
        CLI11::CLI11
        CGAL::CGAL
        ${NLOPT_LIB}
        gmp
        mpfr
)

target_link_libraries(absolute_truth_tester
    PRIVATE
        CLI11::CLI11
        CGAL::CGAL
        ${NLOPT_LIB}
        gmp
        mpfr
)

# Optional: organize binaries in bin/
set_target_properties(residual_stress_estimator mesh_handler_tester set_stress_test
    PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
