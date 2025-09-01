cmake_minimum_required(VERSION 3.14)
project(residual_stress_estimator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -----------------------------
# Dependencies
# -----------------------------

include(FetchContent)
FetchContent_Declare(
    CLI11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.4.1
)
FetchContent_MakeAvailable(CLI11)

find_library(NLOPT_LIB nlopt REQUIRED)

# -----------------------------
# Executable
# -----------------------------

add_executable(residual_stress_estimator src/force_estimator.cpp src/main.cpp)


target_compile_options(residual_stress_estimator PRIVATE -frounding-math)


# -----------------------------
# Link Libraries
# -----------------------------

target_link_libraries(residual_stress_estimator
    PRIVATE
        CLI11::CLI11
        gmp
        mpfr
        ${NLOPT_LIB}
)
