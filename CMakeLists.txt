cmake_minimum_required(VERSION 3.15)
project(residual_stress_estimator LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)

FetchContent_Declare(
  CLI11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG v2.4.1
)
FetchContent_MakeAvailable(CLI11)

find_package(CGAL REQUIRED)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(NLopt REQUIRED)
find_package(glfw3 REQUIRED)


# Only core OCCT components
find_package(OpenCASCADE REQUIRED COMPONENTS
  FoundationClasses
  ModelingData
  ModelingAlgorithms
)

set(TOOLPATH_DIR "$ENV{HOME}/code/research/residual_stress/surfacic_toolpaths")

file(GLOB_RECURSE TOOLPATH_SOURCES CONFIGURE_DEPENDS
  "${TOOLPATH_DIR}/src/*.cpp"
)

set(SRC_FILES
  src/mesh_handler.cpp
  src/fea_solver.cpp
  ${TOOLPATH_SOURCES}
)

function(link_common target)
  target_compile_options(${target} PRIVATE -frounding-math)

  target_include_directories(${target} PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${TOOLPATH_DIR}/include
    ${TOOLPATH_DIR}/src/include
    ${GMP_INCLUDE_DIR}
    ${MPFR_INCLUDE_DIR}
    ${OpenCASCADE_INCLUDE_DIR}
  )

  target_link_libraries(${target} PRIVATE
    CLI11::CLI11
    CGAL::CGAL
    NLopt::nlopt
    glfw
    ${GMP_LIBRARIES}
    ${MPFR_LIBRARIES}
    ${OpenCASCADE_LIBRARIES}
  )
endfunction()

#add_executable(residual_stress_estimator
#  ${SRC_FILES}
#  src/main.cpp
#)
#link_common(residual_stress_estimator)
#
#add_executable(mesh_handler_tester
#  ${SRC_FILES}
#  src/mesh_handler_tester.cpp
#)
#link_common(mesh_handler_tester)

#add_executable(set_stress_test
#  ${SRC_FILES}
#  src/set_stress_test.cpp
#)
#link_common(set_stress_test)

#add_executable(absolute_truth_tester
#  ${SRC_FILES}
#  src/tester.cpp
#)
#link_common(absolute_truth_tester)

add_executable(abaqus_tester
  ${SRC_FILES}
  src/abaqus_test.cpp
)
link_common(abaqus_tester)

set_target_properties(
  #residual_stress_estimator
  #mesh_handler_tester
  #set_stress_test
  #absolute_truth_tester
  abaqus_tester
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
