cmake_minimum_required(VERSION 3.28)
project(fluid_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Dependencies
find_package(OpenMP REQUIRED)
add_subdirectory(external/eigen EXCLUDE_FROM_ALL)
add_subdirectory(external/polyscope)

# libigl: core geometry functionality only
set(LIBIGL_BUILD_TESTS
    OFF
    CACHE BOOL "" FORCE)
set(LIBIGL_BUILD_TUTORIALS
    OFF
    CACHE BOOL "" FORCE)
set(LIBIGL_INSTALL
    OFF
    CACHE BOOL "" FORCE)

set(LIBIGL_GLFW
    OFF
    CACHE BOOL "" FORCE)
set(LIBIGL_IMGUI
    OFF
    CACHE BOOL "" FORCE)
set(LIBIGL_OPENGL
    OFF
    CACHE BOOL "" FORCE)

set(LIBIGL_EMBREE
    OFF
    CACHE BOOL "" FORCE)
set(LIBIGL_COPYLEFT_CGAL
    OFF
    CACHE BOOL "" FORCE)
set(LIBIGL_COPYLEFT_COMISO
    OFF
    CACHE BOOL "" FORCE)
set(LIBIGL_COPYLEFT_TETGEN
    OFF
    CACHE BOOL "" FORCE)
set(LIBIGL_RESTRICTED_TRIANGLE
    OFF
    CACHE BOOL "" FORCE)

add_subdirectory(external/libigl EXCLUDE_FROM_ALL)

# Sources
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM SOURCES "${PROJECT_SOURCE_DIR}/src/main.cpp")

# Core library
add_library(fluid_core STATIC ${SOURCES})

target_include_directories(fluid_core PUBLIC "${PROJECT_SOURCE_DIR}/src")

target_link_libraries(
  fluid_core
  PUBLIC Eigen3::Eigen polyscope igl::core
  PRIVATE OpenMP::OpenMP_CXX)

# Suppress warnings originating from libigl
target_compile_options(
  fluid_core
  PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-invalid-offsetof>)

# Main executable
add_executable(fluid_sim src/main.cpp)
target_compile_options(
  fluid_sim
  PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-invalid-offsetof>)
target_link_libraries(fluid_sim PRIVATE fluid_core)

target_compile_definitions(
  fluid_core PUBLIC FLUID_ASSET_DIR="${PROJECT_SOURCE_DIR}/assets")

# Tests
enable_testing()

add_executable(test_pressure tests/test_pressure.cpp)
target_link_libraries(test_pressure PRIVATE fluid_core)

add_test(NAME pressure COMMAND test_pressure)

# add benchmarks
add_subdirectory(benchmark)
