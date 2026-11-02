cmake_minimum_required(VERSION 3.28)
project(fluid_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Dependencies
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

add_executable(fluid_sim ${SOURCES})

target_compile_options(
  fluid_sim
  PRIVATE # surpress libigl warnings
          $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-invalid-offsetof>)

target_include_directories(fluid_sim PRIVATE "${PROJECT_SOURCE_DIR}/src")

target_link_libraries(fluid_sim PRIVATE Eigen3::Eigen polyscope igl::core)

target_compile_definitions(
  fluid_sim PRIVATE FLUID_ASSET_DIR="${PROJECT_SOURCE_DIR}/assets")
