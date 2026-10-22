cmake_minimum_required(VERSION 3.28)
project(fluid_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(external/eigen EXCLUDE_FROM_ALL)
add_subdirectory(external/polyscope)

add_executable(fluid_sim src/main.cpp)

target_include_directories(fluid_sim PRIVATE "${PROJECT_SOURCE_DIR}/src")

target_link_libraries(fluid_sim PRIVATE Eigen3::Eigen polyscope)

target_compile_definitions(
  fluid_sim PRIVATE FLUID_ASSET_DIR="${PROJECT_SOURCE_DIR}/assets")
