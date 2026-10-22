cmake_minimum_required(VERSION 3.28)
project(fluid_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(external/eigen EXCLUDE_FROM_ALL)
add_subdirectory(external/polyscope)

add_executable(fluid_sim main.cpp)
target_link_libraries(fluid_sim PRIVATE Eigen3::Eigen polyscope)
