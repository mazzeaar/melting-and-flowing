cmake_minimum_required(VERSION 3.28)
project(fluid_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(external/eigen EXCLUDE_FROM_ALL)
add_subdirectory(external/polyscope)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM SOURCES "${PROJECT_SOURCE_DIR}/src/main.cpp")

add_library(fluid_core STATIC ${SOURCES})

target_include_directories(fluid_core PUBLIC "${PROJECT_SOURCE_DIR}/src")

target_link_libraries(fluid_core PUBLIC Eigen3::Eigen)

add_executable(fluid_sim src/main.cpp)

target_link_libraries(fluid_sim PRIVATE fluid_core polyscope)

target_compile_definitions(
  fluid_sim PRIVATE FLUID_ASSET_DIR="${PROJECT_SOURCE_DIR}/assets")

enable_testing()

add_executable(test_pressure tests/test_pressure.cpp)

target_link_libraries(test_pressure PRIVATE fluid_core)

add_test(NAME pressure COMMAND test_pressure)
