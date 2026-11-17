add_subdirectory(${PROJECT_SOURCE_DIR}/external/tomlplusplus
                 ${CMAKE_BINARY_DIR}/tomlplusplus EXCLUDE_FROM_ALL)

add_library(fluid_core_benchmark STATIC ${SOURCES})

target_include_directories(fluid_core_benchmark
                           PUBLIC "${PROJECT_SOURCE_DIR}/src")

target_link_libraries(fluid_core_benchmark PUBLIC Eigen3::Eigen polyscope
                                                  igl::core OpenMP::OpenMP_CXX)

target_compile_definitions(
  fluid_core_benchmark PUBLIC ENABLE_BENCHMARKS
                              FLUID_ASSET_DIR="${PROJECT_SOURCE_DIR}/assets")

target_compile_options(
  fluid_core_benchmark
  PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-invalid-offsetof>)

add_executable(fluid_benchmark main.cpp)

target_link_libraries(fluid_benchmark PRIVATE fluid_core_benchmark
                                              tomlplusplus::tomlplusplus)

target_compile_options(
  fluid_benchmark
  PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-invalid-offsetof>)
