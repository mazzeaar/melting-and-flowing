#include <polyscope/point_cloud.h>
#include <polyscope/polyscope.h>

#include <Eigen/Dense>
#include <array>
#include <core/paths.hpp>
#include <grid/example.hpp>
#include <print>
#include <unordered_map>
#include <vector>

int main() {
  grid::example::func();

  // Eigen: basic vector operations
  Eigen::Vector3d position{1.0, 2.0, 3.0};
  Eigen::Vector3d direction{0.0, 1.0, 0.0};

  std::println("Position: {}", position.transpose());
  std::println("Length: {}", position.norm());
  std::println("Dot product: {}", position.dot(direction));

  // Polyscope: visualize a point cloud
  std::vector<std::array<double, 3>> points{
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}, {1.0, 1.0, 1.0},
  };

  polyscope::init();
  polyscope::registerPointCloud("Points", points);
  polyscope::show();
}
