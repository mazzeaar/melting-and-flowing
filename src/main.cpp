#include <array>

#include "simulation.hpp"
#include "visualization/visualizer.hpp"

int main() {
    auto n = 15;
    auto extents = std::array{n, n, n};
    auto size = std::array{1.0, 1.0, 1.0};
    auto anchor = std::array{0.0, 0.0, 0.0};

    Simulation sim{
        extents,
        size,
        anchor,
    };

    visualization::Visualizer viewer{sim};
    viewer.run();
}
