#include <array>

#include "simulation.hpp"
#include "visualization/visualizer.hpp"

int main() {
    Simulation sim{
        std::array{20, 20, 20},
        std::array{20.0, 20.0, 20.0},
        std::array{0.0, 0.0, 0.0},
    };

    visualization::Visualizer viewer{sim};
    viewer.run();
}
