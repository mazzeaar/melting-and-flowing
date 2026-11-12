#include "simulation/simulation.hpp"
#include "visualization/visualizer.hpp"

int main() {
    auto settings = SimulationSettings{};

    Simulation sim{settings};

    visualization::Visualizer viewer{sim};
    viewer.run();
}
