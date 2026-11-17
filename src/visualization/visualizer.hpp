#pragma once

#include "simulation/mesh.hpp"
#include "simulation/simulation.hpp"

// forward declaration
struct Simulation;

namespace visualization {

struct Visualizer {
    explicit Visualizer(Simulation& sim);
    ~Visualizer();

    Visualizer(const Visualizer&) = delete;
    Visualizer& operator=(const Visualizer&) = delete;

    void run();

   private:
    Simulation& sim_;
    std::unique_ptr<Mesh> mesh_ptr;

    bool running_ = false;
    int steps_per_frame_ = 1;

    void init();
    void update();
    void reset();
    void step();
    void callback();
    void draw_controls();
    void load_mesh(std::filesystem::path);
    void save_mesh();
    void update_particles();
    void update_grid(bool rebuild = false);
};

}  // namespace visualization
