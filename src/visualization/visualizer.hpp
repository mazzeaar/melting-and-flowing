#pragma once

#include "grid_stats_panel.hpp"
#include "runners.hpp"

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
    // UI settings are separate from the worker-owned live simulation.
    SimulationSettings settings_;
    std::filesystem::path mesh_path_;
    std::size_t generation_ = 0;
    std::shared_ptr<const SimulationSnapshot> snapshot_;
    std::shared_ptr<const MeshSnapshot> mesh_snapshot_;
    SimulationRunner simulation_;
    MeshBuilderRunner mesh_builder_;

    GridStatisticsPanel grid_stats_;

    int resolution_upfactor = 2;
    float radius_factor = 1.5f;
    int marker_stride = 1;

    bool running_ = false;

    void init();
    void update();
    void reset();
    void callback();
    void draw_controls();
    void load_mesh(std::filesystem::path);
    void update_particles();
    void update_grid(bool rebuild = false);
    void update_mesh();
};

}  // namespace visualization
