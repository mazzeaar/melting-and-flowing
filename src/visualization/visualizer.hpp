#pragma once

#include <span>

// forward declaration
struct Simulation;

namespace visualization {

struct Visualizer {
    explicit Visualizer(Simulation& sim);

    Visualizer(const Visualizer&) = delete;
    Visualizer& operator=(const Visualizer&) = delete;

    void run();

   private:
    Simulation& sim_;

    bool running_ = false;
    double dt_ = 0.01;
    int steps_per_frame_ = 1;

    void initialize();
    void update();
    void step();
    void callback();
    void draw_controls();
    void load_mesh(std::span<const std::string> files);
};

}  // namespace visualization
