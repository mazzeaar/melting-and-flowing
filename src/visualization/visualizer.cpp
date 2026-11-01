#include "visualizer.hpp"

#include <imgui.h>
#include <polyscope/polyscope.h>

#include "polyscope/messages.h"
#include "simulation.hpp"

namespace visualization {

Visualizer::Visualizer(Simulation& sim) : sim_{sim} {}

void Visualizer::run() {
    initialize();

    polyscope::state::userCallback = [this] { callback(); };
    polyscope::state::filesDroppedCallback = [this](auto files) {
        load_mesh(files);
    };

    polyscope::show();

    polyscope::state::filesDroppedCallback = nullptr;
    polyscope::state::userCallback = nullptr;
}

void Visualizer::load_mesh(std::span<const std::string> files) {
    for (std::filesystem::path file : files) {
        if (file.extension() != ".obj") {
            polyscope::warning(
                "Invalid file format, only *.obj is ok",
                std::format("file extension '{}' is not supported",
                            file.extension().c_str()));
            continue;
        }

        sim_.load_mesh(file);
    }
}

void Visualizer::initialize() {
    polyscope::init();

    update();
}

void Visualizer::callback() {
    draw_controls();

    if (!running_) {
        return;
    }

    for (int i = 0; i < steps_per_frame_; ++i) {
        step();
    }

    update();
}

void Visualizer::draw_controls() {
    ImGui::Text("Melting and Flowing");
    ImGui::Separator();

    ImGui::Text("File: %s", sim_.mesh_path.filename().c_str());
    ImGui::Text("Step: %zu", sim_.step_count);
    ImGui::Text("Time: %.6g", sim_.time);
}

void Visualizer::step() {
    if (dt_ <= 0.0) {
        return;
    }

    sim_.step(dt_);
}

void Visualizer::update() {
    // synchronize simulation data with Polyscope.
    // polyscope has:
    //  - VolumeGrid (need to interpolate velocities to cell centers)
    //  - PointCloud
    //
    // future implementations:
    //   update_temperature();
    //   update_pressure();
    //   update_velocities();
    //   update_surface();
}

}  // namespace visualization
