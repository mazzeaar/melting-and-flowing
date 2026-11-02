#include "visualizer.hpp"

#include <imgui.h>
#include <polyscope/polyscope.h>

#include <exception>

#include "core/exception.hpp"
#include "core/paths.hpp"
#include "mesh/mesh.hpp"
#include "polyscope/messages.h"
#include "polyscope/surface_mesh.h"
#include "simulation.hpp"

namespace visualization {

template <typename F> void handle_errors(F&& fn) {
    try {
        std::forward<F>(fn)();
    } catch (const core::Error& e) {
        polyscope::warning("Error", e.what());
    } catch (const std::exception& e) {
        polyscope::error(std::format("Unexpected error: {}", e.what()));
    } catch (...) {
        polyscope::error("Unknown exception");
    }
}

Visualizer::Visualizer(Simulation& sim) : sim_{sim} {}

Visualizer::~Visualizer() {
    // probably leaks memory
    polyscope::state::filesDroppedCallback = nullptr;
    polyscope::state::userCallback = nullptr;
}

void Visualizer::init() {
    polyscope::options::giveFocusOnShow = true;
    polyscope::init();

    update();

    const auto default_mesh = core::paths::mesh("bunny.obj");
    load_mesh(std::array{default_mesh.string()});

    polyscope::state::userCallback = [this] {
        handle_errors([this] { callback(); });
    };

    polyscope::state::filesDroppedCallback = [this](auto files) {
        handle_errors([this, &files] { load_mesh(files); });
    };
}

void Visualizer::run() {
    init();
    polyscope::show();
}

void Visualizer::load_mesh(std::span<const std::string> files) {
    for (std::filesystem::path file : files) {
        if (file.extension() != ".obj") {
            throw core::Error(
                "Invalid file '{}', only *.obj is supported as of now",
                file.extension().c_str());
        }

        mesh_ptr = std::make_unique<Mesh>(file);
        sim_.load_mesh(mesh_ptr);
        // temporary, delete this
        // just here to display something on the screen
        polyscope::registerSurfaceMesh("Mesh", mesh_ptr->vertices,
                                       mesh_ptr->faces);
    }
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
