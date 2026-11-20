#include "visualizer.hpp"

#include <imgui.h>
#include <polyscope/point_cloud.h>
#include <polyscope/polyscope.h>
#include <polyscope/view.h>
#include <polyscope/volume_grid.h>

#include <Eigen/Core>
#include <algorithm>
#include <exception>
#include <filesystem>
#include <format>
#include <memory>
#include <utility>
#include <vector>

#include "core/exception.hpp"
#include "polyscope/messages.h"
#include "polyscope/surface_mesh.h"
#include "simulation/grid.hpp"
#include "simulation/simulation.hpp"

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

Visualizer::Visualizer(Simulation& sim)
    : settings_{sim.settings}, mesh_path_{sim.mesh_path}, simulation_{sim} {}

Visualizer::~Visualizer() {
    polyscope::state::filesDroppedCallback = nullptr;
    polyscope::state::userCallback = nullptr;
}

void Visualizer::init() {
    polyscope::options::giveFocusOnShow = true;
    polyscope::init();

    load_mesh(mesh_path_);

    polyscope::state::userCallback = [this] {
        handle_errors([this] {
            try {
                callback();
            } catch (...) {
                running_ = false;
                simulation_.set_running(false);
                throw;
            }
        });
    };

    polyscope::state::filesDroppedCallback = [this](auto files) {
        handle_errors([this, &files] {
            for (std::filesystem::path file : files) {
                load_mesh(file);
            }
        });
    };
}

void Visualizer::run() {
    init();
    polyscope::show();
}

void Visualizer::load_mesh(std::filesystem::path file) {
    if (file.extension() != ".obj") {
        throw core::Error(
            "Invalid file '{}', only *.obj is supported as of now",
            file.extension().string());
    }

    mesh_path_ = std::move(file);
    reset();
}

void Visualizer::callback() {
    // Poll only: these calls never wait for a timestep or mesh to finish.
    update();
    draw_controls();
    grid_stats_.draw(simulation_);
}

void Visualizer::draw_controls() {
    auto& settings = settings_;
    bool settings_changed = false;
    bool mesh_changed = false;
    std::optional<float> temperature_edit;
    auto& grid = settings.grid;
    auto& solver = settings.solver;
    auto& physics = settings.physics;
    auto& init = settings.initialization;

    ImGui::TextUnformatted("Melting and Flowing");
    ImGui::Separator();

    ImGui::Text("File: %s", mesh_path_.filename().string().c_str());
    ImGui::Text("Step: %zu", snapshot_ ? snapshot_->step : 0);
    ImGui::Text("Time: %.6g s", snapshot_ ? snapshot_->time : 0.0);
    if (!snapshot_) ImGui::TextDisabled("Loading...");

    ImGui::Spacing();

    if (ImGui::Button(running_ ? "Pause" : "Play")) {
        running_ = !running_;
        simulation_.set_running(running_);
    }

    ImGui::SameLine();
    if (ImGui::Button("Step")) {
        running_ = false;
        simulation_.step();
    }

    ImGui::SameLine();
    if (ImGui::Button("Reset")) {
        reset();
    }

    ImGui::Separator();

    // Grid
    if (ImGui::CollapsingHeader("Grid")) {
        ImGui::TextDisabled("Changes require reset");

        settings_changed |= ImGui::InputFloat3("Origin", grid.origin.data());

        settings_changed |=
            ImGui::InputInt3("Cell counts", grid.cell_counts.data());

        settings_changed |= ImGui::InputFloat("Cell size (m)", &grid.cell_size,
                                              0.001f, 0.01f, "%.6g");

        grid.cell_counts =
            grid.cell_counts.cwiseMax(Eigen::Vector3i::Constant(2));
        grid.cell_size = std::max(grid.cell_size, 1e-6f);
    }

    // Solver
    if (ImGui::CollapsingHeader("Solver", ImGuiTreeNodeFlags_DefaultOpen)) {
        settings_changed |=
            ImGui::InputFloat("dt (s)", &solver.dt, 0.0001f, 0.001f, "%.6g");

        settings_changed |= ImGui::SliderFloat("CFL", &solver.cfl, 0.01f, 1.0f);

        int max_substeps = static_cast<int>(solver.max_substeps);
        if (ImGui::InputInt("Max substeps", &max_substeps)) {
            settings_changed = true;
            solver.max_substeps =
                static_cast<std::size_t>(std::max(1, max_substeps));
        }

        int max_steps = static_cast<int>(solver.max_steps);
        if (ImGui::InputInt("Max steps (0 = unlimited)", &max_steps)) {
            settings_changed = true;
            solver.max_steps = static_cast<std::size_t>(std::max(0, max_steps));
        }

        ImGui::SeparatorText("Diffusion solvers");

        settings_changed |=
            ImGui::InputInt("Max iterations", &solver.max_solver_iterations);

        settings_changed |= ImGui::InputFloat(
            "Tolerance", &solver.solver_tolerance, 1e-7f, 1e-6f, "%.3e");
    }

    // Mesh reconstruction
    if (ImGui::CollapsingHeader("Mesh reconstructure")) {
        mesh_changed |=
            ImGui::InputInt("Resolution upfactor", &resolution_upfactor);
        resolution_upfactor = std::max(1, resolution_upfactor);

        mesh_changed |=
            ImGui::InputFloat("Particle radius factor", &radius_factor);
        radius_factor = std::max(0.01f, radius_factor);

        mesh_changed |=
            ImGui::SliderInt("Marker stride", &marker_stride, 1, 32);
        marker_stride = std::max(marker_stride, 1);
        ImGui::TextDisabled("Use every n-th marker for the surface");
    }

    // Physics
    if (ImGui::CollapsingHeader("Physics", ImGuiTreeNodeFlags_DefaultOpen)) {
        settings_changed |=
            ImGui::InputFloat("Gravity (m/s^2)", &physics.gravity, 0.1f, 1.0f);

        settings_changed |= ImGui::InputFloat("Density (kg/m^3)",
                                              &physics.density, 10.0f, 100.0f);

        settings_changed |= ImGui::SliderFloat(
            "Thermal diffusivity (m^2/s)", &physics.diffusivity, 1e-9f, 1e-3f,
            "%.3e", ImGuiSliderFlags_Logarithmic);

        float temperature = snapshot_ ? snapshot_->grid.temperature.mean()
                                      : init.initial_temperature;
        if (ImGui::InputFloat("Uniform temperature (K)", &temperature, 1.0f,
                              10.0f, "%.2f")) {
            temperature_edit = std::max(temperature, 0.0f);
        }
        ImGui::TextDisabled("Sets the current temperature throughout the grid");

        if (ImGui::TreeNode("Temperature-dependent viscosity")) {
            settings_changed |= ImGui::InputFloat(
                "Reference temperature (K)", &physics.reference_temperature);

            settings_changed |= ImGui::InputFloat("Reference viscosity (Pa s)",
                                                  &physics.reference_viscosity,
                                                  0.0f, 0.0f, "%.3e");

            settings_changed |= ImGui::SliderFloat("Thermal sensitivity (1/K)",
                                                   &physics.thermal_sensitivity,
                                                   0.0f, 0.2f, "%.4f");

            settings_changed |= ImGui::InputFloat(
                "Min viscosity", &physics.min_viscosity, 0.0f, 0.0f, "%.3e");

            settings_changed |= ImGui::InputFloat(
                "Max viscosity", &physics.max_viscosity, 0.0f, 0.0f, "%.3e");

            ImGui::TreePop();
        }
    }

    // Initialization
    if (ImGui::CollapsingHeader("Initialization")) {
        ImGui::TextDisabled("Changes require reset");

        int particle_count = static_cast<int>(init.particle_count);
        if (ImGui::InputInt("Particle count", &particle_count, 1000)) {
            settings_changed = true;
            init.particle_count =
                static_cast<std::size_t>(std::max(1, particle_count));
        }

        settings_changed |= ImGui::InputFloat("Initial temperature (K)",
                                              &init.initial_temperature);

        settings_changed |=
            ImGui::SliderFloat("Mesh padding", &init.mesh_padding, 0.0f, 0.49f);

        settings_changed |= ImGui::SliderInt(
            "Samples / cell axis", &init.samples_per_cell_axis, 1, 10);
    }

    ImGui::Separator();

    if (ImGui::Button("Restore defaults")) {
        settings = SimulationSettings{};
        reset();
    }

    // Maintain invariants.
    solver.dt = std::max(solver.dt, 1e-7f);
    solver.cfl = std::max(solver.cfl, 1e-4f);
    solver.max_substeps = std::max<std::size_t>(solver.max_substeps, 1);
    solver.max_solver_iterations = std::max(solver.max_solver_iterations, 1);
    solver.solver_tolerance = std::max(solver.solver_tolerance, 1e-12f);

    physics.density = std::max(physics.density, 1e-6f);
    physics.diffusivity = std::max(physics.diffusivity, 0.0f);
    physics.reference_viscosity = std::max(physics.reference_viscosity, 1e-12f);
    physics.min_viscosity = std::max(physics.min_viscosity, 1e-12f);
    physics.max_viscosity =
        std::max(physics.max_viscosity, physics.min_viscosity);

    init.initial_temperature = std::max(init.initial_temperature, 0.0f);
    init.samples_per_cell_axis = std::max(init.samples_per_cell_axis, 1);

    if (settings_changed || temperature_edit) {
        // Copy UI inputs into the command; never capture references to widgets.
        simulation_.edit([settings, temperature_edit](Simulation& sim) {
            sim.settings = settings;
            if (temperature_edit) sim.grid.t.setConstant(*temperature_edit);
        });
    }
    if (mesh_changed && snapshot_) {
        mesh_builder_.request(snapshot_, resolution_upfactor, radius_factor,
                              marker_stride);
    }
    if (snapshot_) {
        const auto& grid = snapshot_->grid;
        ImGui::Text("T: %.3g .. %.3g K", grid.temperature.minCoeff(),
                    grid.temperature.maxCoeff());
        ImGui::Text("eta: %.3g .. %.3g Pa s", grid.viscosity.minCoeff(),
                    grid.viscosity.maxCoeff());
    }
}

void Visualizer::update_particles() {
    const auto& points = snapshot_->particles.positions;

    if (polyscope::hasPointCloud("Particles")) {
        auto* cloud = polyscope::getPointCloud("Particles");

        if (cloud->nPoints() == static_cast<std::size_t>(points.rows())) {
            cloud->updatePointPositions(points);
            return;
        }
    }

    auto* cloud = polyscope::registerPointCloud("Particles", points);
    cloud->setPointRadius(0.002);
}

void Visualizer::update_mesh() {
    const auto& vertices = mesh_snapshot_->vertices;
    const auto& faces = mesh_snapshot_->faces;
    if (faces.rows() == 0) {
        polyscope::removeSurfaceMesh("Fluid Surface", false);
        return;
    }

    polyscope::registerSurfaceMesh("Fluid Surface", vertices, faces);
}

void Visualizer::update() {
    auto snapshot = simulation_.latest();
    if (snapshot && snapshot->generation == generation_ &&
        snapshot != snapshot_) {
        const bool rebuild = !snapshot_;
        snapshot_ = std::move(snapshot);
        grid_stats_.sample(snapshot_->statistics);
        if (rebuild && snapshot_->source_mesh) {
            const auto& source = *snapshot_->source_mesh;
            polyscope::registerSurfaceMesh("Mesh", source.vertices,
                                           source.faces)
                ->setEnabled(false);
        }
        update_grid(rebuild);
        update_particles();
        if (rebuild) polyscope::view::resetCameraToHomeView();
        mesh_builder_.request(snapshot_, resolution_upfactor, radius_factor,
                              marker_stride);
        const auto limit = settings_.solver.max_steps;
        if (limit != 0 && snapshot_->step >= limit) running_ = false;
    }

    auto mesh = mesh_builder_.latest();
    if (mesh && mesh->generation == generation_ && mesh != mesh_snapshot_) {
        mesh_snapshot_ = std::move(mesh);
        update_mesh();
    }
}

void Visualizer::update_grid(bool rebuild) {
    const auto& grid = snapshot_->grid;
    const auto cfg = grid.config;
    const auto& state = grid.cell_state;

    const auto ext = cfg.extents();

    const int nx = ext.x();
    const int ny = ext.y();
    const int nz = ext.z();

    const auto min = cfg.min_corner();
    const auto max = cfg.max_corner();

    const glm::vec3 lower{
        min.x(),
        min.y(),
        min.z(),
    };

    const glm::vec3 upper{
        max.x(),
        max.y(),
        max.z(),
    };

    // Polyscope requires X to vary fastest.
    std::vector<float> values(cfg.size());

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const auto index = i + nx * (j + ny * k);

                values[index] = static_cast<float>(state[cfg.cell(i, j, k)]);
            }
        }
    }

    // Rebuild after reset if the grid settings changed.
    const bool existing = !rebuild && polyscope::hasVolumeGrid("Fluid Grid");
    auto* ps_grid = existing
                        ? polyscope::getVolumeGrid("Fluid Grid")
                        : polyscope::registerVolumeGrid(
                              "Fluid Grid", glm::uvec3{nx + 1, ny + 1, nz + 1},
                              lower, upper);
    if (!existing) {
        ps_grid->setEnabled(false);
    }

    // Update the existing scalar quantity without replacing it.
    auto* quantity =
        ps_grid->getStructureQuantity<polyscope::VolumeGridCellScalarQuantity>(
            "Cell State");

    if (quantity) {
        quantity->updateData(values);
    } else {
        quantity = ps_grid->addCellScalarQuantity(
            "Cell State", values, polyscope::DataType::CATEGORICAL);

        // Only enable when first created.
        quantity->setEnabled(false);
    }
}

void Visualizer::reset() {
    running_ = false;
    snapshot_.reset();
    mesh_snapshot_.reset();
    // Discard old render geometry immediately; stale worker results are
    // filtered by generation even when an old timestep/reconstruction is still
    // finishing.
    polyscope::removeSurfaceMesh("Fluid Surface", false);
    polyscope::removeSurfaceMesh("Mesh", false);
    polyscope::removePointCloud("Particles", false);
    polyscope::removeVolumeGrid("Fluid Grid", false);
    simulation_.reset(mesh_path_, settings_, ++generation_);
    grid_stats_.reset(simulation_);
}

}  // namespace visualization
