#include "visualizer.hpp"

#include <igl/marching_cubes.h>
#include <imgui.h>
#include <polyscope/point_cloud.h>
#include <polyscope/polyscope.h>
#include <polyscope/sparse_volume_grid.h>
#include <polyscope/volume_grid.h>

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <format>
#include <memory>
#include <utility>
#include <vector>

#include "core/exception.hpp"
#include "core/paths.hpp"
#include "polyscope/messages.h"
#include "polyscope/surface_mesh.h"
#include "simulation/grid.hpp"
#include "simulation/mesh.hpp"
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
    load_mesh(default_mesh.string());

    polyscope::state::userCallback = [this] {
        handle_errors([this] { callback(); });
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

    mesh_ptr = std::make_unique<Mesh>(file);

    // Fit the mesh and initialise the sim using the current settings.
    sim_.load_mesh(mesh_ptr);
    grid_stats_.reset(sim_);

    // temporary, delete this, later or add as explicit option?
    // just here to display something on the screen
    auto* mesh = polyscope::registerSurfaceMesh("Mesh", mesh_ptr->vertices,
                                                mesh_ptr->faces);
    mesh->setEnabled(false);

    // sync particles to polyscope for rendering
    update_particles();
    update_grid(true);
    update_mesh();
}

void Visualizer::callback() {
    draw_controls();
    grid_stats_.draw(sim_);

    if (!running_) {
        return;
    }

    for (int i = 0; i < steps_per_frame_; ++i) {
        step();
        if (!running_) {
            break;
        }
    }

    update();
}

void Visualizer::draw_controls() {
    auto& settings = sim_.settings;
    auto& grid = settings.grid;
    auto& solver = settings.solver;
    auto& physics = settings.physics;
    auto& init = settings.initialization;

    ImGui::TextUnformatted("Melting and Flowing");
    ImGui::Separator();

    ImGui::Text("File: %s", sim_.mesh_path.filename().c_str());
    ImGui::Text("Step: %zu", sim_.step_count);
    ImGui::Text("Time: %.6g s", sim_.time);

    ImGui::Spacing();

    if (ImGui::Button(running_ ? "Pause" : "Play")) {
        running_ = !running_;
    }

    ImGui::SameLine();
    if (ImGui::Button("Step")) {
        running_ = false;
        step();
        update();
    }

    ImGui::SameLine();
    if (ImGui::Button("Reset")) {
        reset();
    }

    ImGui::SliderInt("Steps / frame", &steps_per_frame_, 1, 100);
    ImGui::Separator();

    // Grid
    if (ImGui::CollapsingHeader("Grid")) {
        ImGui::TextDisabled("Changes require reset");

        ImGui::InputFloat3("Origin", grid.origin.data());

        ImGui::InputInt3("Cell counts", grid.cell_counts.data());

        ImGui::InputFloat("Cell size (m)", &grid.cell_size, 0.001f, 0.01f,
                          "%.6g");

        grid.cell_counts =
            grid.cell_counts.cwiseMax(Eigen::Vector3i::Constant(2));
        grid.cell_size = std::max(grid.cell_size, 1e-6f);
    }

    // Solver
    if (ImGui::CollapsingHeader("Solver", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::InputFloat("dt (s)", &solver.dt, 0.0001f, 0.001f, "%.6g");

        ImGui::SliderFloat("CFL", &solver.cfl, 0.01f, 1.0f);

        int max_substeps = static_cast<int>(solver.max_substeps);
        if (ImGui::InputInt("Max substeps", &max_substeps)) {
            solver.max_substeps =
                static_cast<std::size_t>(std::max(1, max_substeps));
        }

        int max_steps = static_cast<int>(solver.max_steps);
        if (ImGui::InputInt("Max steps (0 = unlimited)", &max_steps)) {
            solver.max_steps = static_cast<std::size_t>(std::max(0, max_steps));
        }

        ImGui::SeparatorText("Diffusion solvers");

        ImGui::InputInt("Max iterations", &solver.max_solver_iterations);

        ImGui::InputFloat("Tolerance", &solver.solver_tolerance, 1e-7f, 1e-6f,
                          "%.3e");
    }

    // Mesh reconstruction
    if (ImGui::CollapsingHeader("Mesh reconstructure")) {
        ImGui::InputInt("Cell counts", &resolution_upfactor);

        ImGui::InputFloat("Particle radius factor", &radius_factor);
    }

    // Physics
    if (ImGui::CollapsingHeader("Physics", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::InputFloat("Gravity (m/s^2)", &physics.gravity, 0.1f, 1.0f);

        ImGui::InputFloat("Density (kg/m^3)", &physics.density, 10.0f, 100.0f);

        ImGui::SliderFloat("Thermal diffusivity (m^2/s)", &physics.diffusivity,
                           1e-9f, 1e-3f, "%.3e", ImGuiSliderFlags_Logarithmic);

        float temperature = sim_.grid.t.mean();
        if (ImGui::InputFloat("Uniform temperature (K)", &temperature, 1.0f,
                              10.0f, "%.2f")) {
            sim_.grid.t.setConstant(std::max(temperature, 0.0f));
        }
        ImGui::TextDisabled("Sets the current temperature throughout the grid");

        if (ImGui::TreeNode("Temperature-dependent viscosity")) {
            ImGui::InputFloat("Reference temperature (K)",
                              &physics.reference_temperature);

            ImGui::InputFloat("Reference viscosity (Pa s)",
                              &physics.reference_viscosity, 0.0f, 0.0f, "%.3e");

            ImGui::SliderFloat("Thermal sensitivity (1/K)",
                               &physics.thermal_sensitivity, 0.0f, 0.2f,
                               "%.4f");

            ImGui::InputFloat("Min viscosity", &physics.min_viscosity, 0.0f,
                              0.0f, "%.3e");

            ImGui::InputFloat("Max viscosity", &physics.max_viscosity, 0.0f,
                              0.0f, "%.3e");

            ImGui::TreePop();
        }
    }

    // Initialization
    if (ImGui::CollapsingHeader("Initialization")) {
        ImGui::TextDisabled("Changes require reset");

        int particle_count = static_cast<int>(init.particle_count);
        if (ImGui::InputInt("Particle count", &particle_count, 1000)) {
            init.particle_count =
                static_cast<std::size_t>(std::max(1, particle_count));
        }

        ImGui::InputFloat("Initial temperature (K)", &init.initial_temperature);

        ImGui::SliderFloat("Mesh padding", &init.mesh_padding, 0.0f, 0.49f);

        ImGui::SliderInt("Samples / cell axis", &init.samples_per_cell_axis, 1,
                         10);
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

    // Keep the displayed viscosity and the next solve in sync, also while
    // paused.
    sim_.solver.update_viscosity(Solver::Config{
        .grid = sim_.grid,
        .particles = sim_.particles.positions,
        .settings = settings.solver,
        .physics = settings.physics,
    });
}

void Visualizer::step() {
    const auto& settings = sim_.settings.solver;
    if (settings.max_steps != 0 && sim_.step_count >= settings.max_steps) {
        running_ = false;
        return;
    }
    if (settings.dt <= 0.0f) {
        return;
    }

    try {
        sim_.step(settings.dt);
        grid_stats_.sample(sim_);
    } catch (...) {
        running_ = false;
        throw;
    }
    if (settings.max_steps != 0 && sim_.step_count >= settings.max_steps) {
        running_ = false;
    }
}

void Visualizer::update_particles() {
    const auto& points = sim_.particles.positions;

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
    // Surface reconstruction after Zhu & Bridson 2005, "Animating Sand as a
    // Fluid", eq. 7-10:

    // Potential speedups:
    // - Look only at cells with SURFACE (would require cell -> particles
    // mapping)
    // - Skip cells where cellstate = FULL

    // Config values (move to settings later)

    const auto& points = sim_.particles.positions;
    const auto& config = sim_.grid.config();
    const Eigen::Vector3f min_corner = config.min_corner();

    const Eigen::Vector3f node_spacing =
        config.cell_size() / static_cast<float>(resolution_upfactor);
    const Eigen::Vector3i node_counts =
        config.extents() * resolution_upfactor + Eigen::Vector3i::Ones();

    const int nx = node_counts.x();
    const int ny = node_counts.y();
    const int nz = node_counts.z();
    const int n_nodes = nx * ny * nz;

    auto kernel = [](float s) {
        return std::max(0.f, std::pow(1.f - std::pow(s, 2.f), 3.f));
    };

    auto node_index = [&](int a, int b, int c) {
        return a + nx * (b + ny * c);
    };

    auto node_position = [&](int a, int b, int c) -> Eigen::Vector3f {
        return min_corner + node_spacing.cwiseProduct(
                                Eigen::Vector3i(a, b, c).cast<float>());
    };

    const float particle_radius = radius_factor * node_spacing.maxCoeff();
    const float kernel_radius = 2.f * particle_radius;

    // Splat every particle into the fine grid nodes within kernel_radius
    Eigen::VectorXf weight_sum = Eigen::VectorXf::Zero(n_nodes);
    Eigen::MatrixXf weighted_pos = Eigen::MatrixXf::Zero(n_nodes, 3);

    for (Eigen::Index i = 0; i < points.rows(); ++i) {
        const Eigen::Vector3f point = points.row(i).transpose();

        // Range of node indices within kernel_radius, clamped to the grid
        const Eigen::Vector3f local = point - min_corner;
        Eigen::Vector3i lo =
            ((local.array() - kernel_radius) / node_spacing.array())
                .ceil()
                .cast<int>()
                .matrix();
        Eigen::Vector3i hi =
            ((local.array() + kernel_radius) / node_spacing.array())
                .floor()
                .cast<int>()
                .matrix();
        lo = lo.cwiseMax(Eigen::Vector3i::Zero());
        hi = hi.cwiseMin(node_counts - Eigen::Vector3i::Ones());

        // Mesh cells around particles
        for (int c = lo.z(); c <= hi.z(); ++c) {
            for (int b = lo.y(); b <= hi.y(); ++b) {
                for (int a = lo.x(); a <= hi.x(); ++a) {
                    const Eigen::Vector3f node = node_position(a, b, c);
                    const float k =
                        kernel((node - point).norm() / kernel_radius);
                    if (k == 0.f) {
                        continue;
                    }

                    const auto n = node_index(a, b, c);
                    weight_sum(n) += k;
                    weighted_pos.row(n) += k * point.transpose();
                }
            }
        }
    }

    Eigen::VectorXf phi(n_nodes);
    Eigen::MatrixXf node_positions(n_nodes, 3);

    for (int c = 0; c < nz; ++c) {
        for (int b = 0; b < ny; ++b) {
            for (int a = 0; a < nx; ++a) {
                const auto n = node_index(a, b, c);
                const Eigen::Vector3f x = node_position(a, b, c);
                node_positions.row(n) = x.transpose();

                // Outermost layer is forced outside to close mesh
                const bool on_boundary = a == 0 || b == 0 || c == 0 ||
                                         a == nx - 1 || b == ny - 1 ||
                                         c == nz - 1;

                // No particle within kernel_radius -> outside
                if (on_boundary || weight_sum(n) == 0.f) {
                    phi(n) = kernel_radius;
                    continue;
                }

                const Eigen::Vector3f x_bar =
                    weighted_pos.row(n).transpose() / weight_sum(n);

                phi(n) = (x - x_bar).norm() - particle_radius;
            }
        }
    }

    Eigen::MatrixXf vertices;
    Eigen::MatrixXi faces;
    igl::marching_cubes(phi, node_positions, static_cast<unsigned>(nx),
                        static_cast<unsigned>(ny), static_cast<unsigned>(nz),
                        0.f, vertices, faces);

    if (faces.rows() == 0) {
        polyscope::removeSurfaceMesh("Fluid Surface", false);
        return;
    }

    polyscope::registerSurfaceMesh("Fluid Surface", vertices, faces);
}

void Visualizer::update() {
    update_grid();
    update_particles();
    update_mesh();
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

void Visualizer::update_grid(bool rebuild) {
    const auto& grid = sim_.grid;
    const auto cfg = grid.config();
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

    if (!mesh_ptr) {
        return;
    }

    sim_.load_mesh(mesh_ptr);
    grid_stats_.reset(sim_);

    auto* mesh = polyscope::registerSurfaceMesh("Mesh", mesh_ptr->vertices,
                                                mesh_ptr->faces);
    mesh->setEnabled(false);

    sim_.time = 0.0;
    sim_.step_count = 0;

    update_grid(true);
    update_particles();
    update_mesh();
}

}  // namespace visualization
