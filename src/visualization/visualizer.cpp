#include "visualizer.hpp"

#include <imgui.h>
#include <polyscope/polyscope.h>
#include <polyscope/sparse_volume_grid.h>
#include <polyscope/volume_grid.h>

#include <exception>
#include <filesystem>

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
            file.extension().c_str());
    }

    mesh_ptr = std::make_unique<Mesh>(file);

    // rescale the mesh to fit the grid, with some padding
    const auto grid_conf = sim_.storage.config;
    auto grid_min = Eigen::Vector3d{grid_conf.lowest_corner.data()};
    auto grid_max = grid_min + Eigen::Vector3d{grid_conf.grid_size.data()};
    auto padding = 0.1;

    mesh_ptr->fit_to_box(grid_min, grid_max, padding);

    // init the sim from the mesh
    sim_.load_mesh(mesh_ptr);

    // temporary, delete this, later or add as explicit option?
    // just here to display something on the screen
    polyscope::registerSurfaceMesh("Mesh", mesh_ptr->vertices, mesh_ptr->faces);

    // sync particles to polyscope for rendering
    update_particles();
    update_grid();
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

void Visualizer::update_particles() {
    const auto points = sim_.particles.pos();

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

void Visualizer::update() {
    update_grid();
    update_particles();
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

void Visualizer::update_grid() {
    const auto cfg = sim_.grid().config();
    const auto state = sim_.grid().view<grid::Field::State>();

    const auto& ext = cfg.grid_extents;

    const int nx = ext.extent(0);
    const int ny = ext.extent(1);
    const int nz = ext.extent(2);

    const glm::vec3 lower{cfg.lowest_corner[0], cfg.lowest_corner[1],
                          cfg.lowest_corner[2]};

    const glm::vec3 upper{lower.x + cfg.grid_size[0],
                          lower.y + cfg.grid_size[1],
                          lower.z + cfg.grid_size[2]};

    // Polyscope requires X to vary fastest.
    std::vector<float> values(nx * ny * nz);

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const auto index = i + nx * (j + ny * k);

                values[index] =
                    static_cast<float>(std::to_underlying(state[i, j, k]));
            }
        }
    }

    auto* ps_grid = polyscope::registerVolumeGrid(
        "Fluid Grid", glm::uvec3{nx + 1, ny + 1, nz + 1}, lower, upper);

    ps_grid->setEnabled(true);
    ps_grid->setEdgeWidth(0.5);

    auto* quantity = ps_grid->addCellScalarQuantity(
        "Cell State", values, polyscope::DataType::CATEGORICAL);

    quantity->setEnabled(true);
}

}  // namespace visualization
