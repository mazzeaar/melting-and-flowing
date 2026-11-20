#include "runners.hpp"

#include <igl/marching_cubes.h>

#include <algorithm>
#include <cmath>
#include <utility>

SimulationRunner::SimulationRunner(Simulation& sim)
    : sim_{sim}, worker_{[this](std::stop_token stop) { run(stop); }} {}

void SimulationRunner::set_running(bool running) {
    {
        std::lock_guard lock(mutex_);
        running_ = running;
        if (running) steps_ = 0;
    }
    wake_.notify_one();
}

void SimulationRunner::step() {
    {
        std::lock_guard lock(mutex_);
        running_ = false;
        ++steps_;
    }
    wake_.notify_one();
}

void SimulationRunner::reset(std::filesystem::path path,
                             SimulationSettings settings,
                             std::size_t generation) {
    {
        std::lock_guard lock(mutex_);
        running_ = false;
        steps_ = 0;
        commands_.clear();
        commands_.push_back([this, path = std::move(path), settings,
                             generation](Simulation& sim) {
            initialized_ = false;
            sim.settings = settings;
            auto mesh = std::make_unique<Mesh>(path);
            sim.load_mesh(mesh);
            source_mesh_ = std::make_shared<MeshSnapshot>(
                MeshSnapshot{generation, 0, mesh->vertices, mesh->faces});
            generation_ = generation;
            initialized_ = true;
        });
    }
    wake_.notify_one();
}

void SimulationRunner::edit(std::function<void(Simulation&)> command) {
    {
        std::lock_guard lock(mutex_);
        commands_.push_back(std::move(command));
    }
    wake_.notify_one();
}

auto SimulationRunner::latest() -> std::shared_ptr<const SimulationSnapshot> {
    std::lock_guard lock(mutex_);
    if (error_) std::rethrow_exception(std::exchange(error_, {}));
    return latest_;
}

void SimulationRunner::run(std::stop_token stop) {
    while (!stop.stop_requested()) {
        std::vector<std::function<void(Simulation&)>> commands;
        bool advance;
        {
            std::unique_lock lock(mutex_);
            if (!wake_.wait(lock, stop, [this] {
                    return running_ || steps_ != 0 || !commands_.empty();
                }))
                return;
            commands.swap(commands_);
            advance = running_ || steps_ != 0;
            if (steps_ != 0) --steps_;
        }

        // No lock during solver work, mesh loading, or snapshot copying.
        try {
            for (auto& command : commands) command(sim_);
            if (!initialized_) {
                std::lock_guard lock(mutex_);
                running_ = false;
                steps_ = 0;
                continue;
            }
            if (stop.stop_requested()) return;

            if (!commands.empty()) {
                sim_.solver.update_viscosity(Solver::Config{
                    .grid = sim_.grid,
                    .particles = sim_.particles.positions,
                    .settings = sim_.settings.solver,
                    .physics = sim_.settings.physics,
                });
            }

            const auto limit = sim_.settings.solver.max_steps;
            if (advance && (limit == 0 || sim_.step_count < limit)) {
                sim_.step(sim_.settings.solver.dt);
            }
            if (limit != 0 && sim_.step_count >= limit) {
                std::lock_guard lock(mutex_);
                running_ = false;
                steps_ = 0;
            }
            publish();
        } catch (...) {
            std::lock_guard lock(mutex_);
            error_ = std::current_exception();
            running_ = false;
            steps_ = 0;
        }
    }
}

void SimulationRunner::publish() {
    auto snapshot = std::make_shared<SimulationSnapshot>(SimulationSnapshot{
        .generation = generation_,
        .step = sim_.step_count,
        .time = sim_.time,
        .grid = {sim_.grid.config(), sim_.grid.p, sim_.grid.t, sim_.grid.eta,
                 sim_.grid.cell_state},
        .velocity = {},
        .particles = {sim_.particles.positions},
        .source_mesh = source_mesh_,
    });

    const auto& cfg = sim_.grid.config();
    const auto ex = cfg.extents();
    snapshot->velocity.vectors.resize(cfg.size(), 3);
    for (int k = 0; k < ex.z(); ++k) {
        for (int j = 0; j < ex.y(); ++j) {
            for (int i = 0; i < ex.x(); ++i) {
                const Eigen::Vector3f center =
                    cfg.min_corner() +
                    cfg.cell_size().cwiseProduct(
                        Eigen::Vector3f(i + 0.5f, j + 0.5f, k + 0.5f));
                snapshot->velocity.vectors.row(cfg.cell(i, j, k)) =
                    interpolate_velocity(sim_.grid, center).transpose();
            }
        }
    }

    std::lock_guard lock(mutex_);
    latest_ = std::move(snapshot);
}

MeshBuilderRunner::MeshBuilderRunner()
    : worker_{[this](std::stop_token stop) { run(stop); }} {}

void MeshBuilderRunner::request(
    std::shared_ptr<const SimulationSnapshot> snapshot, int resolution_upfactor,
    float radius_factor) {
    {
        std::lock_guard lock(mutex_);
        pending_ = Job{std::move(snapshot), std::max(1, resolution_upfactor),
                       std::max(0.01f, radius_factor)};
    }
    wake_.notify_one();
}

auto MeshBuilderRunner::latest() -> std::shared_ptr<const MeshSnapshot> {
    std::lock_guard lock(mutex_);
    if (error_) std::rethrow_exception(std::exchange(error_, {}));
    return latest_;
}

void MeshBuilderRunner::run(std::stop_token stop) {
    while (!stop.stop_requested()) {
        Job job;
        {
            std::unique_lock lock(mutex_);
            if (!wake_.wait(lock, stop,
                            [this] { return pending_.has_value(); })) {
                return;
            }
            job = std::move(*pending_);
            pending_.reset();
        }
        try {
            auto mesh = std::make_shared<MeshSnapshot>(build(job));
            std::lock_guard lock(mutex_);
            latest_ = std::move(mesh);
        } catch (...) {
            std::lock_guard lock(mutex_);
            error_ = std::current_exception();
        }
    }
}

auto MeshBuilderRunner::build(const Job& job) -> MeshSnapshot {
    // Surface reconstruction after Zhu & Bridson 2005, "Animating Sand as a
    // Fluid", eq. 7-10:

    // Potential speedups:
    // - Look only at cells with SURFACE (would require cell -> particles
    // mapping)
    // - Skip cells where cellstate = FULL

    // Config values (move to settings later)

    const auto& points = job.snapshot->particles.positions;
    const auto& config = job.snapshot->grid.config;
    const Eigen::Vector3f min_corner = config.min_corner();

    const Eigen::Vector3f node_spacing =
        config.cell_size() / static_cast<float>(job.resolution_upfactor);
    const Eigen::Vector3i node_counts =
        config.extents() * job.resolution_upfactor + Eigen::Vector3i::Ones();

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

    const float particle_radius = job.radius_factor * node_spacing.maxCoeff();
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

    return {job.snapshot->generation, job.snapshot->step, std::move(vertices),
            std::move(faces)};
}
