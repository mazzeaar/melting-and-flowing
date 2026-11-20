#include "runners.hpp"

#include <algorithm>
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
    float radius_factor, int marker_stride) {
    {
        std::lock_guard lock(mutex_);
        pending_ =
            Job{std::move(snapshot), std::max(1, resolution_upfactor),
                std::max(0.01f, radius_factor), std::max(1, marker_stride)};
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
    surface_reconstructor_.resolution_upfactor = job.resolution_upfactor;
    surface_reconstructor_.radius_factor = job.radius_factor;
    surface_reconstructor_.marker_stride = job.marker_stride;

    Eigen::MatrixXf vertices;
    Eigen::MatrixXi faces;
    surface_reconstructor_.reconstruct(job.snapshot->particles.positions,
                                       job.snapshot->grid.config, vertices,
                                       faces);

    return {job.snapshot->generation, job.snapshot->step, std::move(vertices),
            std::move(faces)};
}
