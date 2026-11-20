#pragma once

#include <condition_variable>
#include <exception>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <thread>
#include <vector>

#include "simulation/simulation.hpp"
#include "snapshot.hpp"

// While this runner exists, only its worker may access Simulation.
// Callers read immutable snapshots and send edits to the worker.
class SimulationRunner {
   public:
    explicit SimulationRunner(Simulation& sim);
    SimulationRunner(const SimulationRunner&) = delete;
    SimulationRunner& operator=(const SimulationRunner&) = delete;

    void set_running(bool running);
    void step();
    void reset(std::filesystem::path path, SimulationSettings settings,
               std::size_t generation);

    // Runs between timesteps, including while paused. Capture inputs by value.
    // Future heat/laser inputs can use this same path to modify the live grid.
    void edit(std::function<void(Simulation&)> command);
    auto latest() -> std::shared_ptr<const SimulationSnapshot>;

   private:
    void run(std::stop_token stop);
    void publish();

    Simulation& sim_;
    bool initialized_ = false;                         // worker only
    std::size_t generation_ = 0;                       // worker only
    std::shared_ptr<const MeshSnapshot> source_mesh_;  // worker only

    std::mutex mutex_;
    std::condition_variable_any wake_;
    bool running_ = false;
    std::size_t steps_ = 0;
    std::vector<std::function<void(Simulation&)>> commands_;
    std::shared_ptr<const SimulationSnapshot> latest_;
    std::exception_ptr error_;
    // Last member: jthread stops/joins before the data above is destroyed.
    std::jthread worker_;
};

class MeshBuilderRunner {
   public:
    MeshBuilderRunner();
    MeshBuilderRunner(const MeshBuilderRunner&) = delete;
    MeshBuilderRunner& operator=(const MeshBuilderRunner&) = delete;

    void request(std::shared_ptr<const SimulationSnapshot> snapshot,
                 int resolution_upfactor, float radius_factor);
    auto latest() -> std::shared_ptr<const MeshSnapshot>;

   private:
    struct Job {
        std::shared_ptr<const SimulationSnapshot> snapshot;
        int resolution_upfactor;
        float radius_factor;
    };
    void run(std::stop_token stop);
    static auto build(const Job& job) -> MeshSnapshot;

    std::mutex mutex_;
    std::condition_variable_any wake_;
    // A newer request replaces the pending one; no unbounded mesh queue.
    std::optional<Job> pending_;
    std::shared_ptr<const MeshSnapshot> latest_;
    std::exception_ptr error_;
    std::jthread worker_;
};
