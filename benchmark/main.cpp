#include <Eigen/Core>
#include <chrono>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <print>
#include <stdexcept>
#include <string>
#include <toml++/toml.hpp>
#include <type_traits>
#include <vector>

#include "simulation/settings.hpp"
#include "simulation/simulation.hpp"
#include "simulation/surface_reconstruction.hpp"

struct BenchConfig {
    SimulationSettings sim_settings;
    std::filesystem::path mesh_file;

    int resolution_upfactor = SurfaceReconstructor{}.resolution_upfactor;
    float radius_factor = SurfaceReconstructor{}.radius_factor;
    int marker_stride = SurfaceReconstructor{}.marker_stride;

    std::size_t num_steps = 0;
    std::size_t warmup_runs = 5;
    std::size_t runs = 5;
};

auto load_config(const std::filesystem::path& path) -> BenchConfig {
    const auto toml = toml::parse_file(path.string());
    BenchConfig config;

    // Override defaults only when a value is provided.
    const auto read = [](auto section, const char* key, auto& value) {
        value = section[key].value_or(value);
    };

    // Read Eigen vectors from TOML arrays.
    const auto read_vec = [](auto section, const char* key, auto& vec) {
        if (const auto* arr = section[key].as_array()) {
            if (arr->size() != 3) {
                throw std::runtime_error(std::string(key) +
                                         " must contain 3 elements");
            }

            for (int i = 0; i < 3; ++i) {
                using Scalar =
                    typename std::remove_cvref_t<decltype(vec)>::Scalar;

                auto value = arr->get(i)->template value<Scalar>();
                if (!value) {
                    throw std::runtime_error(std::string(key) +
                                             " contains an invalid value");
                }

                vec[i] = *value;
            }
        }
    };

    // Benchmark settings.
    const auto mesh = toml["mesh_file"].value<std::string>();
    const auto num_steps = toml["num_steps"].value<std::int64_t>();

    if (!mesh || !num_steps || *num_steps <= 0) {
        throw std::runtime_error(
            "mesh_file and positive num_steps are required");
    }

    config.mesh_file = std::filesystem::absolute(path.parent_path() / *mesh)
                           .lexically_normal();
    config.num_steps = static_cast<std::size_t>(*num_steps);

    read(toml, "warmup_runs", config.warmup_runs);
    read(toml, "runs", config.runs);

    auto& s = config.sim_settings;

    // Grid.
    const auto grid = toml["grid"];
    read_vec(grid, "origin", s.grid.origin);
    read_vec(grid, "cell_counts", s.grid.cell_counts);
    read(grid, "cell_size", s.grid.cell_size);

    // Solver.
    const auto solver = toml["solver"];
    read(solver, "dt", s.solver.dt);
    read(solver, "cfl", s.solver.cfl);
    read(solver, "max_substeps", s.solver.max_substeps);
    read(solver, "max_steps", s.solver.max_steps);
    read(solver, "max_solver_iterations", s.solver.max_solver_iterations);
    read(solver, "solver_tolerance", s.solver.solver_tolerance);

    // Physics.
    const auto physics = toml["physics"];
    read(physics, "gravity", s.physics.gravity);
    read(physics, "density", s.physics.density);
    read(physics, "diffusivity", s.physics.diffusivity);
    read(physics, "reference_temperature", s.physics.reference_temperature);
    read(physics, "reference_viscosity", s.physics.reference_viscosity);
    read(physics, "thermal_sensitivity", s.physics.thermal_sensitivity);
    read(physics, "min_viscosity", s.physics.min_viscosity);
    read(physics, "max_viscosity", s.physics.max_viscosity);

    // Initialization.
    const auto init = toml["initialization"];
    read(init, "particle_count", s.initialization.particle_count);
    read(init, "initial_temperature", s.initialization.initial_temperature);
    read(init, "mesh_padding", s.initialization.mesh_padding);
    read(init, "samples_per_cell_axis", s.initialization.samples_per_cell_axis);

    // Surface reconstruction.
    const auto mesh_section = toml["mesh"];
    read(mesh_section, "resolution_upfactor", config.resolution_upfactor);
    read(mesh_section, "radius_factor", config.radius_factor);
    read(mesh_section, "marker_stride", config.marker_stride);

    if (config.resolution_upfactor <= 0) {
        throw std::runtime_error("mesh.resolution_upfactor must be positive");
    }

    if (config.marker_stride <= 0) {
        throw std::runtime_error("mesh.marker_stride must be positive");
    }

    return config;
}

struct AverageResult {
    using Benchmark = benchmark::SolverBenchmark;
    using Duration = Benchmark::Duration;
    using Series = std::map<std::string, std::vector<double>>;

    std::vector<double> cfl_substeps;
    Series steps;
    Series mesh;
    Series cfl;
    std::size_t runs = 0;

    static double ms(Duration d) {
        return std::chrono::duration<double, std::milli>(d).count();
    }

    // Sums measurements that are recorded exactly once per full timestep.
    static void add_per_step(const benchmark::BenchmarkBase& bench,
                             std::size_t n, Series& out) {
        for (const auto& m : bench.measurements) {
            if (m.starts.size() != n || m.ends.size() != n) {
                throw std::runtime_error("Invalid step measurements: " +
                                         m.name);
            }

            auto& values = out[m.name];

            if (values.empty()) {
                values.resize(n, 0.0);
            }

            for (std::size_t i = 0; i < n; ++i) {
                values[i] += ms(m.ends[i] - m.starts[i]);
            }
        }

        // Ensure all measured runs contain the same operations.
        if (out.size() != bench.measurements.size()) {
            throw std::runtime_error("Inconsistent benchmark operations");
        }
    }

    void add(const Benchmark& bench,
             const benchmark::BenchmarkBase& mesh_bench) {
        const auto n = bench.cfl_substeps.size();

        if (runs == 0) {
            cfl_substeps.resize(n, 0.0);
        } else if (cfl_substeps.size() != n) {
            throw std::runtime_error("Inconsistent benchmark step counts");
        }

        for (std::size_t i = 0; i < n; ++i) {
            if (bench.cfl_substeps[i] == 0) {
                throw std::runtime_error("Zero CFL substeps");
            }

            cfl_substeps[i] += bench.cfl_substeps[i];
        }

        // Full-step measurements.
        add_per_step(bench, n, steps);

        // Surface reconstruction, run once after every full timestep.
        add_per_step(mesh_bench, n, mesh);

        // CFL measurements, averaged within each full timestep.
        std::size_t offset = 0;

        for (std::size_t step = 0; step < n; ++step) {
            const auto count = bench.cfl_substeps[step];

            for (const auto& m : bench.cfl.measurements) {
                if (offset + count > m.starts.size() ||
                    offset + count > m.ends.size()) {
                    throw std::runtime_error("Invalid CFL measurements");
                }

                auto& values = cfl[m.name];

                if (values.empty()) {
                    values.resize(n, 0.0);
                }

                double total = 0.0;

                for (std::size_t i = 0; i < count; ++i) {
                    total += ms(m.ends[offset + i] - m.starts[offset + i]);
                }

                values[step] += total / count;
            }

            offset += count;
        }

        for (const auto& m : bench.cfl.measurements) {
            if (m.starts.size() != offset || m.ends.size() != offset) {
                throw std::runtime_error("Unexpected CFL measurements");
            }
        }

        // Ensure all measured runs contain the same operations.
        if (cfl.size() != bench.cfl.measurements.size()) {
            throw std::runtime_error("Inconsistent benchmark operations");
        }

        ++runs;
    }

    void print() const {
        if (runs == 0) {
            return;
        }

        std::print("{:>6} {:>10}", "Step", "CFL avg");

        for (const auto& [name, _] : steps) {
            std::print(" {:>20}", name);
        }

        for (const auto& [name, _] : mesh) {
            std::print(" {:>20}", name);
        }

        for (const auto& [name, _] : cfl) {
            std::print(" {:>20}", name + " (CFL)");
        }

        std::println("");

        for (std::size_t i = 0; i < cfl_substeps.size(); ++i) {
            std::print("{:>6} {:>10.2f}", i, cfl_substeps[i] / runs);

            for (const auto& [name, values] : steps) {
                std::print(" {:>20.3f}", values[i] / runs);
            }

            for (const auto& [name, values] : mesh) {
                std::print(" {:>20.3f}", values[i] / runs);
            }

            for (const auto& [name, values] : cfl) {
                std::print(" {:>20.3f}", values[i] / runs);
            }

            std::println("");
        }
    }

    void save_csv(const std::filesystem::path& path) const {
        if (runs == 0) {
            throw std::runtime_error("No benchmark results to save");
        }

        if (path.has_parent_path()) {
            std::filesystem::create_directories(path.parent_path());
        }

        std::ofstream file(path);

        if (!file) {
            throw std::runtime_error("Failed to open CSV: " + path.string());
        }

        file << std::setprecision(10);

        // CSV header.
        file << "step,cfl_substeps";

        for (const auto& [name, _] : steps) {
            file << ',' << name << "_ms";
        }

        for (const auto& [name, _] : mesh) {
            file << ',' << name << "_ms";
        }

        for (const auto& [name, _] : cfl) {
            file << ',' << name << "_cfl_ms";
        }

        file << '\n';

        // One row per full timestep.
        for (std::size_t i = 0; i < cfl_substeps.size(); ++i) {
            file << i << ',' << cfl_substeps[i] / runs;

            for (const auto& [name, values] : steps) {
                file << ',' << values[i] / runs;
            }

            for (const auto& [name, values] : mesh) {
                file << ',' << values[i] / runs;
            }

            for (const auto& [name, values] : cfl) {
                file << ',' << values[i] / runs;
            }

            file << '\n';
        }

        if (!file) {
            throw std::runtime_error("Failed to write CSV: " + path.string());
        }
    }
};

int main(int argc, char** argv) {
    if (argc != 2) {
        std::println(stderr, "Usage: fluid_benchmark <config.toml>");
        return 1;
    }

    try {
        const auto config = load_config(argv[1]);

        if (config.runs == 0) {
            throw std::runtime_error("runs must be positive");
        }

        AverageResult averages;

        for (std::size_t run = 0; run < config.warmup_runs + config.runs;
             ++run) {
            const bool warmup = run < config.warmup_runs;

            std::println("{} run {}/{}", warmup ? "Warmup" : "Benchmark",
                         warmup ? run + 1 : run - config.warmup_runs + 1,
                         warmup ? config.warmup_runs : config.runs);

            // Fresh simulation for each run.
            auto sim = Simulation{config.sim_settings};
            auto mesh = std::make_unique<Mesh>(config.mesh_file);

            sim.load_mesh(mesh);

            auto reconstructor = SurfaceReconstructor{};
            reconstructor.resolution_upfactor = config.resolution_upfactor;
            reconstructor.radius_factor = config.radius_factor;
            reconstructor.marker_stride = config.marker_stride;

            Eigen::MatrixXf vertices;
            Eigen::MatrixXi faces;

            for (std::size_t step = 0; step < config.num_steps; ++step) {
                sim.step(config.sim_settings.solver.dt);

                // Timed outside Solver::step so total_ms stays comparable.
                reconstructor.reconstruct(sim.particles.positions,
                                          sim.grid.config(), vertices, faces);
            }

            // Discard warm-up measurements.
            if (!warmup) {
                averages.add(sim.solver.get_benchmark(),
                             reconstructor.benchmark);
            }
        }

        std::println("\nAverage timings over {} runs (ms):", averages.runs);
        averages.print();

        const auto name = std::filesystem::path(argv[1]).stem().string();

        const auto timestamp =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now().time_since_epoch())
                .count();

        const auto output = std::filesystem::path("benchmark/results") / name /
                            (std::to_string(timestamp) + ".csv");

        averages.save_csv(output);
        std::println("\nSaved benchmark results to {}", output.string());

        return 0;
    } catch (const std::exception& e) {
        std::println(stderr, "Benchmark failed: {}", e.what());
        return 1;
    }
}
