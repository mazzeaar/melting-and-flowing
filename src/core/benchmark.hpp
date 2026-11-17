#pragma once

#include <chrono>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace benchmark {

#ifdef ENABLE_BENCHMARKS
inline constexpr bool enable = true;
#else
inline constexpr bool enable = false;
#endif

struct StringHash {
    using is_transparent = void;

    std::size_t operator()(std::string_view s) const noexcept {
        return std::hash<std::string_view>{}(s);
    }
};

struct BenchmarkBase {
    using Clock = std::chrono::steady_clock;
    using Duration = std::chrono::nanoseconds;
    using Map = std::unordered_map<std::string, std::size_t, StringHash,
                                   std::equal_to<>>;

    struct Measurement {
        std::string name;
        std::vector<Clock::time_point> starts;
        std::vector<Clock::time_point> ends;
    };

    std::vector<Measurement> measurements;
    Map indices;

    auto get(std::string_view key) -> Measurement& {
        if (auto it = indices.find(key); it != indices.end()) {
            return measurements[it->second];
        }

        const auto index = measurements.size();
        measurements.push_back({std::string(key), {}, {}});
        indices.emplace(measurements.back().name, index);
        return measurements.back();
    }

    void tick(std::string_view key) {
        if constexpr (enable) {
            get(key).starts.push_back(Clock::now());
        }
    }

    void tock(std::string_view key) {
        if constexpr (enable) {
            get(key).ends.push_back(Clock::now());
        }
    }
};

struct SolverBenchmark : BenchmarkBase {
    std::vector<std::size_t> cfl_substeps;

    BenchmarkBase cfl;

    void tick_cfl(std::string_view key) { cfl.tick(key); }
    void tock_cfl(std::string_view key) { cfl.tock(key); }

    void end_step(std::size_t substeps) {
        if constexpr (enable) {
            cfl_substeps.push_back(substeps);
        }
    }
};
}  // namespace benchmark
