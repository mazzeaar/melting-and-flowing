#include <algorithm>
#include <cmath>

#include "../solver.hpp"

// Step 7: Update viscosity from temperature
// eta_new = eta(T_new)
void solver::Solver::update_viscosity(Config cfg) {
    using float_t = Solver::float_type;
    using namespace grid;

    auto temp = cfg.grid.view<Field::T>();
    auto viscosity = cfg.grid.view<Field::ETA>();

    assert(temp.size() == viscosity.size());

    // TODO: move material parameters into config or smth
    constexpr float_t T_ref = 300.0;  // K
    constexpr float_t eta_ref = 1e5;  // Pa s
    constexpr float_t beta = 0.02;    // 1/K

    constexpr float_t eta_min = 1e-3;
    constexpr float_t eta_max = 1e8;

    for (std::size_t i = 0; i < temp.size(); ++i) {
        const auto T = temp[i];
        const auto eta = eta_ref * std::exp(-beta * (T - T_ref));
        viscosity[i] = std::clamp(eta, eta_min, eta_max);
    }
}
