#include <algorithm>
#include <cmath>

#include "../solver.hpp"

// Step 7: Update viscosity from temperature
// eta_new = eta(T_new)
void Solver::update_viscosity(Config cfg) {
    using float_t = Solver::float_type;

    const auto n = cfg.grid.config().size();
    auto& temp = cfg.grid.t;
    auto& viscosity = cfg.grid.eta;

    assert(temp.size() == viscosity.size());

    const auto T_ref = cfg.physics.reference_temperature;
    const auto beta = cfg.physics.thermal_sensitivity;

    const auto eta_ref = cfg.physics.reference_viscosity;
    const auto eta_min = cfg.physics.min_viscosity;
    const auto eta_max = cfg.physics.max_viscosity;

    for (std::size_t i = 0; i < n; ++i) {
        const auto T = temp[i];
        const auto eta = eta_ref * std::exp(-beta * (T - T_ref));

        viscosity[i] = std::clamp(eta, eta_min, eta_max);
    }
}
