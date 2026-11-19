#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

#include "grid.hpp"

struct GridStatistics {
    double kinetic_energy = 0.0;
    double potential_energy = 0.0;
    double total_energy = 0.0;

    double min_speed = 0.0;
    double max_speed = 0.0;
    double mean_speed = 0.0;
    double rms_speed = 0.0;
    double max_divergence = 0.0;  // max absolute MAC divergence, 1/s
    double rms_divergence = 0.0;
    double mean_divergence = 0.0;  // signed; can cancel, unlike RMS/max

    double min_temperature = 0.0;
    double max_temperature = 0.0;
    double mean_temperature = 0.0;
    double min_viscosity = 0.0;
    double max_viscosity = 0.0;
    double mean_viscosity = 0.0;

    // The current pressure solver uses rho = 1: grid.p stores p/rho, not Pa.
    double min_pressure = 0.0;
    double max_pressure = 0.0;
    double mean_pressure = 0.0;

    // Occupancy estimates: each FULL/SURFACE cell counts as one whole cell.
    double fluid_volume = 0.0;
    double fluid_mass = 0.0;
    Eigen::Vector3d momentum = Eigen::Vector3d::Zero();
    Eigen::Vector3d center_of_mass = Eigen::Vector3d::Zero();

    std::size_t fluid_cells = 0;
    std::size_t surface_cells = 0;
    std::size_t invalid_cells = 0;

    // Read-only, allocation-free pass. Call AFTER a completed simulation step,
    // outside the solver/benchmarks, and only at the desired sampling cadence.
    // Reductions exclude fluid cells with non-finite inputs; occupancy does not.
    static auto compute(const Grid& grid, float density,
                        const Eigen::Vector3f& gravity) -> GridStatistics {
        GridStatistics s;

        const auto& cfg = grid.config();
        const auto ex = cfg.extents();
        const auto dx = cfg.cell_size();
        const double volume = dx.cast<double>().prod();
        const Eigen::Vector3d origin = cfg.min_corner().cast<double>();
        const Eigen::Vector3d h = dx.cast<double>();
        const Eigen::Vector3d g = gravity.cast<double>();
        const double mass = density * volume;

        double speed_min = std::numeric_limits<double>::infinity();
        double temp_min = std::numeric_limits<double>::infinity();
        double temp_max = -std::numeric_limits<double>::infinity();
        double eta_min = std::numeric_limits<double>::infinity();
        double eta_max = -std::numeric_limits<double>::infinity();
        double pressure_min = std::numeric_limits<double>::infinity();
        double pressure_max = -std::numeric_limits<double>::infinity();
        double speed_sq_sum = 0.0;
        double div_sq_sum = 0.0;

        for (int k = 0; k < ex.z(); ++k) {
            for (int j = 0; j < ex.y(); ++j) {
                for (int i = 0; i < ex.x(); ++i) {
                    const auto c = cfg.cell(i, j, k);
                    const int state = grid.cell_state[c];

                    if (state != static_cast<int>(CellState::FULL) &&
                        state != static_cast<int>(CellState::SURFACE))
                        continue;

                    ++s.fluid_cells;
                    s.surface_cells += state == static_cast<int>(CellState::SURFACE);

                    // Velocity arrays have a one-element ghost layer on ALL
                    // axes. Physical cell (i,j,k) lies between these faces.
                    const double u0 = grid.u[cfg.u(i + 1, j + 1, k + 1)];
                    const double u1 = grid.u[cfg.u(i + 2, j + 1, k + 1)];
                    const double v0 = grid.v[cfg.v(i + 1, j + 1, k + 1)];
                    const double v1 = grid.v[cfg.v(i + 1, j + 2, k + 1)];
                    const double w0 = grid.w[cfg.w(i + 1, j + 1, k + 1)];
                    const double w1 = grid.w[cfg.w(i + 1, j + 1, k + 2)];
                    const double temp = grid.t[c];
                    const double eta = grid.eta[c];
                    const double pressure = grid.p[c];

                    if (!std::isfinite(u0) || !std::isfinite(u1) ||
                        !std::isfinite(v0) || !std::isfinite(v1) ||
                        !std::isfinite(w0) || !std::isfinite(w1) ||
                        !std::isfinite(temp) || !std::isfinite(eta) ||
                        !std::isfinite(pressure)) {
                        ++s.invalid_cells;
                        continue;
                    }

                    const Eigen::Vector3d velocity{
                        0.5 * (u0 + u1), 0.5 * (v0 + v1), 0.5 * (w0 + w1)};

                    const double speed_sq = velocity.squaredNorm();
                    const double speed = std::sqrt(speed_sq);
                    speed_min = std::min(speed_min, speed);
                    s.max_speed = std::max(s.max_speed, speed);
                    s.mean_speed += speed;
                    speed_sq_sum += speed_sq;

                    const double divergence = (u1 - u0) / dx.x() +
                                              (v1 - v0) / dx.y() +
                                              (w1 - w0) / dx.z();

                    s.max_divergence = std::max(s.max_divergence, std::abs(divergence));
                    s.mean_divergence += divergence;
                    div_sq_sum += divergence * divergence;

                    const Eigen::Vector3d relative_position =
                        Eigen::Vector3d(i + 0.5, j + 0.5, k + 0.5).cwiseProduct(h);

                    s.kinetic_energy += 0.5 * mass * speed_sq;
                    // Zero potential at the domain origin, using actual gravity.
                    s.potential_energy -= mass * g.dot(relative_position);
                    s.momentum += mass * velocity;
                    s.center_of_mass += origin + relative_position;

                    temp_min = std::min(temp_min, temp);
                    temp_max = std::max(temp_max, temp);
                    s.mean_temperature += temp;
                    eta_min = std::min(eta_min, eta);
                    eta_max = std::max(eta_max, eta);
                    s.mean_viscosity += eta;
                    pressure_min = std::min(pressure_min, pressure);
                    pressure_max = std::max(pressure_max, pressure);
                    s.mean_pressure += pressure;
                }
            }
        }

        s.fluid_volume = s.fluid_cells * volume;
        s.fluid_mass = density * s.fluid_volume;
        const auto valid_cells = s.fluid_cells - s.invalid_cells;
        if (valid_cells > 0) {
            const double inv_count = 1.0 / valid_cells;
            s.min_speed = speed_min;
            s.mean_speed *= inv_count;
            s.rms_speed = std::sqrt(speed_sq_sum * inv_count);
            s.min_temperature = temp_min;
            s.max_temperature = temp_max;
            s.mean_temperature *= inv_count;
            s.min_viscosity = eta_min;
            s.max_viscosity = eta_max;
            s.mean_viscosity *= inv_count;
            s.min_pressure = pressure_min;
            s.max_pressure = pressure_max;
            s.mean_pressure *= inv_count;
            s.mean_divergence *= inv_count;
            s.rms_divergence = std::sqrt(div_sq_sum * inv_count);
            s.center_of_mass *= inv_count;
        }

        s.total_energy = s.kinetic_energy + s.potential_energy;
        return s;
    }
};
