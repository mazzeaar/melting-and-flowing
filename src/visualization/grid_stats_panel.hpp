#pragma once

#include <imgui.h>
#include <implot.h>

#include <algorithm>
#include <array>
#include <vector>

#include "simulation/grid_stats.hpp"
#include "simulation/simulation.hpp"

namespace visualization {

// UI-owned diagnostics: no solver hooks, timings, or per-substep reductions.
class GridStatisticsPanel {
    struct Sample {
        double time;
        std::size_t step;
        GridStatistics stats;
    };

    struct Metric {
        const char* label;
        double (*value)(const GridStatistics&);
    };

    static constexpr std::array metrics{
        Metric{"Kinetic", [](const auto& s) { return s.kinetic_energy; }},
        Metric{"Potential", [](const auto& s) { return s.potential_energy; }},
        Metric{"Mechanical", [](const auto& s) { return s.total_energy; }},
        Metric{"Max absolute", [](const auto& s) { return s.max_divergence; }},
        Metric{"RMS", [](const auto& s) { return s.rms_divergence; }},
        Metric{"Mean signed", [](const auto& s) { return s.mean_divergence; }},
        Metric{"Min", [](const auto& s) { return s.min_speed; }},
        Metric{"Max", [](const auto& s) { return s.max_speed; }},
        Metric{"Mean", [](const auto& s) { return s.mean_speed; }},
        Metric{"RMS", [](const auto& s) { return s.rms_speed; }},
        Metric{"Min", [](const auto& s) { return s.min_temperature; }},
        Metric{"Max", [](const auto& s) { return s.max_temperature; }},
        Metric{"Mean", [](const auto& s) { return s.mean_temperature; }},
        Metric{"Min", [](const auto& s) { return s.min_viscosity; }},
        Metric{"Max", [](const auto& s) { return s.max_viscosity; }},
        Metric{"Mean", [](const auto& s) { return s.mean_viscosity; }},
        Metric{"Min", [](const auto& s) { return s.min_pressure; }},
        Metric{"Max", [](const auto& s) { return s.max_pressure; }},
        Metric{"Mean", [](const auto& s) { return s.mean_pressure; }},
        Metric{"Value", [](const auto& s) { return s.fluid_volume; }},
        Metric{"Value", [](const auto& s) { return s.fluid_mass; }},
        Metric{"Fluid", [](const auto& s) { return double(s.fluid_cells); }},
        Metric{"Surface",
               [](const auto& s) { return double(s.surface_cells); }},
        Metric{"Invalid fluid",
               [](const auto& s) { return double(s.invalid_cells); }},
        Metric{"X", [](const auto& s) { return s.momentum.x(); }},
        Metric{"Y", [](const auto& s) { return s.momentum.y(); }},
        Metric{"Z", [](const auto& s) { return s.momentum.z(); }},
        Metric{"X", [](const auto& s) { return s.center_of_mass.x(); }},
        Metric{"Y", [](const auto& s) { return s.center_of_mass.y(); }},
        Metric{"Z", [](const auto& s) { return s.center_of_mass.z(); }},
    };

    struct Group {
        const char* label;
        int first;
        int count;
    };

    // Each group exposes only its existing views, in the metric order above.
    static constexpr std::array groups{
        Group{"Energy (J)", 0, 3},
        Group{"Divergence (1/s)", 3, 3},
        Group{"Speed (m/s)", 6, 4},
        Group{"Temperature (K)", 10, 3},
        Group{"Viscosity (Pa s)", 13, 3},
        Group{"Pressure p/rho (m^2/s^2)", 16, 3},
        Group{"Occupied volume estimate (m^3)", 19, 1},
        Group{"Mass estimate (kg)", 20, 1},
        Group{"Cell counts", 21, 3},
        Group{"Momentum (kg m/s)", 24, 3},
        Group{"Center of mass (m)", 27, 3},
    };

    struct Selection {
        int group;
        int view;
    };

    static constexpr std::size_t capacity = 1024;
    std::vector<Sample> history_;
    std::size_t next_ = 0;
    int steps_since_sample_ = 0;
    bool enabled_ = true;
    bool auto_fit_ = true;
    int sample_every_ = 5;
    std::array<Selection, 4> selected_{{{0, 2}, {1, 0}, {2, 1}, {3, 2}}};

    const Sample& chronological(std::size_t i) const {
        return history_[(i + (history_.size() == capacity ? next_ : 0)) %
                        history_.size()];
    }

    void plot(Selection& selected) const {
        bool metric_changed = false;
        const float width = std::max(1.0f, ImGui::GetContentRegionAvail().x -
                                               ImGui::GetStyle().ItemSpacing.x);
        ImGui::SetNextItemWidth(width * 0.65f);
        if (ImGui::BeginCombo("##quantity", groups[selected.group].label)) {
            for (int i = 0; i < static_cast<int>(groups.size()); ++i) {
                if (ImGui::Selectable(groups[i].label, selected.group == i) &&
                    selected.group != i) {
                    selected = {i,
                                0};  // The new quantity may have fewer views.
                    metric_changed = true;
                }
            }
            ImGui::EndCombo();
        }
        const auto& group = groups[selected.group];
        ImGui::SameLine();
        ImGui::SetNextItemWidth(width * 0.35f);
        ImGui::BeginDisabled(group.count == 1);
        if (ImGui::BeginCombo("##view",
                              metrics[group.first + selected.view].label)) {
            for (int i = 0; i < group.count; ++i) {
                if (ImGui::Selectable(metrics[group.first + i].label,
                                      selected.view == i)) {
                    selected.view = i;
                    metric_changed = true;
                }
            }
            ImGui::EndCombo();
        }
        ImGui::EndDisabled();
        const auto& current = metrics[group.first + selected.view];
        struct PlotData {
            const GridStatisticsPanel* panel;
            const Metric* metric;
        } data{this, &current};
        const auto getter = [](int i, void* raw) -> ImPlotPoint {
            const auto& d = *static_cast<PlotData*>(raw);
            const auto& sample = d.panel->chronological(i);
            return {sample.time, d.metric->value(sample.stats)};
        };

        const double value =
            current.value(chronological(history_.size() - 1).stats);
        ImGui::Text("Latest %s: %.6g", current.label, value);
        if (auto_fit_ || metric_changed) {
            ImPlot::SetNextAxesToFit();
        }
        // ImPlot is already bundled and initialized by Polyscope. The getter
        // reads the ring directly: no temporary plot arrays or float
        // conversion.
        if (ImPlot::BeginPlot("##plot", ImVec2(-1, 150),
                              ImPlotFlags_NoLegend)) {
            ImPlot::SetupAxes("Simulation time (s)", group.label);
            ImPlot::PlotLineG("Value", getter, &data,
                              static_cast<int>(history_.size()));
            ImPlot::EndPlot();
        }
    }

   public:
    GridStatisticsPanel() { history_.reserve(capacity); }

    void sample(const Simulation& sim, bool force = false) {
        if (!enabled_) {
            return;
        }
        if (!force && ++steps_since_sample_ < sample_every_) {
            return;
        }
        steps_since_sample_ = 0;
        const Sample value{
            sim.time, sim.step_count,
            GridStatistics::compute(
                sim.grid, sim.settings.physics.density,
                Eigen::Vector3f(0, sim.settings.physics.gravity, 0))};
        if (history_.size() < capacity) {
            history_.push_back(value);
        } else {
            history_[next_] = value;
        }
        next_ = (next_ + 1) % capacity;
    }

    // Reset/load starts a new history; a paused grid is never polled every
    // frame.
    void reset(const Simulation& sim) {
        history_.clear();
        next_ = 0;
        steps_since_sample_ = 0;
        sample(sim, true);
    }

    void draw(const Simulation& sim) {
        ImGui::SetNextWindowSize(ImVec2(440, 700), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Grid statistics")) {
            draw_contents(sim);
        }
        // End is required even when the window is collapsed.
        ImGui::End();
    }

   private:
    void draw_contents(const Simulation& sim) {
        if (ImGui::Checkbox("Record physical statistics", &enabled_) &&
            enabled_) {
            reset(sim);
        }
        if (ImGui::InputInt("Sample every N full steps", &sample_every_)) {
            sample_every_ = std::max(sample_every_, 1);
            steps_since_sample_ = 0;
        }
        ImGui::Checkbox("Auto-fit plots", &auto_fit_);
        if (ImGui::Button("Clear statistics")) {
            reset(sim);
        }
        if (enabled_) {
            ImGui::SameLine();
            if (ImGui::Button("Sample now")) {
                sample(sim, true);
            }
        } else {
            ImGui::TextDisabled("Recording off; history frozen; no grid scan.");
        }
        if (history_.empty()) {
            ImGui::TextDisabled("No samples recorded.");
            return;
        }

        const auto& latest = chronological(history_.size() - 1);
        const auto& s = latest.stats;
        ImGui::Text("Sample: step %zu, t = %.6g s", latest.step, latest.time);
        ImGui::Text("Fluid: %zu cells (%zu surface)", s.fluid_cells,
                    s.surface_cells);
        ImGui::Text("Volume ~ %.6g m^3; mass ~ %.6g kg", s.fluid_volume,
                    s.fluid_mass);
        ImGui::Text("Speed: %.4g .. %.4g m/s", s.min_speed, s.max_speed);
        ImGui::Text("|div| max / RMS: %.4g / %.4g 1/s", s.max_divergence,
                    s.rms_divergence);
        ImGui::Text("T: %.4g .. %.4g K (mean %.4g)", s.min_temperature,
                    s.max_temperature, s.mean_temperature);
        ImGui::Text("eta: %.4g .. %.4g Pa s", s.min_viscosity, s.max_viscosity);
        if (s.invalid_cells != 0) {
            ImGui::TextUnformatted("WARNING: non-finite fluid data!");
            ImGui::Text("%zu cells excluded from field reductions.",
                        s.invalid_cells);
        }
        ImGui::TextDisabled(
            "FULL + SURFACE cells; each counts as a whole cell.");
        ImGui::TextDisabled(
            "Energy = kinetic + gravitational (no thermal energy).");
        ImGui::TextDisabled(
            "p/rho is the solver's density-normalized pressure.");
        ImGui::TextDisabled(
            "Cell-center speeds; face differences for divergence.");
        ImGui::Text("History: %zu / %zu samples, t = %.4g .. %.4g s",
                    history_.size(), capacity, chronological(0).time,
                    latest.time);
        ImGui::TextDisabled("Disable auto-fit to pan/zoom the plots.");
        for (int i = 0; i < static_cast<int>(selected_.size()); ++i) {
            ImGui::PushID(i);
            plot(selected_[i]);
            ImGui::PopID();
        }
    }
};

}  // namespace visualization
