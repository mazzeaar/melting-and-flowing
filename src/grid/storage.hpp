#pragma once

#include <format>
#include <print>

#include "core.hpp"
#include "core/exception.hpp"
#include "core/log.hpp"
#include "interface.hpp"

// forward declaration
struct Mesh;

namespace grid {

struct Storage {
    friend class GridView;
    Config config;

    // just using flat storage for now
    // all external memory interactions must only happen through view()
    std::vector<Config::float_type> flat_storage;
    std::vector<CellState> state;

    Storage(Config config) : config{config} {
        constexpr auto n_fields = std::to_underlying(Field::COUNT);
        auto buffer_size = config.grid_mapping.required_span_size();
        auto total_size = n_fields * buffer_size;

        state = std::vector<CellState>(total_size, CellState::AIR);
        flat_storage = std::vector<Config::float_type>(total_size, 0.0);
        auto kb = total_size / (sizeof(Config::float_type) * 1000);

        auto format_extents = [](const auto& ext, int offset = 0) {
            return std::format("[{}, {}, {}]", ext.extent(0) + offset,
                               ext.extent(1) + offset, ext.extent(2) + offset);
        };

        Log("Initialised FieldStorage:\n"
            "  field count: {}\n"
            "  allocated: {}kb\n"
            "  number of cells: {}\n"
            "  extended grid size: {}\n"
            "  physical grid size: {}\n",
            n_fields, kb, config.grid_mapping.required_span_size(),
            format_extents(config.grid_extents),
            format_extents(config.grid_extents, -1));
        if (kb > 40'000) {
            Log("\n\ncan your L3 really handle {}mb?\n\n", kb / 1000);
        }
    }

    template <Field F> constexpr auto view(this auto&& self) {
        static constexpr auto field_offset = std::to_underlying(F);
        auto field_size = self.config.grid_mapping.required_span_size();

        using enum Field;
        if constexpr (std::ranges::contains(std::array{U, V, W, P}, F)) {
            return FieldInterface<F>{
                {self.flat_storage.data() + (field_size * field_offset),
                 field_size},
                self.config,
            };
        } else if constexpr (F == Velocities) {
            return FieldInterface<F>{
                {self.flat_storage.data(), 3 * field_size},
                self.config,
            };
        } else if constexpr (F == State) {
            return FieldInterface<F>{{self.state.data(), field_size},
                                     self.config};
        } else {
            static_assert(
                false,
                "state? STATE?, hit the road dude, this aint a valid enum");
        }
    }

    // maybe this also needs a way to copy in buffers?
    // problem is: swap does not work due to the combined memory requirement of
    // the 3 velocity fields
};

class GridView {
    Storage* storage;

   public:
    constexpr GridView(Storage& s) : storage{&s} {}
    constexpr auto config() const -> const Config& { return storage->config; }

    // maybe we should add
    //   - mapping()
    //   - extents()
    //   - other stuff?
    // here as well idk

    void initialise_cell_states(const Mesh& mesh);
    void update_cell_states(
        const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& pos);

    void advect_particles(Eigen::MatrixXd& positions) const;

    template <Field F> constexpr auto view(this auto&& self) {
        return self.storage->template view<F>();
    }

    template <Field F, class T>
    constexpr auto resize_like(std::vector<T>& buff) {
        // geometric fields are distributed evenly throughout the grid, hence
        // all the same buffer size
        if constexpr (detail::GeometricField<F>) {
            buff.resize(storage->config.grid_mapping.required_span_size());
        } else if constexpr (F == Field::Velocities) {
            buff.resize(3 * storage->config.grid_mapping.required_span_size());
        } else {
            // assert(false) will always crash...
            static_assert(F != F,  // NOLINT
                          "handle your field here, or in "
                          "detail::GeometricField");
        }
    }
};
}  // namespace grid
