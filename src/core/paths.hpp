#pragma once

#include <filesystem>

namespace fs = std::filesystem;

namespace core::paths {
inline auto assets() -> fs::path {
  static constexpr std::string_view assets_path = FLUID_ASSET_DIR;
  return fs::path{};
}

inline auto mesh(const fs::path& filename) -> fs::path {
  return assets() / "meshes" / filename;
}
}  // namespace core::paths
