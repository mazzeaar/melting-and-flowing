#pragma once

#include <format>
#include <print>
#include <source_location>
#include <utility>

namespace logging {
template <class... Args>
inline constexpr void log(std::source_location loc,
                          std::format_string<Args...> fmt, Args&&... args) {
    auto file = std::string_view{loc.file_name()};
    if (auto pos = file.rfind("/src/"); pos != std::string_view::npos) {
        file.remove_prefix(pos + 5);
    }

    std::print("[{}:{}] ", file, loc.line());
    std::println(fmt, std::forward<Args>(args)...);
}

}  // namespace logging

#ifdef NDEBUG
#    define Log(...) ((void)0)
#else
#    define Log(...) logging::log(std::source_location::current(), __VA_ARGS__)
#endif
