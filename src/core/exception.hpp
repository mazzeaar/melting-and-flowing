#pragma once

#include <format>
#include <stdexcept>
#include <utility>

namespace core {

class Error : public std::runtime_error {
   public:
    using std::runtime_error::runtime_error;

    template <class... Args>
    explicit Error(std::format_string<Args...> fmt, Args&&... args)
        : std::runtime_error(std::format(fmt, std::forward<Args>(args)...)) {}
};
}  // namespace core
