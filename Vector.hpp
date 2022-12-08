/// ---------------------------------------------------------------------------
///   HPC Course Fall 2022
///   Indiana University Bloomington
///
///   Luke D'Alessandro
///
/// ---------------------------------------------------------------------------
#pragma once

#include <vector>

/// The Vector class template just wraps a std::vector.
///
/// This allows end-users to uniformly use `operator()` syntax to access both
/// Matrix and Vector data.
template <class T>
struct Vector
{
    int shape[1];
    std::vector<T> _data;

    Vector(int n) : _data(n) {
        shape[0] = n;
    }

    /// These two permit indexing with the () operator.
    /// @{
    auto operator()(int i) const -> T const& {
        return _data[i];
    }

    auto operator()(int i) -> T& {
        return _data[i];
    }
    /// @}

    /// Everything below just forwards to std::vector directly.
    /// @{
    auto size() const {
        return _data.size();
    }

    auto begin() const {
        return _data.begin();
    }

    auto end() const {
        return _data.end();
    }

    auto data() const {
        return _data.data();
    }

    auto data() {
        return _data.data();
    }
    /// @}
};