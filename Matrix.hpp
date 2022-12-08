/// ---------------------------------------------------------------------------
///   HPC Course Fall 2022
///   Indiana University Bloomington
///
///   Luke D'Alessandro
///
/// ---------------------------------------------------------------------------
#pragma once

#include <vector>

/// Our column major layout function.
inline auto column_major(int i, int j, int n_rows) -> int {
    return i + j * n_rows;
}

/// Our row major layout function.
inline auto row_major(int i, int j, int n_columns) -> int {
    return i * n_columns + j;
}

template <class T>
struct Matrix
{
    int n_rows;
    int n_columns;
    std::vector<T> _data;

    /// Basic constructor, set n_rows and n_columns and allocate data.
    Matrix(int n_rows, int n_columns)
            : n_rows(n_rows)
            , n_columns(n_columns)
            , _data(n_rows * n_columns)
    {
    }

    /// 2D indexing operation for const matrices.
    auto operator()(int i, int j) const -> T const& {
        int k = column_major(i, j, n_rows);
        return _data[k];
    }

    /// 2D indexing operation for matrices.
    auto operator()(int i, int j) -> T& {
        int k = column_major(i, j, n_rows);
        return _data[k];
    }
};