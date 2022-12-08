/// ---------------------------------------------------------------------------
///   HPC Course Fall 2022
///   Indiana University Bloomington
///
///   Luke D'Alessandro
///
/// ---------------------------------------------------------------------------
#pragma once

#include <fmt/format.h>
#include <cstdio>                               // std::sscanf
#include <cstring>                              // std::strchr
#include <string>

/// A convenience class for reading information out of an mmio file.
///
/// Usage:
///
///     // path is a c-string or c++ string from the command line
///     std::string path = ...;
///     char const* path = ...;
///
///     MMIOReader mmio(path); // construct a reader for the file
///
///     Matrix<double> A(mmio.n_rows, mmio.n_columns); // allocate a matrix
///
///     for (auto&& [i, j, w] : mmio) { // for each edge (uses structured bindings)
///         if (this_rank_owns(i)) {
///             A[global_to_local_row(i)][j] = w;
///         }
///     }
///
struct MMIOReader
{
    int n_rows      = 0;
    int n_columns   = 0;
    int n_non_zeros = 0;

    bool is_symmetric = false;
    bool has_weights  = false;

    char const* _data  = nullptr;
    char const* _edges = nullptr;
    char const* _end   = nullptr;

    /// Create a reader for the passed file.
    MMIOReader(std::string path);
    ~MMIOReader();

    /// Readers can't be copied
    MMIOReader(MMIOReader const&) = delete;
    MMIOReader(MMIOReader&&) noexcept = default;

    /// This type encodes the i,j,w we read from the graph as edges.
    struct edge_record {
        int i = -1;
        int j = -1;
        double w = 1.0;
    };

    /// This type allows you to iterate through the edges.
    struct const_iterator
    {
        char const* _i;
        bool _has_weights;
        bool _is_symmetric;
        int _output_symmetric = 0;

        /// Read the current line as an edge record.
        auto operator*() const -> edge_record
        {
            // These are used in the strto* functions to parse values (you can
            // read about strtol,strtod on cppreference).
            char const *i = _i;
            char *e = nullptr;

            // The record we'll return.
            edge_record out;

            // Read the i, j, and possibly w values.
            out.i = std::strtol(i, &e, 10) - 1; i = e; // update i to the next value
            out.j = std::strtol(i, &e, 10) - 1; i = e; // update i to the next value
            if (_has_weights) {
                out.w = std::strtod(i, &e);
            }

            // Maybe swap this edge.
            if (_is_symmetric and _output_symmetric) {
                std::swap(out.i, out.j);
            }

            // Return the record.
            return out;
        }

        /// Advance to the next line.
        auto operator++() -> const_iterator&
        {
            if (_is_symmetric) {
                // this bounces back and forth between 0 and 1, and if its one
                // it means we don't want to advance to the next edge, we want
                // to first output the symmetric edge
                if ((_output_symmetric = 1 - _output_symmetric) == 1) {
                    return *this;
                }
            }
            _i = std::strchr(_i, '\n') + 1;
            return *this;
        }

        bool operator!=(const_iterator const& b) const {
            return _i != b._i;
        }
    };

    /// An iterator to the first edge.
    auto begin() const -> const_iterator
    {
        return const_iterator {
                ._i = _edges,
                ._has_weights = has_weights,
                ._is_symmetric = is_symmetric
        };
    }

    /// An iterator one-past-the-end of the edges.
    auto end() const -> const_iterator
    {
        return const_iterator {
                ._i = _end,
                ._has_weights = has_weights,
                ._is_symmetric = is_symmetric
        };
    }
};