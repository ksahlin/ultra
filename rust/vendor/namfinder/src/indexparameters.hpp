#ifndef INDEXPARAMETERS_HPP
#define INDEXPARAMETERS_HPP

#include <cstdint>
#include <algorithm>
#include <iostream>
#include <limits>
#include "exceptions.hpp"

/* Settings that influence index creation */
class IndexParameters {
public:
    const int k;
    const int s;
    const int l;
    const int u;
    const int max_dist;
    const int filter_cutoff;
    const int t_syncmer;
    const unsigned w_min;
    const unsigned w_max;

    static const int DEFAULT = std::numeric_limits<int>::min();
    IndexParameters(int k, int s, int l, int u, int max_dist, int filter_cutoff)
        : k(k)
        , s(s)
        , l(l)
        , u(u)
        , t_syncmer((k - s) / 2 + 1)
        , w_min(std::max(1, l))
        , w_max(std::max(1, u))
        , max_dist(255) // Maximum allowed offset to pack in 8 bits
        , filter_cutoff(filter_cutoff)
    {
        verify();
    }

//    static IndexParameters from_read_length(int read_length, int k = DEFAULT, int s = DEFAULT, int l = DEFAULT, int u = DEFAULT, int max_seed_len = DEFAULT, int filter_cutoff = DEFAULT);
    static IndexParameters read(std::istream& os);
    std::string filename_extension() const;
    void write(std::ostream& os) const;
    bool operator==(const IndexParameters& other) const;
    bool operator!=(const IndexParameters& other) const { return !(*this == other); }

private:
    void verify() const {
        if (k <= 7 || k > 32) {
            throw BadParameter("k not in [8,32]");
        }
        if (s > k) {
            throw BadParameter("s is larger than k");
        }
        if ((k - s) % 2 != 0) {
            throw BadParameter("(k - s) must be an even number to create canonical syncmers. Please set s to e.g. k-2, k-4, k-6, ...");
        }
        if (max_dist > 255) {
            throw BadParameter("maximum seed length (-m <max_dist>) is larger than 255");
        }
    }
};

std::ostream& operator<<(std::ostream& os, const IndexParameters& parameters);

#endif
