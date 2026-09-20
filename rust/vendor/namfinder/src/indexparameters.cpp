#include "indexparameters.hpp"
#include <vector>
#include <cmath>
#include <sstream>
#include <iostream>
#include "io.hpp"


/* Pre-defined index parameters that work well for a certain
 * "canonical" read length (and similar read lengths)  */
struct Profile {
    int canonical_read_length;
    int r_threshold;
    int k;
    int s_offset;
    int l;
    int u;
};

static auto max{std::numeric_limits<int>::max()};


void IndexParameters::write(std::ostream& os) const {
    write_int_to_ostream(os, k);
    write_int_to_ostream(os, s);
    write_int_to_ostream(os, l);
    write_int_to_ostream(os, u);
    write_int_to_ostream(os, max_dist);
    write_int_to_ostream(os, filter_cutoff);

}

IndexParameters IndexParameters::read(std::istream& is) {
    int k = read_int_from_istream(is);
    int s = read_int_from_istream(is);
    int l = read_int_from_istream(is);
    int u = read_int_from_istream(is);
    int max_dist = read_int_from_istream(is);
    int filter_cutoff = read_int_from_istream(is);
    return IndexParameters(k, s, l, u, max_dist, filter_cutoff);
}

bool IndexParameters::operator==(const IndexParameters& other) const {
    return
        this->k == other.k
        && this->s == other.s
        && this->l == other.l
        && this->u == other.u
        && this->t_syncmer == other.t_syncmer
        && this->w_min == other.w_min
        && this->w_max == other.w_max
        && this->max_dist == other.max_dist
        && this->filter_cutoff == other.filter_cutoff;
}

/*
 * Return a parameter-specific filename extension such as ".r100.sti"
 * If any of the parameters deviate from the defaults for the current
 * canonical read length, the returned extension is just ".sti".
 */
std::string IndexParameters::filename_extension() const {
    std::stringstream sstream;
//    if (*this == from_read_length(canonical_read_length)) {
//        // nothing was overridden
//        sstream << ".r" << canonical_read_length;
//    }
    sstream << ".sti";
    return sstream.str();
}

std::ostream& operator<<(std::ostream& os, const IndexParameters& parameters) {
    os << "IndexParameters("
        << ", k=" << parameters.k
        << ", s=" << parameters.s
        << ", l=" << parameters.l
        << ", u=" << parameters.u
        << ", max_dist=" << parameters.max_dist
        << ", t_syncmer=" << parameters.t_syncmer
        << ", w_min=" << parameters.w_min
        << ", w_max=" << parameters.w_max
        << ", filter_cutoff=" << parameters.filter_cutoff
        << ")";
    return os;
}
