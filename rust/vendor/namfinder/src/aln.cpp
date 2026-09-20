#include "aln.hpp"

#include <iostream>
#include <math.h>
#include <sstream>
#include "revcomp.hpp"
#include "timer.hpp"
#include "nam.hpp"
#include "output.hpp"
//#include "aligner.hpp"

#include "logger.hpp"

static Logger& logger = Logger::get();

using namespace klibpp;



static inline bool score(const Nam &a, const Nam &b) {
    return a.score > b.score;
}

/*
 * Determine whether the NAM represents a match to the forward or
 * reverse-complemented sequence by checking in which orientation the
 * first and last strobe in the NAM match
 *
 * - If first and last strobe match in forward orientation, return true.
 * - If first and last strobe match in reverse orientation, update the NAM
 *   in place and return true.
 * - If first and last strobe do not match consistently, return false.
 */
bool reverse_nam_if_needed(Nam& n, const Read& read, const References& references, int k) {
    auto read_len = read.size();
    std::string ref_start_kmer = references.sequences[n.ref_id].substr(n.ref_s, k);
    std::string ref_end_kmer = references.sequences[n.ref_id].substr(n.ref_e-k, k);

    std::string seq, seq_rc;
    if (n.is_rc) {
        seq = read.rc;
        seq_rc = read.seq;
    } else {
        seq = read.seq;
        seq_rc = read.rc;
    }
    std::string read_start_kmer = seq.substr(n.query_s, k);
    std::string read_end_kmer = seq.substr(n.query_e-k, k);
    if (ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer) {
        return true;
    }

    // False forward or false reverse (possible due to symmetrical hash values)
    //    we need two extra checks for this - hopefully this will remove all the false hits we see (true hash collisions should be very few)
    int q_start_tmp = read_len - n.query_e;
    int q_end_tmp = read_len - n.query_s;
    // false reverse hit, change coordinates in nam to forward
    read_start_kmer = seq_rc.substr(q_start_tmp, k);
    read_end_kmer = seq_rc.substr(q_end_tmp - k, k);
    if (ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer) {
        n.is_rc = !n.is_rc;
        n.query_s = q_start_tmp;
        n.query_e = q_end_tmp;
        return true;
    }
    return false;
}








static inline bool compareByQueryCoord(const Nam &a, const Nam &b)
{
    // first sort on ref ID, then on query, then on reference
    return (a.ref_id < b.ref_id) ||
           ( (a.ref_id == b.ref_id) && (a.query_s < b.query_s) ) ||
           ((a.ref_id == b.ref_id) && (a.query_s == b.query_s ) && (a.ref_s < b.ref_s)) ;
}


void align_SE_read(
    const KSeq &record,
    std::string &outstring,
    AlignmentStatistics &statistics,
    const mapping_params &map_param,
    const IndexParameters& index_parameters,
    const References& references,
    const StrobemerIndex& index
) {
    Timer strobe_timer;
    auto query_randstrobes = randstrobes_query(index_parameters.k, index_parameters.w_min, index_parameters.w_max, record.seq, index_parameters.s, index_parameters.t_syncmer, index_parameters.max_dist);
    statistics.tot_construct_strobemers += strobe_timer.duration();

    // Find NAMs
    Timer nam_timer;

//    logger.debug() << "index_parameters.filter_cutoff: " << std::to_string(index_parameters.filter_cutoff)  << "index.filter_cutoff: " << std::to_string(index.filter_cutoff) << std::endl;

    auto [nonrepetitive_fraction, nams] = find_nams(query_randstrobes, index);
    statistics.tot_find_nams += nam_timer.duration();

//    if (map_param.R > 1) {
//        Timer rescue_timer;
//        if (nams.empty() || nonrepetitive_fraction < 0.7) {
//            statistics.tried_rescue += 1;
//            nams = find_nams_rescue(query_randstrobes, index, map_param.rescue_cutoff);
//        }
//        statistics.tot_time_rescue += rescue_timer.duration();
//    }

    Timer nam_sort_timer;
    std::sort(nams.begin(), nams.end(), score);
    statistics.tot_sort_nams += nam_sort_timer.duration();

    // Take first L NAMs for output
    unsigned int cut_nam_vec_at = (map_param.L < nams.size()) ? map_param.L : nams.size();
    std::vector<Nam> nams_cut(nams.begin(), nams.begin() + cut_nam_vec_at);

    //Sort hits based on start choordinate on query sequence
    if (!map_param.sort_on_scores) {
//        logger.debug() << "Sorting output on scores. sort_on_scores: " << std::endl;
        std::sort(nams_cut.begin(), nams_cut.end(), compareByQueryCoord);
    }

    output_nams(outstring, nams_cut, record.name, references);
//    output_hits_paf(outstring, nams, record.name, references, index_parameters.k,
//                        record.seq.length());

}


