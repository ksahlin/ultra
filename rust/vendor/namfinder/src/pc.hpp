#ifndef pc_hpp
#define pc_hpp

#include <thread>
#include <condition_variable>
#include <mutex>
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <queue>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <optional>

#include "index.hpp"
#include "aln.hpp"
#include "refs.hpp"
#include "fastq.hpp"

class InputBuffer {

public:

    InputBuffer(std::string fname1, int chunk_size)
    : ks1(open_fastq(fname1)),
    chunk_size(chunk_size) { }

    std::mutex mtx;

    input_stream_t ks1;
    input_stream_t ks2;
    std::optional<klibpp::KSeq> lookahead1;
    bool finished_reading{false};
    int chunk_size;
    size_t chunk_index{0};
    bool is_interleaved{false};

    void rewind_reset();
    size_t read_records(std::vector<klibpp::KSeq> &records3,
            int read_count=-1);
};


class OutputBuffer {

public:
    OutputBuffer(std::ostream& out) : out(out) { }

    std::mutex mtx;
    std::ostream &out;
    std::unordered_map<size_t, std::string> chunks;
    size_t next_chunk_index{0};

    void output_records(std::string chunk, size_t chunk_index);
};


void perform_task(InputBuffer &input_buffer, OutputBuffer &output_buffer,
                  AlignmentStatistics& statistics, int& done,
                  const mapping_params &map_param, const IndexParameters& index_parameters, const References& references, const StrobemerIndex& index);

#endif
