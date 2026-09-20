#ifndef CMDLINE_HPP
#define CMDLINE_HPP

#include <vector>
#include <string>
#include <utility>

struct CommandLineOptions {
    int n_threads { 3 };
    int chunk_size { 10000 };

    // Input/output
    std::string output_file_name;
    bool write_to_stdout { true };
    std::string logfile_name { "" };
    bool use_index { false };
    bool sort_on_scores{false};

    // Seeding
    bool max_seed_len_set { false };
    bool k_set { false };
    bool s_set { false };
    bool l_set { false };
    bool u_set { false };
    bool c_set { false };
    int k { 20 };
    int l { 0 };
    int u { 7 };
    int s { 16 };
    int c { 8 };
    int max_seed_len {255};

    int L { 1000 };
    int C { 1000 };
    int filter_cutoff {1000};

    // Reference and read files
    std::string ref_filename; // This is either a fasta file or an index file - if fasta, indexing will be run
    std::string reads_filename1;
    std::string reads_filename2;
    bool is_SE { true };
    bool is_interleaved { false };
};

CommandLineOptions parse_command_line_arguments(int argc, char **argv);

#endif
