#include "cmdline.hpp"

#include <args.hxx>
#include "arguments.hpp"
#include "version.hpp"

class Version {};

CommandLineOptions parse_command_line_arguments(int argc, char **argv) {

    args::ArgumentParser parser("namfinder " + version_string());
    parser.helpParams.showTerminator = false;
    parser.helpParams.helpindent = 20;
    parser.helpParams.width = 90;
    parser.helpParams.programName = "namfinder";
    parser.helpParams.shortSeparator = " ";

    args::HelpFlag help(parser, "help", "Print help and exit", {'h', "help"});
    args::ActionFlag version(parser, "version", "Print version and exit", {"version"}, []() { throw Version(); });

    // Threading
    args::ValueFlag<int> threads(parser, "INT", "Number of threads [3]", {'t', "threads"});
    args::ValueFlag<int> chunk_size(parser, "INT", "Number of reads processed by a worker thread at once [10000]", {"chunk-size"}, args::Options::Hidden);

    args::Group io(parser, "Input/output:");
    args::ValueFlag<std::string> o(parser, "PATH", "redirect output to file [stdout]", {'o'});
    args::Flag v(parser, "v", "Verbose output", {'v'});
    args::Flag S(parser, "S", "Sort output NAMs for each query based on score. Default is to sort first by ref ID, then by query coordinate, then by reference coordinate.", {'S'});

    args::ValueFlag<int> N(parser, "INT", "Retain at most INT secondary alignments (is upper bounded by -M and depends on -S) [0]", {'N'});
    args::ValueFlag<std::string> index_statistics(parser, "PATH", "Print statistics of indexing to PATH", {"index-statistics"});
    args::Flag i(parser, "index", "Do not map reads; only generate the strobemer index and write it to disk. If read files are provided, they are used to estimate read length", {"create-index", 'i'});
    args::Flag use_index(parser, "use_index", "Use a pre-generated index previously written with --create-index.", { "use-index" });

    args::Group seeding_group(parser, "Seeding:");
    auto seeding = SeedingArguments{parser};

    args::Group search(parser, "Search parameters:");

    args::ValueFlag<int> C(parser, "INT", "Mask (do not process) strobemer hits with count larger than C [1000]", {'C'});
    args::ValueFlag<int> L(parser, "INT", "Print at most L NAMs per query [1000]. Will print the NAMs with highest score S = n_strobemer_hits * query_span.", {'L'});

    args::Positional<std::string> ref_filename(parser, "reference", "Reference in FASTA format", args::Options::Required);
    args::Positional<std::string> reads1_filename(parser, "reads1", "Reads 1 in FASTA or FASTQ format, optionally gzip compressed");
    args::Positional<std::string> reads2_filename(parser, "reads2", "Reads 2 in FASTA or FASTQ format, optionally gzip compressed");

    try {
        parser.ParseCLI(argc, argv);
    }
    catch (const args::Completion& e) {
        std::cout << e.what();
        exit(EXIT_SUCCESS);
    }
    catch (const args::Help&) {
        std::cout << parser;
        exit(EXIT_SUCCESS);
    }
    catch (const Version& e) {
        std::cout << version_string() << std::endl;
        exit(EXIT_SUCCESS);
    }
    catch (const args::Error& e) {
        std::cerr << parser;
        std::cerr << "Error: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    CommandLineOptions opt;

    // Threading
    if (threads) { opt.n_threads = args::get(threads); }
    if (chunk_size) { opt.chunk_size = args::get(chunk_size); }

    // Input/output
    if (o) { opt.output_file_name = args::get(o); opt.write_to_stdout = false; }
    if (S) {opt.sort_on_scores = true;}

    if (index_statistics) { opt.logfile_name = args::get(index_statistics); }
//    if (i) { opt.only_gen_index = true; }
    if (use_index) { opt.use_index = true; }

    // Seeding
    if (seeding.k) { opt.k = args::get(seeding.k); opt.k_set = true; }
    if (seeding.l) { opt.l = args::get(seeding.l); opt.l_set = true; }
    if (seeding.u) { opt.u = args::get(seeding.u); opt.u_set = true; }
    if (seeding.s) { opt.s = args::get(seeding.s); opt.s_set = true; }
    if (seeding.m) { opt.max_seed_len = args::get(seeding.m); opt.max_seed_len_set = true; }

//    if (seeding.c) { opt.c = args::get(seeding.c); opt.c_set = true; }



    if (C) { opt.filter_cutoff = args::get(C); }
    if (L) { opt.L = args::get(L); }

    // Reference and read files
    opt.ref_filename = args::get(ref_filename);
    opt.reads_filename1 = args::get(reads1_filename);

      opt.reads_filename2 = std::string();
      opt.is_SE = true;
//    }


    return opt;
}
