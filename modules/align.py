import os
import sys

from collections import defaultdict
from time import time
import dill as pickle 
import gffutils
import pysam


import signal
from multiprocessing import Pool
import multiprocessing as mp

from modules import colinear_solver 
from modules import help_functions
from modules import classify_read_with_mams
from modules import classify_alignment2
from modules import sam_output
from modules import mummer_wrapper


# def align_single(read_data, auxillary_data, refs_lengths, args,  batch_number):
def align_single(reads, auxillary_data, refs_lengths, args,  batch_number):
    mems_path =  os.path.join( args.outfolder, "mummer_mems.txt" )
    mems_path_rc =  os.path.join( args.outfolder, "mummer_mems_rc.txt" )
    nlog_n_instance_counter = 0
    quadratic_instance_counter = 0
    if batch_number == -1:
        alignment_outfile = pysam.AlignmentFile( os.path.join(args.outfolder, "torkel.sam"), "w", reference_names=list(refs_lengths.keys()), reference_lengths=list(refs_lengths.values()) ) #, template=samfile)
    else:  
        alignment_outfile = pysam.AlignmentFile( os.path.join(args.outfolder, "torkel_batch_{0}.sam".format(batch_number)), "w", reference_names=list(refs_lengths.keys()), reference_lengths=list(refs_lengths.values()) ) #, template=samfile)


    exon_id_to_choordinates, ref_exon_sequences, splices_to_transcripts, transcripts_to_splices, all_splice_pairs_annotations, all_splice_sites_annotations, parts_to_exons = auxillary_data
    classifications = defaultdict(str)
    read_accessions_with_mappings = set()
    processed_read_counter = 0
    for (read_acc, mems), (_, mems_rc) in zip(mummer_wrapper.get_mummer_records(mems_path), mummer_wrapper.get_mummer_records(mems_path_rc)):
        if read_acc not in reads: # if parallelization not all reads in mummer file are in read batches
            continue
        else:
            read_seq = reads[read_acc]
        # print(read_acc, len(mems), len(mems_rc))
    # for curr_index, (read_acc, read_seq, mems, mems_rc) in enumerate(read_data):
        processed_read_counter += 1
        if processed_read_counter % 100 == 0:
            print('Processed {0} reads in batch {1}'.format(processed_read_counter, batch_number))
        # do the chaining here immediately!
        all_chainings = []
        for chr_id, all_mems_to_chromosome in mems.items():
            if len(all_mems_to_chromosome) < 80:
                solution, mem_solution_value, unique = colinear_solver.read_coverage(all_mems_to_chromosome)
                quadratic_instance_counter += 1 
            else:
                solution, mem_solution_value, unique = colinear_solver.n_logn_read_coverage(all_mems_to_chromosome)
                nlog_n_instance_counter += 1
            # assert mem_solution_value == mem_solution_value2
            # if solution2 != solution and unique:
            #     print("BUG", mem_solution_value, mem_solution_value2)
            #     print(solution)
            #     print(solution2)
            #     sys.exit()

            all_chainings.append( (chr_id, solution, mem_solution_value, False) )
        
        for chr_id, all_mems_to_chromosome in mems_rc.items():
            if len(all_mems_to_chromosome) < 80:
                solution, mem_solution_value, unique = colinear_solver.read_coverage(all_mems_to_chromosome)
            else:
                solution, mem_solution_value, unique = colinear_solver.n_logn_read_coverage(all_mems_to_chromosome)
            # assert mem_solution_value == mem_solution_value2
            # if solution2 != solution and unique:
            #     print("BUG")
            #     print(solution)
            #     print(solution2)
            #     sys.exit()

            all_chainings.append( (chr_id, solution, mem_solution_value, True) )

        is_secondary =  False
        is_rc =  False
        if not all_chainings:
            sam_output.main(read_acc, '*', 'unaligned', [], '*', '*', '*', alignment_outfile, is_rc, is_secondary)
            continue

        all_chainings = sorted(all_chainings, key=lambda x: x[2], reverse=True)

        best_chaining_score = all_chainings[0][2]
        for chr_id, mem_solution, chaining_score, is_rc in all_chainings:
            if chaining_score/float(best_chaining_score) < args.dropoff:
                # print(chr_id, chaining_score, best_chaining_score, "NOT CONSIDERED")
                continue
            # print(chr_id, chaining_score, best_chaining_score)

            if is_rc:
                read_seq = help_functions.reverse_complement(read_seq)
            else:
                read_seq = read_seq
            # print("mem solution:", chaining_score, mem_solution)
            non_covered_regions, mam_value, mam_solution, unique_exon_choordinates =  classify_read_with_mams.main(mem_solution, ref_exon_sequences, parts_to_exons, exon_id_to_choordinates, read_seq, args.overlap_threshold, is_rc)
            if mam_value > 0:
                chained_parts_seq = []
                prev_ref_stop = -1
                predicted_exons = []
                for mam in mam_solution:
                    predicted_exons.append( (mam.x, mam.y) )
                    seq = ref_exon_sequences[mam.ref_chr_id][(mam.x, mam.y)] 
                    if mam.x < prev_ref_stop:
                        chained_parts_seq.append(seq[prev_ref_stop - mam.x: ])
                        help_functions.eprint(read_acc,chr_id, 'mem score:', chaining_score,  best_chaining_score, "mam score:",mam_value)
                        help_functions.eprint("Overlapping exons in solution!",  mam.x, prev_ref_stop, mam)
                        help_functions.eprint(mam_solution)
                        # sys.exit()
                    else:
                        chained_parts_seq.append(seq)
                    prev_ref_stop = mam.y

                created_ref_seq = "".join([part for part in chained_parts_seq])
                predicted_splices = [ (e1[1],e2[0]) for e1, e2 in zip(predicted_exons[:-1],predicted_exons[1:])]

                read_aln, ref_aln, cigar_string, cigar_tuples = help_functions.parasail_alignment(read_seq, created_ref_seq)
                # print(read_acc, "alignment to:", chr_id, "best solution val mems:", mem_solution, 'best mam value:', mam_value, 'read length:', len(read_seq), "final_alignment_stats:" )
                # print(read_aln)
                # print(ref_aln)
                classification, annotated_to_transcript_id = classify_alignment2.main(chr_id, predicted_splices, splices_to_transcripts, transcripts_to_splices, all_splice_pairs_annotations, all_splice_sites_annotations)

                # checing for internal (splice site) non covered regions
                if len(non_covered_regions) >= 3 and (max(non_covered_regions[1:-1]) > args.non_covered_cutoff):
                    classification = 'Unclassified_insufficient_junction_coverage'


                classifications[read_acc] = (classification, mam_value / float(len(read_seq)))
                if chaining_score < best_chaining_score: 
                    is_secondary =  True
                else:
                    is_secondary =  False

                sam_output.main(read_acc, chr_id, classification, predicted_exons, read_aln, ref_aln, annotated_to_transcript_id, alignment_outfile, is_rc, is_secondary)
                read_accessions_with_mappings.add(read_acc)
            else:
                if read_acc not in read_accessions_with_mappings:
                    sam_output.main(read_acc, '*', 'unaligned', [], '*', '*', '*', alignment_outfile, is_rc, is_secondary)

    alignment_outfile.close()

    print("Number of instances solved with quadratic collinear chainer solution:", quadratic_instance_counter)
    print("Number of instances solved with n*log n collinear chainer solution:", nlog_n_instance_counter)
    # print(alignment_outfile.filename, dir(alignment_outfile))
    return classifications, alignment_outfile.filename




def align_single_helper(arguments):
    # for v in arguments:
    #     args = v
    # print(len(arguments))
    return align_single(*arguments)
    # print(len(args))
    # read_data, auxillary_data = args[0], args[1]
    # align_single()
    # return 1



def align_parallel(read_data, auxillary_data, refs_lengths, args):
    ####### parallelize alignment #########
    # pool = Pool(processes=mp.cpu_count())
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, original_sigint_handler)
    mp.set_start_method('spawn')
    print(mp.get_context())
    print("Environment set:", mp.get_context())
    print("Using {0} cores.".format(args.nr_cores))
    start_multi = time()
    pool = Pool(processes=int(args.nr_cores))
    # exon_id_to_choordinates, ref_exon_sequences, splices_to_transcripts, transcripts_to_splices, all_splice_pairs_annotations, all_splice_sites_annotations, parts_to_exons = import_data(args)
    # auxillary_data = exon_id_to_choordinates, ref_exon_sequences, splices_to_transcripts, transcripts_to_splices, all_splice_pairs_annotations, all_splice_sites_annotations, parts_to_exons
    # read_data = batch([ (acc, reads[acc], mems[acc], mems_rc[acc], ) for i, acc in enumerate(mems)], batch_size)
    # alignment_outfiles = []
    # for i in range(args.nr_cores):
    #     # alignment_outfile = pysam.AlignmentFile( os.path.join(args.outfolder, "torkel_batch_{0}.sam".format(i)), "w", reference_names=list(refs_lengths.keys()), reference_lengths=list(refs_lengths.values()) ) #, template=samfile)
    #     # alignment_outfiles.append(alignment_outfile)
    #     alignment_outfiles.append(alignment_outfile)
    try:

        res = pool.map_async(align_single_helper, [ (d, auxillary_data, refs_lengths, args, i) for i,d in enumerate(read_data)] )
        results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        pool.terminate()
        sys.exit()
    else:
        pool.close()
    pool.join()

    classifications = defaultdict(str)
    alignment_outfiles = []
    for r, alignment_outfile_name in results:
        alignment_outfiles.append(alignment_outfile_name)
        for acc in r:
            classifications[acc] = r[acc]
        # print(r)
    print("Time elapesd multiprocessing:", time() - start_multi)  
    return classifications, alignment_outfiles


