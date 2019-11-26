import os
import sys

from collections import defaultdict
from time import time
import dill as pickle 
import gffutils
import pysam


import signal
from multiprocessing import Pool
from operator import attrgetter

from modules import colinear_solver 
from modules import help_functions
from modules import classify_read_with_mams
from modules import classify_alignment2
from modules import sam_output
from modules import mummer_wrapper

def annotate_guaranteed_optimal_bound(mems, is_rc):
    # annotate each instance based on total coverage of read
    upper_bound = {}

    for chr_id, all_mems_to_chromosome in mems.items():
        starts = [("start", m.c) for m in all_mems_to_chromosome]
        stops = [("stop", m.d) for m in all_mems_to_chromosome]
        all_starts_stops = sorted(list(starts) + list(stops), key = lambda x: x[1])
        assert all_starts_stops[0][0] == 'start'
        active_start = all_starts_stops[0][1]
        nr_actives = 1
        intervals = []
        for site, pos in all_starts_stops[1:]:
            if nr_actives == 0:
                active_start = pos

            if site == 'stop':
                nr_actives -= 1
            else:
                nr_actives += 1

            if nr_actives == 0:
                assert site == 'stop'
                intervals.append( (active_start, pos) )

            assert nr_actives >= 0

        tot_cov = 0
        for start, stop in intervals:
            tot_cov += stop - start

        upper_bound[chr_id] = (tot_cov, is_rc, all_mems_to_chromosome)

    # for chr_id, all_mems_to_chromosome in mems.items():
    #     # sort on first c then d since sorting is stable in python 3
    #     # mems_read_sorted = sorted(all_mems_to_chromosome, key= lambda x: x.c) 

    #     mems_read_sorted = sorted(all_mems_to_chromosome, key=attrgetter('d', 'c'))
    #     tot_cov = mems_read_sorted[0].val
    #     curr_active_min = mems_read_sorted[0].c
    #     for m1,m2 in zip(mems_read_sorted[:-1], mems_read_sorted[1:]):
    #         # print(m1.d, m2.d, m2.d - m1.d, curr_active_min)
    #         # assert m2.d >= m1.d
    #         if m2.c > m1.d:
    #             tot_cov += m2.val
    #             curr_active_min = m2.c
    #         else:
    #             if m2.c >= curr_active_min:
    #                 tot_cov += m2.d - m1.d
    #             else:
    #                 tot_cov += m2.d - m1.d
    #                 tot_cov += curr_active_min - m2.c
    #                 curr_active_min = m2.c

    #     upper_bound[chr_id] = (tot_cov, is_rc, all_mems_to_chromosome)

    return upper_bound


def align_single(reads, auxillary_data, refs_lengths, args,  batch_number):
    mems_path =  os.path.join( args.outfolder, "mummer_mems.txt" )
    mems_path_rc =  os.path.join( args.outfolder, "mummer_mems_rc.txt" )
    nlog_n_instance_counter = 0
    quadratic_instance_counter = 0
    if batch_number == -1:
        alignment_outfile = pysam.AlignmentFile( os.path.join(args.outfolder, "torkel.sam"), "w", reference_names=list(refs_lengths.keys()), reference_lengths=list(refs_lengths.values()) ) #, template=samfile)
        warning_log_file = open(os.path.join(args.outfolder, "torkel.stderr"), "w")

    else:  
        alignment_outfile = pysam.AlignmentFile( os.path.join(args.outfolder, "torkel_batch_{0}.sam".format(batch_number)), "w", reference_names=list(refs_lengths.keys()), reference_lengths=list(refs_lengths.values()) ) #, template=samfile)
        warning_log_file = open(os.path.join(args.outfolder, "torkel_batch_{0}.stderr".format(batch_number)), "w")

    exon_id_to_choordinates, ref_exon_sequences, splices_to_transcripts, transcripts_to_splices, all_splice_pairs_annotations, all_splice_sites_annotations, parts_to_exons = auxillary_data
    classifications = defaultdict(str)
    read_accessions_with_mappings = set()
    processed_read_counter = 0

    for (read_acc, mems), (_, mems_rc) in zip(mummer_wrapper.get_mummer_records(mems_path,reads), mummer_wrapper.get_mummer_records(mems_path_rc, reads)):
        multiple = False
        if read_acc not in reads: # if parallelization not all reads in mummer file are in read batches
            continue
        else:
            read_seq = reads[read_acc]
        # print("instance sizes fw:", [ (chr_id, len(mm)) for chr_id, mm in mems.items()])
        # print("instance sizes rc:", [ (chr_id, len(mm)) for chr_id, mm in mems_rc.items()])
        upper_bound = annotate_guaranteed_optimal_bound(mems, False)
        upper_bound_rc = annotate_guaranteed_optimal_bound(mems_rc, True)
        # print(upper_bound)
        # print(upper_bound_rc)
        # TODO: Calculate unsorted mem coverage over the read here! Use this coverage to ignore doing chaining for some chromosomes!
        
        # print(read_acc, len(mems), len(mems_rc))
    # for curr_index, (read_acc, read_seq, mems, mems_rc) in enumerate(read_data):
        processed_read_counter += 1
        if processed_read_counter % 5000 == 0:
            print('Processed {0} reads in batch {1}'.format(processed_read_counter, batch_number))
        # do the chaining here immediately!
        all_chainings = []
        best_solution_value = 0
        for chr_id, (upper_bound_cov, is_rc, all_mems_to_chromosome) in sorted(list(upper_bound.items()) + list(upper_bound_rc.items()), key = lambda x: x[1][0], reverse = True ): # mems.items():
            if upper_bound_cov < best_solution_value:
                # print("Breaking for", chr_id, is_rc, upper_bound_cov, "best:", best_solution_value, "read length:", len(read_seq))
                break
            
            # print("Processing", chr_id, is_rc, upper_bound_cov, "best:", best_solution_value, "read length:", len(read_seq))
            if len(all_mems_to_chromosome) < 80:
                solutions, mem_solution_value = colinear_solver.read_coverage(all_mems_to_chromosome)
                quadratic_instance_counter += 1 
            else:
                solutions, mem_solution_value = colinear_solver.n_logn_read_coverage(all_mems_to_chromosome)
                nlog_n_instance_counter += 1

            if mem_solution_value > best_solution_value:
                best_solution_value = mem_solution_value
                # print("best now:", mem_solution_value)

            if len(solutions) > 1:
                # print("More than 1 solution on chromosome")
                # for sol in solutions:
                #     print(mem_solution_value, [ (m.x, m.y) for m in sol])
                multiple = True

            for sol in solutions:
                all_chainings.append( (chr_id, sol, mem_solution_value, is_rc) )

        # sys.exit()

        # for chr_id, all_mems_to_chromosome in mems_rc.items():
        #     if len(all_mems_to_chromosome) < 80:
        #         solutions, mem_solution_value = colinear_solver.read_coverage(all_mems_to_chromosome)
        #     else:
        #         solutions, mem_solution_value = colinear_solver.n_logn_read_coverage(all_mems_to_chromosome)

        #     if len(solutions) > 1:
        #         print("More than 1 solution on chromosome")
        #         for sol in solutions:
        #             print(mem_solution_value, [ (m.x, m.y) for m in sol])
        #         multiple = True

        #     for sol in solutions:
        #         all_chainings.append( (chr_id, sol, mem_solution_value, True) )
        
        # print("Finished solving colinear_solver fw and rv", len(all_mems_to_chromosome) )
        # sys.exit()
        is_secondary =  False
        is_rc =  False
        if not all_chainings:
            sam_output.main(read_acc, '*', 'unaligned', [], '*', '*', '*', alignment_outfile, is_rc, is_secondary, 0)
            continue

        all_chainings = sorted(all_chainings, key=lambda x: x[2], reverse=True)
        # if multiple:
        #     print(all_chainings)
        best_chaining_score = all_chainings[0][2]
        read_alignments = []
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
            non_covered_regions, mam_value, mam_solution, unique_exon_choordinates =  classify_read_with_mams.main(mem_solution, ref_exon_sequences, parts_to_exons, exon_id_to_choordinates, read_seq, args.overlap_threshold, is_rc, warning_log_file)
            # print("finished Mam solution:", mam_solution)
            if mam_value > 0:
                chained_parts_seq = []
                prev_ref_stop = -1
                predicted_exons = []
                for mam in mam_solution:
                    predicted_exons.append( (mam.x, mam.y) )
                    seq = ref_exon_sequences[mam.ref_chr_id][(mam.x, mam.y)] 
                    if mam.x < prev_ref_stop:
                        # chained_parts_seq.append(seq[prev_ref_stop - mam.x: ])
                        warning_log_file.write("Overlapping exons in solution with {0} bases. {1}, {2}, {3}, {4}.\n".format(prev_ref_stop - mam.x, chr_id, mam.x, prev_ref_stop, mam))
                        warning_log_file.write("{0},{1}, mem score: {2}, best mem score:{3}, mam score:{4}\n".format(read_acc, chr_id, chaining_score,  best_chaining_score, mam_value))
                        warning_log_file.write(str(mam_solution) + '\n\n')
                        # sys.exit()
                    # else:
                    chained_parts_seq.append(seq)
                    prev_ref_stop = mam.y
                created_ref_seq = "".join([part for part in chained_parts_seq])
                predicted_splices = [ (e1[1],e2[0]) for e1, e2 in zip(predicted_exons[:-1],predicted_exons[1:])]

                if len(created_ref_seq) > 20000 or len(read_seq) > 20000:
                    print("lenght ref: {0}, length query:{1}".format(len(created_ref_seq), len(read_seq)))
                    read_aln, ref_aln = help_functions.edlib_alignment(read_seq, created_ref_seq)
                else:
                    read_aln, ref_aln, cigar_string, cigar_tuples, alignment_score = help_functions.parasail_alignment(read_seq, created_ref_seq)
                # print(read_acc, "alignment to:", chr_id, "best solution val mems:", mem_solution, 'best mam value:', mam_value, 'read length:', len(read_seq), "final_alignment_stats:" )
                # print(read_aln)
                # print(ref_aln)
                classification, annotated_to_transcript_id = classify_alignment2.main(chr_id, predicted_splices, splices_to_transcripts, transcripts_to_splices, all_splice_pairs_annotations, all_splice_sites_annotations)

                # checing for internal (splice site) non covered regions
                if len(non_covered_regions) >= 3 and (max(non_covered_regions[1:-1]) > args.non_covered_cutoff):
                    classification = 'Unclassified_insufficient_junction_coverage'
                classifications[read_acc] = (classification, mam_value / float(len(read_seq)))

                read_alignments.append( (alignment_score, read_acc, chr_id, classification, predicted_exons, read_aln, ref_aln, annotated_to_transcript_id, is_rc) )


        ##################  Process alignments and decide primary
        if len(read_alignments) == 0:
            sam_output.main(read_acc, '*', 'unaligned', [], '*', '*', '*', alignment_outfile, is_rc, is_secondary, 0)
        else:
            sorted_wrt_alignement_score = sorted(read_alignments, key = lambda x: x[0], reverse = True)
            best_aln_sw_score = sorted_wrt_alignement_score[0][0]
            more_than_one_alignment = True if len(sorted_wrt_alignement_score) > 1 else False
            # if len(sorted_wrt_alignement_score) > 1:
            #     if best_aln_sw_score >  sorted_wrt_alignement_score[1][0]:
            #         map_score = 60
            #     else:                    
            #         map_score = 60

            for i, (alignment_score, read_acc, chr_id, classification, predicted_exons, read_aln, ref_aln, annotated_to_transcript_id, is_rc) in enumerate(sorted_wrt_alignement_score):
                if i == 0:
                    is_secondary =  False
                    assert alignment_score == best_aln_sw_score
                    if more_than_one_alignment:
                        if alignment_score == sorted_wrt_alignement_score[1][0]:
                            map_score = 0
                        else:
                            map_score = 60                           
                    else:
                        map_score = 60
                else:
                    is_secondary =  True
                    map_score = 0


                sam_output.main(read_acc, chr_id, classification, predicted_exons, read_aln, ref_aln, annotated_to_transcript_id, alignment_outfile, is_rc, is_secondary, map_score)
                read_accessions_with_mappings.add(read_acc)


    alignment_outfile.close()
    warning_log_file.close()
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


