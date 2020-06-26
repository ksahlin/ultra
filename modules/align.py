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
from modules import mem_wrapper

def annotate_guaranteed_optimal_bound(mems, is_rc, max_intron_chr, max_global_intron):
    """
        Calculate the maximum coverage (mem-score) that a read can get per chromosome
        and annotate this value to each instance. We can use this annotation to avoid expensive MAM calculation 
        because we can basically continue to next read if the theoretically possibel maximum chaining value 
        for a chromosome is smaller than a solution already computed
    """
    # print("max_intron_chr", max_intron_chr)
    # split mems into separate instances if there is more than max_allowed_intron nt between two consecutive hits!
    # need to reindex j based on splitting the solution per chromosome..
    all_mem_sols = {} 
    for chr_id, all_mems_to_chromosome in mems.items():
        max_allowed_intron = min(max_intron_chr[chr_id] + 20000, max_global_intron)
        # print("max_allowed_intron",chr_id, max_allowed_intron)
        chr_instance_index = 0
        j_reindex = 0
        if len(all_mems_to_chromosome) > 1:
            m1 = all_mems_to_chromosome[0]
            m1 = m1._replace(j = j_reindex)
            # m1.j = j_reindex 
            curr_instance = [m1]
            for m1,m2 in zip(all_mems_to_chromosome[:-1], all_mems_to_chromosome[1:]):
                if m2.x - m1.y > max_allowed_intron:
                    all_mem_sols[(chr_id, chr_instance_index)] = curr_instance
                    chr_instance_index  += 1 
                    j_reindex = 0
                    m2 = m2._replace(j = j_reindex)
                    curr_instance = [m2]

                else:
                    j_reindex += 1 
                    # print(m2)
                    m2 = m2._replace(j = j_reindex)
                    # print(m2)
                    # print()
                    curr_instance.append(m2)
            all_mem_sols[(chr_id, chr_instance_index)] = curr_instance
                       
        else:
            all_mem_sols[(chr_id, chr_instance_index)] = all_mems_to_chromosome


    # Annotate maximum possible solution to each instance
    upper_bound = {}

    for (chr_id, chr_instance_index), all_mems_to_solution in all_mem_sols.items():
        starts = [("start", m.c) for m in all_mems_to_solution]
        stops = [("stop", m.d) for m in all_mems_to_solution]
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
            tot_cov += stop - start  + 1 # MEM choordinates are inclusive

        upper_bound[(chr_id, chr_instance_index)] = (tot_cov, is_rc, all_mems_to_solution)
    # for i in upper_bound:
    #     print("sols:", i, upper_bound[i][0], len(upper_bound[i][2]))

    return upper_bound


def align_single(reads, auxillary_data, refs_lengths, args,  batch_number):
    mems_path =  os.path.join( args.outfolder, "mummer_mems_batch_{0}.txt".format(batch_number) )
    mems_path_rc =  os.path.join( args.outfolder, "mummer_mems_batch_{0}_rc.txt".format(batch_number) )
    nlog_n_instance_counter = 0
    quadratic_instance_counter = 0
    max_global_intron = args.max_intron
    min_acc = args.min_acc
    if batch_number == -1:
        alignment_outfile = pysam.AlignmentFile( os.path.join(args.outfolder, "torkel.sam"), "w", reference_names=list(refs_lengths.keys()), reference_lengths=list(refs_lengths.values()) ) #, template=samfile)
        warning_log_file = open(os.path.join(args.outfolder, "torkel.stderr"), "w")

    else:  
        alignment_outfile = pysam.AlignmentFile( os.path.join(args.outfolder, "torkel_batch_{0}.sam".format(batch_number)), "w", reference_names=list(refs_lengths.keys()), reference_lengths=list(refs_lengths.values()) ) #, template=samfile)
        warning_log_file = open(os.path.join(args.outfolder, "torkel_batch_{0}.stderr".format(batch_number)), "w")

    exon_id_to_choordinates, ref_exon_sequences, ref_flank_sequences, splices_to_transcripts, \
    transcripts_to_splices, all_splice_pairs_annotations, \
    all_splice_sites_annotations, parts_to_exons, \
    exon_to_gene, gene_to_small_exons, max_intron_chr = auxillary_data

    classifications = defaultdict(str)
    read_accessions_with_mappings = set()
    processed_read_counter = 0

    for (read_acc, mems), (_, mems_rc) in zip(mem_wrapper.get_mem_records(mems_path,reads), mem_wrapper.get_mem_records(mems_path_rc, reads)):
        # multiple = False
        if read_acc not in reads: # if parallelization not all reads in mummer file are in read batches
            continue
        else:
            read_seq = help_functions.remove_read_polyA_ends(reads[read_acc], args.reduce_read_ployA, 1)
        # print("instance sizes fw:", [ (chr_id, len(mm)) for chr_id, mm in mems.items()])
        # print("instance sizes rc:", [ (chr_id, len(mm)) for chr_id, mm in mems_rc.items()])
        # print()
        # print()
        # print(read_acc)
        upper_bound = annotate_guaranteed_optimal_bound(mems, False, max_intron_chr, max_global_intron)
        upper_bound_rc = annotate_guaranteed_optimal_bound(mems_rc, True, max_intron_chr, max_global_intron)
        # print()
        processed_read_counter += 1
        if processed_read_counter % 5000 == 0:
            print('Processed {0} reads in batch {1}'.format(processed_read_counter, batch_number))
        # do the chaining here immediately!
        all_chainings = []
        best_solution_value = 0
        for (chr_id, chr_instance_index) , (upper_bound_cov, is_rc, all_mems_to_chromosome) in sorted(list(upper_bound.items()) + list(upper_bound_rc.items()), key = lambda x: x[1][0], reverse = True ): # mems.items():
            # print((chr_id, chr_instance_index), upper_bound_cov, all_mems_to_chromosome)
            if upper_bound_cov < best_solution_value*args.dropoff:
                # print("Breaking for", chr_id, is_rc, upper_bound_cov, "best:", best_solution_value, "read length:", len(read_seq))
                break
            
            # print("Processing", chr_id, is_rc, upper_bound_cov, "best:", best_solution_value, "read length:", len(read_seq))
            # for mem in all_mems_to_chromosome:
            #     print(mem.exon_part_id, mem.x, mem.y, mem.c, mem.d, '\t', mem.val)
            # print(len(all_mems_to_chromosome))


            max_allowed_intron = min(max_intron_chr[chr_id] + 20000, max_global_intron)
            # print("max_allowed_intron", max_allowed_intron, max_intron_chr[chr_id])
            if len(all_mems_to_chromosome) < 90:
                solutions, mem_solution_value = colinear_solver.read_coverage(all_mems_to_chromosome, max_allowed_intron)
                quadratic_instance_counter += 1 
            else:
                solutions, mem_solution_value = colinear_solver.n_logn_read_coverage(all_mems_to_chromosome)
                nlog_n_instance_counter += 1

            if mem_solution_value > best_solution_value:
                best_solution_value = mem_solution_value
                # print("best now:", mem_solution_value)

            # if len(solutions) > 1:
            #     # print("More than 1 solution on chromosome")
            #     # for sol in solutions:
            #     #     print(mem_solution_value, [ (m.x, m.y) for m in sol])
            #     # multiple = True

            for sol in solutions:
                all_chainings.append( (chr_id, sol, mem_solution_value, is_rc) )
                # for mem in sol:
                #     print(mem.exon_part_id, mem.x, mem.y, "(sol)")
            # print(all_chainings)
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
            # print("mem solution:", is_rc, chaining_score, mem_solution)
            non_covered_regions, mam_value, mam_solution = classify_read_with_mams.main(mem_solution, ref_exon_sequences, ref_flank_sequences, parts_to_exons, \
                                                                                                                    exon_id_to_choordinates, exon_to_gene, gene_to_small_exons, \
                                                                                                                    read_seq, warning_log_file, min_acc)
            # print("finished Mam solution:",mam_value, mam_solution)
            # for zzz2 in mam_solution:
            #     print(zzz2)
            # print(non_covered_regions)
            mam_sol_exons_length = sum([ mam.y - mam.x for mam in mam_solution])
            # print(max_intron_size)
            if mam_value > 0:
                chained_exon_seqs = []
                prev_ref_stop = -1
                predicted_exons = []
                tot_exons_len = 0
                prev_y_coord = -1
                covered = 0
                for mam in mam_solution:
                    if (mam.x, mam.y) in ref_exon_sequences[mam.ref_chr_id]:
                        seq = ref_exon_sequences[mam.ref_chr_id][(mam.x, mam.y)] 
                        covered += mam.d - mam.c + 1
                    elif (mam.x, mam.y) in ref_flank_sequences[mam.ref_chr_id]:
                        seq = ref_flank_sequences[mam.ref_chr_id][(mam.x, mam.y)] 
                        covered += mam.d - mam.c + 1

                    else:
                        print("Bug encountered, {0} is not in {1}".format((mam.x, mam.y), mam_solution))

                    if prev_y_coord == mam.x: #adjacent segments means its a flank and we should not add an new exon (i.e., intron split)
                        predicted_exons[-1] = (predicted_exons[-1][0], mam.y)  # update the last exon
                    else:
                        predicted_exons.append( (mam.x, mam.y) )
                    prev_y_coord = mam.y
                    # if mam.x < prev_ref_stop:
                        # chained_exon_seqs.append(seq[prev_ref_stop - mam.x: ])
                        # warning_log_file.write("Overlapping exons in solution with {0} bases. {1}, {2}, {3}, {4}.\n".format(prev_ref_stop - mam.x, chr_id, mam.x, prev_ref_stop, mam))
                        # warning_log_file.write("{0},{1}, mem score: {2}, best mem score:{3}, mam score:{4}\n".format(read_acc, chr_id, chaining_score,  best_chaining_score, mam_value))
                        # warning_log_file.write(str(mam_solution) + '\n\n')
                        # sys.exit()
                    # else:
                    chained_exon_seqs.append(seq)
                    prev_ref_stop = mam.y
                created_ref_seq = "".join([exon for exon in chained_exon_seqs])
                predicted_splices = [ (e1[1],e2[0]) for e1, e2 in zip(predicted_exons[:-1],predicted_exons[1:])]

                classification, annotated_to_transcript_id = classify_alignment2.main(chr_id, predicted_splices, splices_to_transcripts, transcripts_to_splices, all_splice_pairs_annotations, all_splice_sites_annotations)
                largest_intron_size = max([m2.x - m1.y for m1,m2 in zip(mam_solution[:-1], mam_solution[1:]) ]) if len(mam_solution) > 1 else 0
                if largest_intron_size > max_allowed_intron and classification != 'FSM':
                    # print()
                    # print(read_acc)
                    # print(classification, alignment_score/len(read_seq), "Score: {0}, old T: {1}, new T: {2}".format(alignment_score, 2*args.alignment_threshold*len(read_seq), len(read_seq)*8*args.alignment_threshold))
                    continue


                if len(created_ref_seq) > 20000 or len(read_seq) > 20000 or (1000 < len(read_seq) < mam_sol_exons_length/10):
                    # print("lenght ref: {0}, length query:{1}".format(len(created_ref_seq), len(read_seq)))
                    read_aln, ref_aln, edit_distance = help_functions.edlib_alignment(read_seq, created_ref_seq, aln_mode = "HW")
                    match_score = sum([2 for n1,n2 in zip(read_aln, ref_aln) if n1 == n2 ])
                    # diff_score = sum([2 for n1,n2 in zip(read_aln, ref_aln) if n1 != n2 ])
                    alignment_score = match_score - 2*edit_distance
                    # print(read_seq)

                else:
                    read_aln, ref_aln, cigar_string, cigar_tuples, alignment_score = help_functions.parasail_alignment(read_seq, created_ref_seq)
                    # print(read_aln)
                    # print(ref_aln)
                    # print("alignment_score:", alignment_score)

                # matches = sum([1 for n1,n2 in zip(read_aln, ref_aln) if n1 == n2 ])
                # substitutions = sum([1 for n1,n2 in zip(read_aln, ref_aln) if n1 != n2 and n1 != "-" and n2 != "-" ])
                # deletions = sum([1 for n1,n2 in zip(read_aln, ref_aln) if n1 == "-" ])
                # insertions = sum([1 for n1,n2 in zip(read_aln, ref_aln) if n2 == "-" ])

                # print(read_acc, "alignment to:", chr_id, "best solution val mems:", mem_solution, 'best mam value:', mam_value, 'read length:', len(read_seq), "final_alignment_stats:" )
                # print(read_aln)
                # print(ref_aln)
                # print("to chr", chr_id,  alignment_score, 2*args.alignment_threshold*len(read_seq))
                # print("lenght ref: {0}, length query:{1}".format(len(created_ref_seq), len(read_seq)))
                # print(created_ref_seq)
                # print(read_seq)
                # if alignment_score/len(read_seq) < 5:
                #     print()
                #     print(read_acc)
                #     print(read_aln)
                #     print(ref_aln)
                #     print(alignment_score/len(read_seq), "Score: {0}, old T: {1}, new T: {2}".format(alignment_score, 2*args.alignment_threshold*len(read_seq), len(read_seq)*8*args.alignment_threshold))
                
                # if alignment_score < 8*args.alignment_threshold*len(read_seq):
                #     continue

                classification, annotated_to_transcript_id = classify_alignment2.main(chr_id, predicted_splices, splices_to_transcripts, transcripts_to_splices, all_splice_pairs_annotations, all_splice_sites_annotations)
                
                largest_intron_size = max([m2.x - m1.y for m1,m2 in zip(mam_solution[:-1], mam_solution[1:]) ]) if len(mam_solution) > 1 else 0
                if alignment_score < 8*args.alignment_threshold*len(read_seq) and classification != 'FSM':
                    # print()
                    # print(read_acc)
                    # print(read_aln)
                    # print(ref_aln)
                    # print(classification, alignment_score/len(read_seq), "Score: {0}, old T: {1}, new T: {2}".format(alignment_score, 2*args.alignment_threshold*len(read_seq), len(read_seq)*8*args.alignment_threshold))
                    continue
                # print(classification)
                # checing for internal (splice site) non covered regions
                if len(non_covered_regions) >= 3 and (max(non_covered_regions[1:-1]) > args.non_covered_cutoff):
                    classification = 'Insufficient_junction_coverage_unclassified'
                coverage = covered / float(len(read_seq)) 
                read_alignments.append( (alignment_score, read_acc, chr_id, classification, predicted_exons, read_aln, ref_aln, annotated_to_transcript_id, is_rc, coverage) )
                


            # elif 10*len(read_seq) < mam_sol_exons_length:
            #     print("length ref: {0}, length query:{1}".format(mam_sol_exons_length, len(read_seq)))
            #     print(read_acc, "to chr", chr_id)
            #     print(read_seq)

        ##################  Process alignments and decide primary
        if len(read_alignments) == 0:
            sam_output.main(read_acc, '*', 'unaligned', [], '*', '*', '*', alignment_outfile, is_rc, is_secondary, 0)
        else:
            # sorted_wrt_alignement_score = sorted(read_alignments, key = lambda x: x[0], reverse = True)
            sorted_wrt_alignement_score = sorted(read_alignments, key = lambda x: (-x[0], x[3]))
            best_aln_sw_score = sorted_wrt_alignement_score[0][0]
            more_than_one_alignment = True if len(sorted_wrt_alignement_score) > 1 else False
            # if len(sorted_wrt_alignement_score) > 1:
            #     if best_aln_sw_score >  sorted_wrt_alignement_score[1][0]:
            #         map_score = 60
            #     else:                    
            #         map_score = 60

            for i, (alignment_score, read_acc, chr_id, classification, predicted_exons, read_aln, ref_aln, annotated_to_transcript_id, is_rc, coverage) in enumerate(sorted_wrt_alignement_score):
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
                    classifications[read_acc] = (classification, coverage )
                else:
                    is_secondary =  True
                    map_score = 0


                sam_output.main(read_acc, chr_id, classification, predicted_exons, read_aln, ref_aln, annotated_to_transcript_id, alignment_outfile, is_rc, is_secondary, map_score, aln_score = alignment_score)
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


