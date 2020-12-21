import os
import sys

from collections import defaultdict
from time import time
from array import array
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


############### TMP #######
from types import ModuleType, FunctionType
from gc import get_referents

# Custom objects know their class.
# Function objects seem to know way too much, including modules.
# Exclude modules as well.
BLACKLIST = type, ModuleType, FunctionType


def getsize(obj):
    """sum size of object & members."""
    if isinstance(obj, BLACKLIST):
        raise TypeError('getsize() does not take argument of type: '+ str(type(obj)))
    seen_ids = set()
    size = 0
    objects = [obj]
    while objects:
        need_referents = []
        for obj in objects:
            if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size
######################

def import_data(args):
    ref_segment_sequences = help_functions.pickle_load( os.path.join(args.outfolder, 'ref_segment_sequences.pickle') )
    ref_exon_sequences = help_functions.pickle_load( os.path.join(args.outfolder, 'ref_exon_sequences.pickle') )
    ref_flank_sequences = help_functions.pickle_load( os.path.join(args.outfolder, 'ref_flank_sequences.pickle') )
    splices_to_transcripts = help_functions.pickle_load( os.path.join(args.outfolder, 'splices_to_transcripts.pickle') )
    transcripts_to_splices = help_functions.pickle_load( os.path.join(args.outfolder, 'transcripts_to_splices.pickle') )
    all_splice_pairs_annotations = help_functions.pickle_load( os.path.join(args.outfolder, 'all_splice_pairs_annotations.pickle') )
    all_splice_sites_annotations = help_functions.pickle_load( os.path.join(args.outfolder, 'all_splice_sites_annotations.pickle') )
    parts_to_segments = help_functions.pickle_load( os.path.join(args.outfolder, 'parts_to_segments.pickle') )
    segment_to_gene = help_functions.pickle_load( os.path.join(args.outfolder, 'segment_to_gene.pickle') )
    gene_to_small_segments = help_functions.pickle_load( os.path.join(args.outfolder, 'gene_to_small_segments.pickle') )
    max_intron_chr = help_functions.pickle_load( os.path.join(args.outfolder, 'max_intron_chr.pickle') )

    chr_to_id = help_functions.pickle_load( os.path.join(args.outfolder, 'chr_to_id.pickle') )
    id_to_chr = help_functions.pickle_load( os.path.join(args.outfolder, 'id_to_chr.pickle') )

    tiling_ref_segment_sequences = help_functions.pickle_load( os.path.join(args.outfolder, 'tiling_ref_segment_sequences.pickle') )
    tiling_parts_to_segments = help_functions.pickle_load( os.path.join(args.outfolder, 'tiling_parts_to_segments.pickle') )
    tiling_segment_to_gene = help_functions.pickle_load( os.path.join(args.outfolder, 'tiling_segment_to_gene.pickle') )
    tiling_gene_to_small_segments = help_functions.pickle_load( os.path.join(args.outfolder, 'tiling_gene_to_small_segments.pickle') )
    tiling_structures = [tiling_segment_to_gene, tiling_parts_to_segments, tiling_gene_to_small_segments, tiling_ref_segment_sequences]

    # print("ref_segment_sequences:", getsize(ref_segment_sequences)//1000000)
    # print("ref_exon_sequences:", getsize(ref_exon_sequences)//1000000)
    # print("ref_flank_sequences:", getsize(ref_flank_sequences)//1000000)
    # print("splices_to_transcripts:", getsize(splices_to_transcripts)//1000000)
    # print("transcripts_to_splices:", getsize(transcripts_to_splices)//1000000)
    # print("all_splice_pairs_annotations:", getsize(all_splice_pairs_annotations)//1000000)
    # print("all_splice_sites_annotations:", getsize(all_splice_sites_annotations)//1000000)
    # print("parts_to_segments:", getsize(parts_to_segments)//1000000)
    # print("segment_to_gene:", getsize(segment_to_gene)//1000000)
    # print("gene_to_small_segments:", getsize(gene_to_small_segments)//1000000)
    # print("max_intron_chr:", getsize(max_intron_chr)//1000000)
    # print("chr_to_id:", getsize(chr_to_id)//1000000)
    # print("id_to_chr:", getsize(id_to_chr)//1000000)

    # print("tiling_ref_segment_sequences:", getsize(tiling_ref_segment_sequences)//1000000)
    # print("tiling_parts_to_segments:", getsize(tiling_parts_to_segments)//1000000)
    # print("tiling_segment_to_gene:", getsize(tiling_segment_to_gene)//1000000)
    # print("tiling_gene_to_small_segments:", getsize(tiling_gene_to_small_segments)//1000000)

    return ref_segment_sequences, ref_flank_sequences, splices_to_transcripts, \
            transcripts_to_splices, all_splice_pairs_annotations, \
            all_splice_sites_annotations, parts_to_segments,\
            segment_to_gene, gene_to_small_segments, max_intron_chr, \
            ref_exon_sequences, chr_to_id, id_to_chr, tiling_structures


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


def find_exons(chr_id, mam_solution, ref_exon_sequences, ref_segment_sequences, \
                 ref_flank_sequences, all_splice_pairs_annotations):
    # print("BEFO:", [(m.x,m.y) for m in mam_solution])
    valid_introns_sites = all_splice_pairs_annotations[chr_id]
    # print("valid_introns_sites", valid_introns_sites)
    # identify known intron sites
    parts = []
    prev_split_index = 0
    if len(mam_solution) > 1:
        for i, (m1,m2) in enumerate(zip(mam_solution[:-1], mam_solution[1:])):
            if (m1.y, m2.x) in valid_introns_sites:
                parts.append( mam_solution[prev_split_index:i+1] )
                prev_split_index = i+1
        if mam_solution[prev_split_index:]:
            parts.append(mam_solution[prev_split_index:])
        # for p in parts:
        #     print(p)
    else:
        parts.append(mam_solution)
    
    # print("parts",parts)
    # print()
    # identify likely intron sites too large to be a deletion
    parts2 = []
    for part in parts:
        prev_split_index = 0
        if len(part) > 1:
            for i, (m1,m2) in enumerate(zip(part[:-1], part[1:])):
                if m2.x - m1.y > 10:
                    parts2.append( part[prev_split_index:i+1] )
                    prev_split_index = i+1
            if part[prev_split_index:]:
                # print("GAWWWWD", part[prev_split_index:])
                parts2.append(part[prev_split_index:])
        else:
            parts2.append(part)

    # for p in parts2:
    #     print("parts", p)

    # finally look for exact exons or enclosing exons if found
    exons = []
    for part in parts2:
        # check if all segments are adjacent 
        is_adjacent = False
        if len(part) == 1:
            is_adjacent = True
        else: 
            is_adjacent = all([ m1.y == m2.x for (m1,m2) in zip(part[:-1], part[1:]) ])

        if is_adjacent:
            for j, mam in enumerate(part):
                exons.append((mam.x, mam.y, mam.c, mam.d, mam.ref_chr_id))
        else:
            # check if there is an combination of segments and/or adjacent exons over the total region spanned by the exons 

            all_points = []
            segm = {}
            start_points = {}
            end_points = {}
            for m in part:
                all_points.append(m.x)
                all_points.append(m.y)
                segm[(chr_id, m.x, m.y)] = m
                start_points[m.x] = m
                end_points[m.y] = m
            cover = {}
            for i, p1 in enumerate(sorted(all_points)):
                cover[p1] = []
                for j, p2 in enumerate(sorted(all_points)):
                    key = array('L', [chr_id, p1, p2])
                    if key.tobytes() in ref_exon_sequences or  (chr_id, p1, p2) in segm:
                        cover[p1].append(p2)


            # print("all_points", all_points)
            # print("all covers", cover)
            paths = help_functions.find_all_paths(cover, part[0].x, part[-1].y)
            # print(start_points)
            # print(end_points)

            if len(paths) > 0:
                # print("paths[0]" , paths[0])
                path = paths[0]
                for (p1, p2) in zip(path[:-1], path[1:]):
                    if (chr_id, p1, p2) in segm:
                        mam = segm[(chr_id, p1, p2)]
                        exons.append((mam.x, mam.y, mam.c, mam.d, mam.ref_chr_id))
                    else:
                        if p1 in start_points:
                            m = start_points[p1]
                            c = m.c
                            x = m.x
                        else:
                            m = end_points[p1]
                            c = m.d
                            x = m.y 

                        if p2 in start_points:
                            m = start_points[p2]
                            d = m.c
                            y = m.x
                        else:
                            m = end_points[p2]
                            d = m.d
                            y = m.y                            
                        exons.append((x, y, c, d, chr_id))

            else:  #if not, simply give up and put the small spurious intron caused by optimal solution is not containing adjacent segments
                for j, mam in enumerate(part):
                    exons.append((mam.x, mam.y, mam.c, mam.d, mam.ref_chr_id))  

    # print("LOOL NEW",[ (x,y) for x, y, c, d, seq_id  in  exons])

    chained_exon_seqs = []
    predicted_exons = []
    prev_y_coord = -1
    covered = 0
    # for mam in mam_solution:
    for x, y, c, d, seq_id in exons:
        key_tmp = array('L', [seq_id, x, y])
        key = key_tmp.tobytes()
        if key in ref_exon_sequences:
            seq = ref_exon_sequences[key] 
            covered += d - c + 1
        else: 
            if key in ref_segment_sequences:
                seq = ref_segment_sequences[key] 
                covered += d - c + 1
            elif key in ref_flank_sequences:
                seq = ref_flank_sequences[key] 
                covered += d - c + 1
            else:
                print("Bug encountered, {0} is not in {1}".format((x, y), mam_solution))

        if prev_y_coord == x: #adjacent segments means its a flank and we should not add an new exon (i.e., intron split)
            predicted_exons[-1] = (predicted_exons[-1][0], y)  # update the last exon
        else:
            predicted_exons.append( (x, y) )

        prev_y_coord = y
        chained_exon_seqs.append(seq)
    created_ref_seq = "".join([exon for exon in chained_exon_seqs])
    predicted_splices = [ (e1[1],e2[0]) for e1, e2 in zip(predicted_exons[:-1],predicted_exons[1:])]
    return exons, created_ref_seq, predicted_exons, predicted_splices, covered


def get_exact_alignment(read_seq, created_ref_seq, mam_sol_exons_length):
    # if statement: faster and more memory efficient alignments with edlib if sequences are very long. 
    # However, much worse quality.
    if len(created_ref_seq) > 20000 or len(read_seq) > 20000 or (1000 < len(read_seq) < mam_sol_exons_length/10): 
        # print("lenght ref: {0}, length query:{1}".format(len(created_ref_seq), len(read_seq)))
        read_aln, ref_aln, edit_distance = help_functions.edlib_alignment(read_seq, created_ref_seq, aln_mode = "HW")
        match_score = sum([2 for n1,n2 in zip(read_aln, ref_aln) if n1 == n2 ])
        # diff_score = sum([2 for n1,n2 in zip(read_aln, ref_aln) if n1 != n2 ])
        alignment_score = match_score - 2*edit_distance
        # print(read_seq)
    else:
        read_aln, ref_aln, cigar_string, cigar_tuples, alignment_score = help_functions.parasail_alignment(read_seq, created_ref_seq)
        # print(read_acc)
        # print(read_aln)
        # print(ref_aln)

    return read_aln, ref_aln, alignment_score 


def run_tiling_solution(mem_solution, tiling_ref_segment_sequences, ref_flank_sequences, tiling_parts_to_segments, \
                        tiling_segment_to_gene, tiling_gene_to_small_segments, \
                        read_seq, warning_log_file, min_acc, \
                        chr_id, ref_exon_sequences, \
                        all_splice_pairs_annotations,
                        splices_to_transcripts, transcripts_to_splices, \
                        all_splice_sites_annotations, mam_sol_exons_length, \
                        alignment_score, read_aln, ref_aln, current_classification, non_covered_regions, covered,
                        predicted_exons, annotated_to_transcript_id):

    non_covered_regions_tiling, mam_value_tiling, mam_solution_tiling = classify_read_with_mams.main(mem_solution, tiling_ref_segment_sequences, ref_flank_sequences, tiling_parts_to_segments, \
                                                                                                            tiling_segment_to_gene, tiling_gene_to_small_segments, \
                                                                                                            read_seq, warning_log_file, min_acc, is_tiling_instance = True)

    # print("TILING finished Mam solution Tiling!!:",mam_value_tiling, mam_solution_tiling)
    # for zzz2 in mam_solution_tiling:
    #     print(zzz2)
    exons, tiling_created_ref_seq, tiling_predicted_exons, tiling_predicted_splices, tiling_covered = find_exons(chr_id, mam_solution_tiling, ref_exon_sequences, \
                                                                            tiling_ref_segment_sequences, ref_flank_sequences, all_splice_pairs_annotations)

    tiling_classification, tiling_annotated_to_transcript_id = classify_alignment2.main(chr_id, tiling_predicted_splices, splices_to_transcripts, transcripts_to_splices, all_splice_pairs_annotations, all_splice_sites_annotations)

    tiling_read_aln, tiling_ref_aln, tiling_alignment_score = get_exact_alignment(read_seq, tiling_created_ref_seq, mam_sol_exons_length)

    if tiling_alignment_score > alignment_score or tiling_classification == 'FSM':
        # print("TILING")        
        # print("old aln:", read_aln)
        # print("old aln:", ref_aln)
        # print(tiling_classification)            
        # print("BEFORE", non_covered_regions)
        # print("new aln:", tiling_read_aln)
        # print("new aln:", tiling_ref_aln)
        # print("NOW", non_covered_regions_tiling) 
        # print("read was", current_classification, "had aln score", alignment_score)
        # print("tiling read is ", tiling_classification, "has aln score", tiling_alignment_score)
        return  tiling_classification, tiling_alignment_score, non_covered_regions_tiling, \
                tiling_read_aln, tiling_ref_aln, tiling_predicted_exons, tiling_annotated_to_transcript_id, tiling_covered
    else:
        return current_classification, alignment_score, non_covered_regions, \
                read_aln, ref_aln, predicted_exons, annotated_to_transcript_id, covered


def align_single(reads, refs_lengths, args,  batch_number):
    # print("reads:", getsize(reads)//1000000)
    auxillary_data = import_data(args)
    # import time
    # time.sleep(10000)
    mems_path =  os.path.join( args.outfolder, "mummer_mems_batch_{0}.txt".format(batch_number) )
    mems_path_rc =  os.path.join( args.outfolder, "mummer_mems_batch_{0}_rc.txt".format(batch_number) )
    nlog_n_instance_counter = 0
    quadratic_instance_counter = 0
    max_global_intron = args.max_intron
    min_acc = args.min_acc
    ref_segment_sequences, ref_flank_sequences, splices_to_transcripts, \
    transcripts_to_splices, all_splice_pairs_annotations, \
    all_splice_sites_annotations, parts_to_segments, \
    segment_to_gene, gene_to_small_segments, max_intron_chr, \
    ref_exon_sequences, chr_to_id, id_to_chr, tiling_structures = auxillary_data
    # print(ref_flank_sequences.keys())
    if batch_number == -1:
        alignment_outfile = pysam.AlignmentFile( os.path.join(args.outfolder, "reads.sam"), "w", reference_names=list([id_to_chr[id_] for id_ in refs_lengths.keys()]), reference_lengths=list(refs_lengths.values()) ) #, template=samfile)
        warning_log_file = open(os.path.join(args.outfolder, "uLTRA.stderr"), "w")

    else:  
        alignment_outfile = pysam.AlignmentFile( os.path.join(args.outfolder, "reads_batch_{0}.sam".format(batch_number)), "w", reference_names=list([id_to_chr[id_] for id_ in refs_lengths.keys()]), reference_lengths=list(refs_lengths.values()) ) #, template=samfile)
        warning_log_file = open(os.path.join(args.outfolder, "uLTRA_batch_{0}.stderr".format(batch_number)), "w")



    tiling_segment_to_gene, \
    tiling_parts_to_segments, tiling_gene_to_small_segments, \
    tiling_ref_segment_sequences = tiling_structures # unpacking tiling structures

    classifications = defaultdict(str)
    read_accessions_with_mappings = set()
    processed_read_counter = 0

    for (read_acc, mems), (_, mems_rc) in zip(mem_wrapper.get_mem_records(mems_path,reads), mem_wrapper.get_mem_records(mems_path_rc, reads)):
        # multiple = False
        if read_acc not in reads: # if parallelization not all reads in mummer file are in read batches
            continue
        else:
            read_seq_mod = help_functions.remove_read_polyA_ends(reads[read_acc], args.reduce_read_ployA, 1)
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
                # print("Breaking for", chr_id, is_rc, upper_bound_cov, "best:", best_solution_value, "read length:", len(read_seq_mod))
                break
            
            # print("Processing", chr_id, is_rc, upper_bound_cov, "best:", best_solution_value, "read length:", len(read_seq_mod))
            # for mem in all_mems_to_chromosome:
            #     print(mem.exon_part_id, mem.x, mem.y, mem.c, mem.d, '\t', mem.val)
            # print(len(all_mems_to_chromosome))


            max_allowed_intron = min(max_intron_chr[chr_id] + 20000, max_global_intron)
            # print("max_allowed_intron", max_allowed_intron, max_intron_chr[chr_id])
            # print("LEM MEMS", len(all_mems_to_chromosome))
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
                #     print(mem.exon_part_id, mem.x, mem.y, mem.c, mem.d, '\t', mem.val, "(sol)")
            # print(all_chainings)
        is_secondary =  False
        is_rc =  False
        if not all_chainings:
            sam_output.main(read_acc, read_seq_mod, '*', 'unaligned', [], '*', '*', '*', alignment_outfile, is_rc, is_secondary, 0)
            continue

        # all_chainings = sorted(all_chainings, key=lambda x: x[2], reverse=True) 
        all_chainings = sorted(all_chainings, key=lambda x:(-x[2],(x[1][-1].y - x[1][0].x))) 
        # if multiple:
        #     print(all_chainings)
        best_chaining_score = all_chainings[0][2]
        read_alignments = []
        mam_solutions = set()
        for i_nr_sol, (chr_id, mem_solution, chaining_score, is_rc) in enumerate(all_chainings):
            # print(i_nr_sol, chr_id, chaining_score)
            # if read_acc == "100:823|c8740d7a-53bd-4690-aa53-de3ebd003d20": #len(all_chainings) > 20: 
            #     print(read_acc, len(all_chainings), chaining_score)
            #     print(all_chainings)
            # for c in all_chainings:
            #     print(c[1][-1].y - c[1][0].x, c[1][0].x, c[1][-1].y)
            #     sys.exit()
            if chaining_score/float(best_chaining_score) < args.dropoff or i_nr_sol >= args.max_loc:
                # print(chr_id, chaining_score, best_chaining_score, "NOT CONSIDERED")
                # print(i_nr_sol, "lol")
                break
            # print(chr_id, chaining_score, best_chaining_score)

            if is_rc:
                read_seq = help_functions.reverse_complement(read_seq_mod)
            else:
                read_seq = read_seq_mod
            # print("mem solution:", is_rc, chaining_score, mem_solution)
            # print()
            # print()
            # print(is_rc, mem_solution[0].x)
            # print(mem_solution)
            # print(read_seq)
            non_covered_regions, mam_value, mam_solution = classify_read_with_mams.main(mem_solution, ref_segment_sequences, ref_flank_sequences, parts_to_segments, \
                                                                                                                    segment_to_gene, gene_to_small_segments, \
                                                                                                                    read_seq, warning_log_file, min_acc)
            
            # We can enter the if statement below because sometimes the MEM chaining finder will 
            # return multiple optimal chainings that lead to the same mam_solution
            if mam_solution in mam_solutions:
                continue

            mam_solutions.add(mam_solution)
            # print("finished Mam solution:",mam_value, mam_solution)
            # for zzz2 in mam_solution:
            #     print(zzz2)
            # print(non_covered_regions)
            mam_sol_exons_length = sum([ mam.y - mam.x for mam in mam_solution])
            # print(max_intron_size)
            if mam_value > 0:
                exons, created_ref_seq, predicted_exons, predicted_splices, covered = find_exons(chr_id, mam_solution, ref_exon_sequences, \
                                                                                            ref_segment_sequences, ref_flank_sequences, all_splice_pairs_annotations)


                classification, annotated_to_transcript_id = classify_alignment2.main(chr_id, predicted_splices, splices_to_transcripts, transcripts_to_splices, all_splice_pairs_annotations, all_splice_sites_annotations)
                largest_intron_size = max([m2.x - m1.y for m1,m2 in zip(mam_solution[:-1], mam_solution[1:]) ]) if len(mam_solution) > 1 else 0
                if largest_intron_size > max_allowed_intron and classification != 'FSM':
                    # print(i_nr_sol, mem_solution[0].x, "HERE", largest_intron_size)
                    # print()
                    # print(read_acc)
                    # print(classification, alignment_score/len(read_seq), "Score: {0}, old T: {1}, new T: {2}".format(alignment_score, 2*args.alignment_threshold*len(read_seq), len(read_seq)*8*args.alignment_threshold))
                    continue

                read_aln, ref_aln, alignment_score = get_exact_alignment(read_seq, created_ref_seq, mam_sol_exons_length)

                # if classification != 'FSM' and classification != 'NO_SPLICE':
                if classification != 'FSM' and len(non_covered_regions) >= 3 and (max(non_covered_regions[1:-1]) > args.non_covered_cutoff):
                # if "NIC" in classification or (classification != 'FSM' and len(non_covered_regions) >= 3 and (max(non_covered_regions[1:-1]) > args.non_covered_cutoff)):
                    classification, alignment_score, non_covered_regions, \
                    read_aln, ref_aln, predicted_exons, annotated_to_transcript_id, covered = run_tiling_solution(mem_solution, tiling_ref_segment_sequences, ref_flank_sequences, tiling_parts_to_segments, \
                                        tiling_segment_to_gene, tiling_gene_to_small_segments, \
                                        read_seq, warning_log_file, min_acc, chr_id, ref_exon_sequences, \
                                        all_splice_pairs_annotations,
                                        splices_to_transcripts, transcripts_to_splices, \
                                        all_splice_sites_annotations, mam_sol_exons_length, alignment_score, \
                                        read_aln, ref_aln, classification, non_covered_regions, covered, predicted_exons, annotated_to_transcript_id)

                # print(read_aln)
                # print(ref_aln)
                # print()
                # if readl_acc == "m64014_190506_005857/153027831/ccs" or read_acc == "m64014_190506_005857/168625399/ccs":
                #     print(classification, alignment_score/len(read_seq), "Score: {0}, T: {1}".format(alignment_score, 2*args.alignment_threshold*len(read_seq)))
                # print("SEQ identity:")
                # print(alignment_score)
                if alignment_score < 2*args.alignment_threshold*len(read_seq) and classification != 'FSM': # match score * aln_threshold
                    # print(i_nr_sol, mem_solution[0].x, "HERE2", chaining_score, alignment_score, 2*args.alignment_threshold*len(read_seq))
                    # print()
                    # print(read_acc)
                    # print(read_aln)
                    # print(ref_aln)
                    # print(classification, alignment_score/len(read_seq), "Score: {0}, T: {1}".format(alignment_score, 2*args.alignment_threshold*len(read_seq)))
                    continue
                # print(classification)
                # checing for internal (splice site) non covered regions
                if len(non_covered_regions) >= 3 and (max(non_covered_regions[1:-1]) > args.non_covered_cutoff):
                    classification = 'Insufficient_junction_coverage_unclassified'
                coverage = covered / float(len(read_seq)) 
                # print(len(read_alignments), read_alignments)
                genome_start = mam_solution[0].x
                genome_stop = mam_solution[-1].y
                read_alignments.append( (alignment_score, genome_start, genome_stop, read_acc, chr_id, classification, predicted_exons, read_aln, ref_aln, annotated_to_transcript_id, is_rc, coverage) )
            
            # else:
            #     print(i_nr_sol, mem_solution[0].x, "OK", chaining_score, mam_value, mem_solution[-1].y - mem_solution[0].x, mem_solution[0].x)    
        # if read_acc == "100:823|c8740d7a-53bd-4690-aa53-de3ebd003d20": #len(all_chainings) > 20: 
        #     print(read_acc, len(all_chainings), chaining_score, chr_id, max_allowed_intron)
        #     print(read_seq)
        #     # print(all_chainings)
        #     for c in all_chainings:
        #         print(is_rc, c[2], c[1][-1].y - c[1][0].x, c[1][0].x, c[1][-1].y)
        #     print("Nr alignments:", len(read_alignments))
        #     for r in read_alignments:
        #         print(r[0],r[1], r[2], r[3],r[4])
        #     # sys.exit()

            # elif 10*len(read_seq) < mam_sol_exons_length:
            #     print("length ref: {0}, length query:{1}".format(mam_sol_exons_length, len(read_seq)))
            #     print(read_acc, "to chr", chr_id)
            #     print(read_seq)

        ##################  Process alignments and decide primary
        if len(read_alignments) == 0:
            sam_output.main(read_acc, read_seq, '*', 'unaligned', [], '*', '*', '*', alignment_outfile, is_rc, is_secondary, 0)
        else:
            # sorted_wrt_alignement_score = sorted(read_alignments, key = lambda x: x[0], reverse = True)
            sorted_wrt_alignement_score = sorted(read_alignments, key = lambda x: (-x[0], (x[2] - x[1]), x[5]))
            # sorted_wrt_alignement_score = sorted(read_alignments, key = lambda x: (-(x[0] - 0.002*(x[2] - x[1])), len(x[5]) ))
            # if len(read_alignments) > 1 and "FSM" in set( [r[5] for r in read_alignments]):
            #     print(read_acc)
            #     for r in sorted_wrt_alignement_score:
            #         print(r[0], r[5], r[1], (r[2] - r[1]))
            #         print(r[7])
            #         print(r[8])
            #     for c in all_chainings:
            #         print(c[0], c[2], c[1][0].x, c[1][-1].y)
            # sorted_wrt_alignement_score = sorted(read_alignments, key = lambda x: (x[5] != "FSM", -x[0], (x[2] - x[1])))
            best_aln_sw_score = sorted_wrt_alignement_score[0][0]
            more_than_one_alignment = True if len(sorted_wrt_alignement_score) > 1 else False

            # if read_acc == "100:823|c8740d7a-53bd-4690-aa53-de3ebd003d20": #len(all_chainings) > 20: 
            #     print()
            #     print( "SORTED", best_aln_sw_score)
            #     for r in sorted_wrt_alignement_score:
            #         print(r[0],r[1], (r[2]-r[1]), r[5],r[6])
            #     sys.exit()


            for i, (alignment_score, genome_start, genome_stop, read_acc, chr_id, classification, predicted_exons, read_aln, ref_aln, annotated_to_transcript_id, is_rc, coverage) in enumerate(sorted_wrt_alignement_score):
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


                sam_output.main(read_acc, read_seq, id_to_chr[chr_id], classification, predicted_exons, read_aln, ref_aln, annotated_to_transcript_id, alignment_outfile, is_rc, is_secondary, map_score, aln_score = alignment_score)
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



def align_parallel(read_data, refs_lengths, args):
    ####### parallelize alignment #########
    # pool = Pool(processes=mp.cpu_count())
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, original_sigint_handler)

    start_multi = time()
    pool = Pool(processes=int(args.nr_cores))
    try:
        res = pool.map_async(align_single_helper, [ (d, refs_lengths, args, i) for i,d in enumerate(read_data)] )
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


