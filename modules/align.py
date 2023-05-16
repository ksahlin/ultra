import os
import sys


from collections import defaultdict
from time import time
from array import array
import dill as pickle 
import gffutils
import pysam

import signal

from modules import colinear_solver 
from modules import help_functions
from modules import classify_read_with_mams
from modules import classify_alignment2
from modules import sam_output
from modules import seed_wrapper

from collections import namedtuple
mem = namedtuple('Mem', ['x', 'y', 'c', 'd', 'val', 'j', "exon_part_id"])
globals()[mem.__name__] = mem # Global needed for multiprocessing

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

    if args.index:
        index_folder = args.index
    else:
        index_folder = args.outfolder
        
    ref_segment_sequences = help_functions.pickle_load( os.path.join(index_folder, 'ref_segment_sequences.pickle') )
    ref_exon_sequences = help_functions.pickle_load( os.path.join(index_folder, 'ref_exon_sequences.pickle') )
    ref_flank_sequences = help_functions.pickle_load( os.path.join(index_folder, 'ref_flank_sequences.pickle') )
    splices_to_transcripts = help_functions.pickle_load( os.path.join(index_folder, 'splices_to_transcripts.pickle') )
    transcripts_to_splices = help_functions.pickle_load( os.path.join(index_folder, 'transcripts_to_splices.pickle') )
    all_splice_pairs_annotations = help_functions.pickle_load( os.path.join(index_folder, 'all_splice_pairs_annotations.pickle') )
    all_splice_sites_annotations = help_functions.pickle_load( os.path.join(index_folder, 'all_splice_sites_annotations.pickle') )
    parts_to_segments = help_functions.pickle_load( os.path.join(index_folder, 'parts_to_segments.pickle') )
    segment_to_gene = help_functions.pickle_load( os.path.join(index_folder, 'segment_to_gene.pickle') )
    gene_to_small_segments = help_functions.pickle_load( os.path.join(index_folder, 'gene_to_small_segments.pickle') )
    max_intron_chr = help_functions.pickle_load( os.path.join(index_folder, 'max_intron_chr.pickle') )

    chr_to_id = help_functions.pickle_load( os.path.join(index_folder, 'chr_to_id.pickle') )
    id_to_chr = help_functions.pickle_load( os.path.join(index_folder, 'id_to_chr.pickle') )

    # tiling_ref_segment_sequences = help_functions.pickle_load( os.path.join(args.outfolder, 'tiling_ref_segment_sequences.pickle') )
    # tiling_parts_to_segments = help_functions.pickle_load( os.path.join(args.outfolder, 'tiling_parts_to_segments.pickle') )
    # tiling_segment_to_gene = help_functions.pickle_load( os.path.join(args.outfolder, 'tiling_segment_to_gene.pickle') )
    # tiling_gene_to_small_segments = help_functions.pickle_load( os.path.join(args.outfolder, 'tiling_gene_to_small_segments.pickle') )
    # tiling_structures = [tiling_segment_to_gene, tiling_parts_to_segments, tiling_gene_to_small_segments, tiling_ref_segment_sequences]

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

    # return ref_segment_sequences, ref_flank_sequences, splices_to_transcripts, \
    #         transcripts_to_splices, all_splice_pairs_annotations, \
    #         all_splice_sites_annotations, parts_to_segments,\
    #         segment_to_gene, gene_to_small_segments, max_intron_chr, \
    #         ref_exon_sequences, chr_to_id, id_to_chr, tiling_structures

    return ref_segment_sequences, ref_flank_sequences, splices_to_transcripts, \
            transcripts_to_splices, all_splice_pairs_annotations, \
            all_splice_sites_annotations, parts_to_segments,\
            segment_to_gene, gene_to_small_segments, max_intron_chr, \
            ref_exon_sequences, chr_to_id, id_to_chr

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

        if prev_y_coord >= x: #adjacent segments means its a flank and we should not add an new exon (i.e., intron split)
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


def get_mems_from_input(hits):
    # speed this parsing up by providing already sorted
    # hits per ref_id and ref coord from namfinder

    read_mems = defaultdict(list)
    read_mems_rev = defaultdict(list)

    for line in hits:
        vals =  line.split() #11404_11606           1     11405       202
        exon_part_id = vals[0]
        chr_id, ref_coord_start, ref_coord_end = exon_part_id.split('^')
        # chr_id, ref_coord_start, ref_coord_end = ['1', '1', '1'] # CURRENT DUMMY LINE FOR TESTING OUTSIDE ULTRA'S FORMAT
        chr_id = int(chr_id)
        mem_len = int(vals[3])

        mem_ref_exon_part_start = int(vals[1]) - 1 # convert to 0-indexed reference as in python
        mem_read_start = int(vals[2]) - 1
        ref_coord_start = int(ref_coord_start) # has already been 0-indexed when constructing parts
        mem_genome_start = ref_coord_start + mem_ref_exon_part_start
        
        info_tuple = ( mem_genome_start, mem_genome_start + mem_len - 1,
                        mem_read_start, mem_read_start + mem_len - 1, 
                        mem_len, exon_part_id) # however, for MEM length last coordinate is inclusive of the hit in MEM solvers, not as in python end-indexing

        read_mems[chr_id].append(info_tuple)

    for chr_id in list(read_mems.keys()):
        coordinate_sorted_tuples = sorted(read_mems[chr_id], key = lambda x: x[1])
        sorted_mems = [ mem(x,y,c,d,val,j,e_id) for j, (x, y, c, d, val, e_id) in enumerate(coordinate_sorted_tuples) ]
        read_mems[chr_id] = sorted_mems

    return read_mems



def align_single(process_id, input_queue, output_sam_buffer, classification_and_aln_cov, args):

    # set counters
    nlog_n_instance_counter = 0
    quadratic_instance_counter = 0
    max_global_intron = args.max_intron
    min_acc = args.min_acc
    classifications = defaultdict(str)
    processed_read_counter = 0

    # read in index
    auxillary_data = import_data(args)
    ref_segment_sequences, ref_flank_sequences, splices_to_transcripts, \
    transcripts_to_splices, all_splice_pairs_annotations, \
    all_splice_sites_annotations, parts_to_segments, \
    segment_to_gene, gene_to_small_segments, max_intron_chr, \
    ref_exon_sequences, chr_to_id, id_to_chr = auxillary_data

    warning_log_file = open(os.path.join(args.outfolder, "uLTRA_batch_{0}.stderr".format(process_id)), "w")
    
    classification_list = [0, 0, 0, 0, 0, 0, 0, 0] # entries: [aln_cov, 'FSM', 'unaligned', 'NO_SPLICE', 'Insufficient_junction_coverage_unclassified', 'ISM/NIC_known', 'NIC_novel', 'NNC']
    # code: aln_cov (0), 'FSM' (1), 'unaligned' (2), 'NO_SPLICE' (3), 'Insufficient_junction_coverage_unclassified' (4), 'ISM/NIC_known' (5), 'NIC_novel' (6), 'NNC' (7)
    class_to_offset = {'FSM' : 1, 'unaligned' : 2, 'NO_SPLICE' : 3, 'Insufficient_junction_coverage_unclassified' : 4, 'ISM/NIC_known' : 5, 'NIC_novel' : 6 , 'NNC': 7}

    while True:
        batch = input_queue.get()
        # check for stop
        if batch is None:
            # add the signal back for other consumers
            input_queue.put(batch)
            # stop running
            break

        alignments_output = []

        for b in batch[1]:
            (read_acc, seq, hits, hits_rc) = b
            mems = get_mems_from_input(hits)
            mems_rc = get_mems_from_input(hits_rc)
            read_seq_mod = help_functions.remove_read_polyA_ends(seq, args.reduce_read_ployA, 1)

            upper_bound = annotate_guaranteed_optimal_bound(mems, False, max_intron_chr, max_global_intron)
            upper_bound_rc = annotate_guaranteed_optimal_bound(mems_rc, True, max_intron_chr, max_global_intron)
            # print()
            processed_read_counter += 1
            if processed_read_counter % 5000 == 0:
                print('Processed {0} reads in consumer process {1}. Queue size: {2}'.format(processed_read_counter, process_id, input_queue.qsize()))
            # do the chaining here immediately!
            all_chainings = []
            best_solution_value = 0
            for (chr_id, chr_instance_index) , (upper_bound_cov, is_rc, all_mems_to_chromosome) in sorted(list(upper_bound.items()) + list(upper_bound_rc.items()), key = lambda x: x[1][0], reverse = True ): # mems.items():
                if upper_bound_cov < best_solution_value*args.dropoff:
                    break

                max_allowed_intron = min(max_intron_chr[chr_id] + 20000, max_global_intron)

                if len(all_mems_to_chromosome) < 90:
                    solutions, mem_solution_value = colinear_solver.read_coverage(all_mems_to_chromosome, max_allowed_intron)
                    quadratic_instance_counter += 1 
                else:
                    solutions, mem_solution_value = colinear_solver.n_logn_read_coverage(all_mems_to_chromosome)
                    nlog_n_instance_counter += 1

                if mem_solution_value > best_solution_value:
                    best_solution_value = mem_solution_value

                for sol in solutions:
                    all_chainings.append( (chr_id, sol, mem_solution_value, is_rc) )

            is_secondary =  False
            is_rc =  False
            if not all_chainings:
                sam_aln_entry = sam_output.main(read_acc, read_seq_mod, '*', 'unaligned', [], '*', '*', '*', is_rc, is_secondary, 0)
                alignments_output.append(sam_aln_entry)

                continue

            all_chainings = sorted(all_chainings, key=lambda x:(-x[2],(x[1][-1].y - x[1][0].x))) 
            best_chaining_score = all_chainings[0][2]
            read_alignments = []
            mam_solutions = set()
            for i_nr_sol, (chr_id, mem_solution, chaining_score, is_rc) in enumerate(all_chainings):
                if chaining_score/float(best_chaining_score) < args.dropoff or i_nr_sol >= args.max_loc:
                    break

                if is_rc:
                    read_seq = help_functions.reverse_complement(read_seq_mod)
                else:
                    read_seq = read_seq_mod

                non_covered_regions, mam_value, mam_solution = classify_read_with_mams.main(mem_solution, ref_segment_sequences, ref_flank_sequences, parts_to_segments, \
                                                                                                                        segment_to_gene, gene_to_small_segments, \
                                                                                                                        read_seq, warning_log_file, min_acc)
                
                # We can enter the if statement below because sometimes the MEM chaining finder will 
                # return multiple optimal chainings that lead to the same mam_solution
                if mam_solution in mam_solutions:
                    continue

                mam_solutions.add(mam_solution)
                mam_sol_exons_length = sum([ mam.y - mam.x for mam in mam_solution])
                if mam_value > 0:
                    exons, created_ref_seq, predicted_exons, predicted_splices, covered = find_exons(chr_id, mam_solution, ref_exon_sequences, \
                                                                                                ref_segment_sequences, ref_flank_sequences, all_splice_pairs_annotations)


                    classification, annotated_to_transcript_id = classify_alignment2.main(chr_id, predicted_splices, splices_to_transcripts, transcripts_to_splices, all_splice_pairs_annotations, all_splice_sites_annotations)
                    largest_intron_size = max([m2.x - m1.y for m1,m2 in zip(mam_solution[:-1], mam_solution[1:]) ]) if len(mam_solution) > 1 else 0
                    if largest_intron_size > max_allowed_intron and classification != 'FSM':
                        continue

                    read_aln, ref_aln, alignment_score = get_exact_alignment(read_seq, created_ref_seq, mam_sol_exons_length)

                    if alignment_score < 2*args.alignment_threshold*len(read_seq) and classification != 'FSM': # match score * aln_threshold
                        continue

                    # checing for internal (splice site) non covered regions
                    if len(non_covered_regions) >= 3 and (max(non_covered_regions[1:-1]) > args.non_covered_cutoff):
                        classification = 'Insufficient_junction_coverage_unclassified'
                    coverage = covered / float(len(read_seq)) 
                    genome_start = mam_solution[0].x
                    genome_stop = mam_solution[-1].y
                    read_alignments.append( (alignment_score, genome_start, genome_stop, read_acc, chr_id, classification, predicted_exons, read_aln, ref_aln, annotated_to_transcript_id, is_rc, coverage) )


            ##################  Process alignments and decide primary
            if len(read_alignments) == 0:
                sam_aln_entry = sam_output.main(read_acc, read_seq, '*', 'unaligned', [], '*', '*', '*', is_rc, is_secondary, 0)
                alignments_output.append(sam_aln_entry)
            else:
                sorted_wrt_alignement_score = sorted(read_alignments, key = lambda x: (-x[0], (x[2] - x[1]), x[5]))
                best_aln_sw_score = sorted_wrt_alignement_score[0][0]
                more_than_one_alignment = True if len(sorted_wrt_alignement_score) > 1 else False

                for i, (alignment_score, genome_start, genome_stop, read_acc, chr_id, classification, predicted_exons, read_aln, ref_aln, annotated_to_transcript_id, is_rc, coverage) in enumerate(sorted_wrt_alignement_score):
                    if i == 0:
                        is_secondary =  False
                        assert alignment_score == best_aln_sw_score
                        if more_than_one_alignment:
                            if alignment_score == sorted_wrt_alignement_score[1][0]:
                                map_score = 0
                            else:
                                map_score = 60  # TODO: predict this value with formula instead.                           
                        else:
                            map_score = 60
                        # classifications[read_acc] = (classification, coverage )
                        classification_list[0] += coverage
                        classification_list[class_to_offset[classification]] += 1
                    else:
                        is_secondary =  True
                        map_score = 0


                    sam_aln_entry = sam_output.main(read_acc, read_seq, id_to_chr[chr_id], classification, predicted_exons, read_aln, ref_aln, annotated_to_transcript_id, is_rc, is_secondary, map_score, aln_score = alignment_score)
                    alignments_output.append(sam_aln_entry)
        
        output_sam_buffer.put(alignments_output)

    classification_and_aln_cov.put(classification_list)
    print('Process:', classification_list)

    warning_log_file.close()
    print("Number of instances solved with quadratic collinear chainer solution:", quadratic_instance_counter)
    print("Number of instances solved with n*log n collinear chainer solution:", nlog_n_instance_counter)
    # return classifications




# def align_single_helper(arguments):
#     return align_single(*arguments)




# def align_parallel(read_data, args):
#     ####### parallelize alignment #########
#     # pool = Pool(processes=mp.cpu_count())
#     original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
#     signal.signal(signal.SIGINT, original_sigint_handler)

#     start_multi = time()
#     pool = Pool(processes=int(args.nr_cores))
#     try:
#         res = pool.map_async(align_single_helper, [ (d, args, i) for i,d in enumerate(read_data)] )
#         results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
#     except KeyboardInterrupt:
#         print("Caught KeyboardInterrupt, terminating workers")
#         pool.terminate()
#         sys.exit()
#     else:
#         pool.close()
#     pool.join()

#     # classifications = defaultdict(str)
#     # for r in results:
#     #     for acc in r:
#     #         classifications[acc] = r[acc]
#     #     # print(r)
#     # print("Time elapesd multiprocessing:", time() - start_multi)  
#     # return classifications


