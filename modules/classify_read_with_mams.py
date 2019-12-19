
import sys
import re
import math
from itertools import groupby 
# import parasail
import edlib

from collections import namedtuple, defaultdict

from modules import colinear_solver 
from modules import help_functions


mam = namedtuple('Mam', ['x', 'y', 'c', 'd', 'val', "min_segment_length", "exon_id", "ref_chr_id"])
globals()[mam.__name__] = mam # Global needed for multiprocessing


def cigar_to_seq(cigar, query, ref):
    cigar_tuples = []
    result = re.split(r'[=DXSMI]+', cigar)
    i = 0
    for length in result[:-1]:
        i += len(length)
        type_ = cigar[i]
        i += 1
        cigar_tuples.append((int(length), type_ ))

    r_index = 0
    q_index = 0
    q_aln = []
    r_aln = []
    for length_ , type_ in cigar_tuples:
        if type_ == "=" or type_ == "X":
            q_aln.append(query[q_index : q_index + length_])
            r_aln.append(ref[r_index : r_index + length_])

            r_index += length_
            q_index += length_
        
        elif  type_ == "I":
            # insertion w.r.t. reference
            r_aln.append('-' * length_)
            q_aln.append(query[q_index: q_index + length_])
            #  only query index change
            q_index += length_

        elif type_ == 'D':
            # deletion w.r.t. reference
            r_aln.append(ref[r_index: r_index + length_])
            q_aln.append('-' * length_)
            #  only ref index change
            r_index += length_
        
        else:
            print("error")
            print(cigar)
            sys.exit()

    return  "".join([s for s in q_aln]), "".join([s for s in r_aln]), cigar_tuples



def edlib_alignment(query, target, mode = "HW", task = 'locations', k=-1):
    result = edlib.align(query, target, task=task, mode=mode, k=k)
    if result['editDistance'] == -1:
        return [0,0], -1
    
    if task == 'path':
        locations = result['locations']
        ref_start, ref_stop = locations[0][0], locations[0][1]
        cigar_string = result["cigar"]
        query_alignment, target_alignment, cigar_tuples = cigar_to_seq(cigar_string, query, target[ref_start: ref_stop+1 ])
    
    # print(locations)
    # print(query_alignment)
    # print(target_alignment)
    # print(result['editDistance'])

    # query_rc = reverse_complement(query)
    # target_rc = reverse_complement(target)
    # result = edlib.align(query_rc, target_rc, task="path", mode="NW")
    # cigar_string = result["cigar"]
    # locations = result['locations']
    # query_rc_alignment, target_rc_alignment, cigar_tuples = cigar_to_seq(cigar_string, query_rc, target_rc)
    # print(locations)
    # print(query_rc_alignment)
    # print(target_rc_alignment)
    # print(result['editDistance'])

    return result['locations'], result['editDistance'] #, query_alignment, target_alignment


def calc_complessed_score(read_alignment, ref_alignment, m, n):
    """
        Raw score: R = aI +  bX - cO -  dG
        lambda=1.37 and K=0.711
        E=mn2**-S
    """
    states = ['I' if n1 == n2 else 'G' if n1 == '-' or n2 == '-' else 'X' for n1, n2 in zip(read_alignment,ref_alignment) ]
    compressed_profile = [ (element, len(list(i))) for element, i in groupby(states)] 
    print(compressed_profile)
    # return evalue


def calc_evalue(read_alignment, ref_alignment, m, n):
    """
        Raw score: R = aI +  bX - cO -  dG
        lambda=1.37 and K=0.711
        E=mn2**-S
    """
    a, b, c, d = 1, -1, 1,  1
    lambda_=1.37
    K=0.711

    states = ['I' if n1 == n2 else 'G' if n1 == '-' or n2 == '-' else 'X' for n1, n2 in zip(read_alignment,ref_alignment) ]
    I = states.count('I')
    X = states.count('X')
    G = states.count('G')
    O =  len([s1 for s1, s2 in zip(states[:-1],states[1:]) if s1 != s2 and s2 == 'G'])
    if states[0] == 'G': # starts with a gap
        O += 1
    raw_score = a*I +  b*X - c*O -  d*G
    if raw_score < 0:
        raw_score = 0
    bit_score = (lambda_*raw_score - math.log(K) )/ math.log(2)
    evalue = m*n*2**(-bit_score)
    # print(read_alignment)
    # print(ref_alignment)
    print(I,X,G,O)
    print(raw_score, bit_score, evalue)
    return evalue

def contains(sub, pri):
    M, N = len(pri), len(sub)
    i, LAST = 0, M-N+1
    while True:
        try:
            found = pri.index(sub[0], i, LAST) # find first elem in sub
        except ValueError:
            return False
        if pri[found:found+N] == sub:
            return True
        else:
            i = found+1

def main(solution, ref_exon_sequences, parts_to_exons, exon_id_to_choordinates, exon_to_gene, gene_to_small_exons, read_seq, overlap_threshold, is_rc, warning_log_file):
    # chained_parts_seq = []
    # chained_parts_ids = []
    prev_ref_stop = -1
    predicted_transcript = []
    predicted_exons = []
    covered_regions = []
    mam_instance = []
    unique_part_locations = set()
    for mem in solution:
        ref_chr_id, ref_start, ref_stop =  mem.exon_part_id.split('_')
        ref_start, ref_stop = int(ref_start), int(ref_stop)
        # get all exons associated with the part
        # print(parts_to_exons)
        unique_part_locations.add((ref_chr_id, ref_start, ref_stop))
    # print()
    # print(unique_part_locations)
    # print()


    ### FASTER METHOD COMPARED TO OLD VERSION: ONLY ALIGN TO READ ONCE PER UNIQUE EXON LOCATION #####
    ###################################################################################################

    # compress unique exons to only do alignment once 
    unique_exon_choordinates =  defaultdict(set)
    for (ref_chr_id, ref_start, ref_stop) in unique_part_locations:
        exon_ids = parts_to_exons[ref_chr_id][(ref_start, ref_stop)]
        for exon_id in exon_ids:
            e_start, e_stop = exon_id_to_choordinates[exon_id]
            unique_exon_choordinates[ (ref_chr_id, e_start, e_stop) ].add(exon_id)

        # also add all small exons that may be smaller than minimum MEM size
        unique_genes = set(gene_id for exon_id in exon_ids for gene_id in exon_to_gene[exon_id])
        small_exons = set(small_exon_id for gene_id in unique_genes for small_exon_id in gene_to_small_exons[gene_id]) 
        for small_exon_id in small_exons:
            e_start, e_stop = exon_id_to_choordinates[small_exon_id]
            if (ref_chr_id,e_start, e_stop) not in unique_exon_choordinates:
                # print("adding small exon,", e_stop - e_start)
                unique_exon_choordinates[ (ref_chr_id, e_start, e_stop) ].add(small_exon_id)



    # print()
    # print('unique_exon_choordinates', unique_exon_choordinates)
    # print()
    # sys.exit()

    for (ref_chr_id, e_start, e_stop), all_exon_ids in sorted(unique_exon_choordinates.items(), key=lambda x: x[0][1]):
        if e_stop - e_start >= 5:
            # if is_rc:
            #     ref_seq = help_functions.reverse_complement(refs[ref_chr_id][e_start: e_stop])
            # else:
            # ref_seq = refs[ref_chr_id][e_start: e_stop]
            ref_seq = ref_exon_sequences[ref_chr_id][(e_start, e_stop)]
            # print(ref_seq == ref_seq2)
            # assert ref_seq == ref_seq2
            # print(exon_id, e_stop - e_start)
            # align them to the read and get the best approxinate match
            if e_stop - e_start >= 9:
                locations, edit_distance = edlib_alignment(ref_seq, read_seq, mode="HW", k = 0.4*min(len(read_seq), len(ref_seq)) )
                if edit_distance >= 0:
                    # calc_complessed_score(read_alignment, ref_alignment, len(read_seq), len(ref_seq))
                    # e_score = calc_evalue(read_alignment, ref_alignment, len(read_seq), len(ref_seq))
                    start, stop = locations[0]
                    min_segment_length = (stop - start)
                    score = min_segment_length - edit_distance
                    if (min_segment_length - edit_distance)/float(min_segment_length) > 0.6:
                        start, stop = locations[0]
                        covered_regions.append((start,stop, score, exon_id, ref_chr_id))
                        for exon_id in all_exon_ids: break # only need one of the redundant exon_ids
                        mam_tuple = mam(e_start, e_stop, start, stop, 
                                score, min_segment_length,  exon_id, ref_chr_id) 
                        mam_instance.append(mam_tuple)
            
            else: # small exons between 5-9bp needs exact match otherwise too much noise
                locations, edit_distance = edlib_alignment(ref_seq, read_seq, mode="HW", k = 0 )
                # print("HEEERE", ref_seq, e_start, e_stop,ref_chr_id)
                if edit_distance == 0:
                    # print("perfect matches:",ref_seq, locations)
                    score = len(ref_seq)
                    # calc_complessed_score(read_alignment, ref_alignment, len(read_seq), len(ref_seq))
                    # e_score = calc_evalue(read_alignment, ref_alignment, len(read_seq), len(ref_seq))
                    for exon_id in all_exon_ids: break # only need one of the redundant exon_ids
                    for start, stop in locations:
                        covered_regions.append((start,stop, score, exon_id, ref_chr_id))
                        mam_tuple = mam(e_start, e_stop, start, stop, 
                                score, score,  exon_id, ref_chr_id) 
                        mam_instance.append(mam_tuple)

        else:
            pass
            # warning_log_file.write("not aligning exons smaller than 5bp: {0}, {1}, {2}, {3}.\n".format(ref_chr_id, e_start, e_stop, ref_exon_sequences[ref_chr_id][(e_start, e_stop)])) # TODO: align these and take all locations

        if  e_stop - e_start >= 0.8*len(read_seq): # read is potentially contained within exon 
            # print()
            # print("aligning read to exon")
            locations, edit_distance = edlib_alignment(read_seq, ref_seq, mode="HW", k = 0.4*min(len(read_seq), len(ref_seq)) )
            # print(exon_id, e_start, e_stop, len(ref_seq), len(read_seq),edit_distance)
            # print()
            if edit_distance >= 0:
                # min_segment_length = min( len(ref_seq) ,len(read_seq) )
                # score = min_segment_length - edit_distance #/len(read_seq)
                
                start, stop = locations[0]
                min_segment_length = (stop - start)
                score = min_segment_length -  edit_distance #/len(read_seq)
                # print("LOOK:", min_segment_length, edit_distance, score, locations)
                # if e_score < 1.0:
                if (min_segment_length -  edit_distance)/float(min_segment_length) > 0.6:
                    start, stop = 0, len(read_seq) - 1
                    covered_regions.append((start,stop, score, exon_id, ref_chr_id))
                    # for exon_id in all_exon_ids:
                    #     mam_tuple = mam(e_start, e_stop, start, stop, 
                    #             score, min_segment_length,  exon_id, ref_chr_id)
                    #     mam_instance.append(mam_tuple)
                    for exon_id in all_exon_ids: break
                    mam_tuple = mam(e_start, e_stop, start, stop, 
                            score, min_segment_length,  exon_id, ref_chr_id)
                    mam_instance.append(mam_tuple)
    ###################################################################################################
    ###################################################################################################
    ###################################################################################################



    ################################### OLD METHOD ###################################################
    ###################################################################################################

    # for (ref_chr_id, ref_start, ref_stop) in unique_part_locations:
    #     exon_ids = parts_to_exons[ref_chr_id][(ref_start, ref_stop)]
    #     print()
    #     print((ref_chr_id, ref_start, ref_stop))
    #     print(exon_ids)
    #     for exon_id in exon_ids:
    #         e_start, e_stop = exon_id_to_choordinates[exon_id]
    #         if e_stop - e_start > 10:
    #             # if is_rc:
    #             #     ref_seq = help_functions.reverse_complement(refs[ref_chr_id][e_start: e_stop])
    #             # else:
    #             ref_seq = refs[ref_chr_id][e_start: e_stop]
    #             # print(exon_id, e_stop - e_start)
    #             # align them to the read and get the best approxinate match

    #             locations, edit_distance, read_alignment, ref_alignment = edlib_alignment(ref_seq, read_seq)
    #             print(exon_id, e_start, e_stop, len(ref_seq), len(read_seq),edit_distance,(len(ref_seq) - edit_distance)/len(ref_seq))
    #             if edit_distance >= 0:
    #                 calc_complessed_score(read_alignment, ref_alignment, len(read_seq), len(ref_seq))
    #                 # e_score = calc_evalue(read_alignment, ref_alignment, len(read_seq), len(ref_seq))
    #                 score = (len(ref_seq) - edit_distance)/len(ref_seq)
    #                 # print(e_score, score, edit_distance )

    #                 if score > 0.6:
    #                     start, stop = locations[0]
    #                     covered_regions.append((start,stop, score, exon_id, ref_chr_id))
    #                     mam_tuple = mam(e_start, e_stop, start, stop, 
    #                             (stop - start + 1)*score, score,  exon_id, ref_chr_id)
    #                     mam_instance.append(mam_tuple)
    #         else:
    #             print("not aligning exons smaller than 10bp") # TODO: align these and take all locations

    #         if  e_stop - e_start >= 0.8*len(read_seq): # read is potentially contained within exon 
    #             print()
    #             print("aligning read to exon")
    #             locations, edit_distance, ref_alignment, read_alignment = edlib_alignment(read_seq, ref_seq)
    #             print(exon_id, e_start, e_stop, len(ref_seq), len(read_seq),edit_distance)
    #             print()
    #             if edit_distance >= 0:
    #                 score = (len(read_seq) - edit_distance)/len(read_seq)

    #                 if score > 0.6:
    #                     start, stop = 0, len(read_seq) - 1
    #                     covered_regions.append((start,stop, score, exon_id, ref_chr_id))
    #                     mam_tuple = mam(e_start, e_stop, start, stop, 
    #                             (stop - start + 1)*score, score,  exon_id, ref_chr_id)
    #                     mam_instance.append(mam_tuple)
    
    ###################################################################################################
    ###################################################################################################
    ###################################################################################################
      
    # print('mam instance')
    # help_functions.eprint(sorted(mam_instance, key= lambda x: x.x))
    # best find a colinear chaining with best cover (score of each segment is length of mam*identity)
    if mam_instance:
        mam_solution, value, unique = colinear_solver.read_coverage_mam_score(mam_instance, overlap_threshold)
    else:
        return 'None', -1, 0, []
    # classify read here!! 
    # pay attention to that if a lolution where a read has been aligned to a exon is accepted, 
    # its probably contained within the exon and should not be counted as FSM

    # output mam_solution per exon
    # print(covered_regions)
    # print("Solutiuon mam")
    # print(len(read_seq))
    # print("val", value, "the mam_solution:", mam_solution)
    covered = sum([mam.d-mam.c + 1 for mam in mam_solution])
    if len(mam_solution) > 0:
        non_covered_regions = []
        non_covered_regions.append( mam_solution[0].c )
        if len(mam_solution) > 1:
            for mam1, mam2 in zip(mam_solution[:-1],mam_solution[1:]):
                non_covered_regions.append( mam2.c - mam1.d -1 )
            # non_covered_regions = [mam2.c-mam1.d for mam1, mam2 in zip(mam_solution[:-1],mam_solution[1:])]
        # add beginning and end
        non_covered_regions.append( len(read_seq)  - mam_solution[-1].d )


    else:
        non_covered_regions = []

    # print(len(read_seq), covered)
    # print(non_covered_regions)
    # if max(non_covered_regions) > 10:
    #     print(non_covered_regions)
    #     sys.exit()
    # read_coverage = 
    # print()
    # print()
    # print()
    return non_covered_regions, value, mam_solution, unique_exon_choordinates
    # check if exon combination is present in some annotation


    #     ref_start, ref_stop = int(ref_start), int(ref_stop)
    #     predicted_transcript.append( (ref_start, ref_stop) )
    #     # print("lll",ref_chr_id, ref_start, ref_stop)
    #     seq = refs_sequences[ref_chr_id][(ref_start, ref_stop)]
    #     if ref_start < prev_ref_stop:
    #         chained_parts_seq.append(seq[prev_ref_stop - ref_start: ])
    #     else:
    #         chained_parts_seq.append(seq)

    #     if ref_start > prev_ref_stop + 1 : # new separate exon on the genome
    #         chained_parts_ids.append(parts_to_exons[ref_chr_id][ (ref_start, ref_stop)])
    #         predicted_exons.append((prev_ref_stop, ref_start))
    #     else: # part of an exon
    #         chained_parts_ids[-1] = chained_parts_ids[-1] ^ parts_to_exons[ref_chr_id][(ref_start, ref_stop)]  #here we get mora annotations than we want , need a function to check that an exon has all parts covered!

    #     prev_ref_stop = ref_stop

    # predicted_exons.append((prev_ref_stop, -1))

    # created_ref_seq = "".join([part for part in chained_parts_seq])
    # # predicted_transcript = tuple( mem.exon_part_id for mem in solution)
    # # print( "predicted:", predicted_transcript)
    # print(predicted_exons)
    # predicted_exons = [ (e1[1],e2[0]) for (e1, e2) in zip(predicted_exons[:-1], predicted_exons[1:] )]
    # predicted_splices = [ (e1[1],e2[0]) for e1, e2 in zip(predicted_exons[:-1],predicted_exons[1:])]

    # print(predicted_exons)
    # print(predicted_splices)


    # read_aln, ref_aln, cigar_string, cigar_tuples = help_functions.parasail_alignment(read_seq, created_ref_seq)
    # # print('read', read_seq)
    # # print('ref',created_ref_seq)
    # print(read_aln)
    # print(ref_aln, tuple( mem.exon_part_id for mem in solution))
    # print(cigar_string)
    # print(read_aln == ref_aln)
    # # print([parts_to_exons[ mem.exon_part_id] for mem in solution])
    # print(chained_parts_ids)

    # # FSM
    # transcript = ''
    # if tuple(predicted_transcript) in parts_to_transcript_annotations[chr_id]:
    #     transcript = ",".join( tr for tr in parts_to_transcript_annotations[chr_id][tuple(predicted_transcript)])
    #     print()
    #     print('Found, FSM to:', transcript)
    #     print()
    #     return "FSM", transcript
    # elif tuple(predicted_splices) in splices_to_transcripts[chr_id]:
    #     transcript = ",".join( tr for tr in splices_to_transcripts[chr_id][tuple(predicted_splices)])  
    #     print()
    #     print('Found, FSM but not classified by parts to:', transcript)
    #     print()
    #     return "FSM", transcript

    # # else:

    # #     print('Did not find FSM', predicted_transcript)
    # #     for ann_tr in parts_to_transcript_annotations[chr_id]:
    # #         print(parts_to_transcript_annotations[chr_id][ann_tr] ,ann_tr)

    # # ISM
    # hits = [all_parts_pairs_annotations[chr_id][part_pair] for part_pair in all_parts_pairs_annotations[chr_id]]
    # # print(hits)
    # in_all_pairs = set.intersection(*hits)
    # for transcript_id in in_all_pairs:
    #     transcript_parts = transcripts_to_parts_annotations[chr_id][transcript_id]
    #     if contains(predicted_transcript, transcript_parts):
    #         # print("Found, ISM to", transcript_id )
    #         transcript = transcript_id
    #         return "ISM", transcript
    #     else:
    #         print(predicted_transcript, transcript)


    # # else:
    # all_sites_annotations_chr  = all_part_sites_annotations[chr_id] 
    # is_nic = True
    # for start, stop in predicted_transcript:
    #     if start not in all_sites_annotations_chr or stop not in all_sites_annotations_chr:
    #         is_nic = False
    # if is_nic:
    #     all_pairs_annotations_chr = all_parts_pairs_annotations[chr_id]
    #     is_nic_comb = True
    #     for start, stop in predicted_transcript:
    #         if (start, stop) not in all_pairs_annotations_chr:
    #             is_nic_comb = False


    #     if is_nic_comb:
    #         print()
    #         print('Found, NIC (new combination of exons):', tuple(predicted_transcript) )
    #         print()             
    #         # for ann_tr in parts_to_transcript_annotations[chr_id]:
    #         #     print(parts_to_transcript_annotations[chr_id][ann_tr] ,ann_tr)
    #         for ann_tr in splices_to_transcripts[chr_id]:
    #             print(splices_to_transcripts[chr_id][ann_tr] ,ann_tr)

    #         return  "NIC_comb", transcript

    #     else:
    #         print()
    #         print('Found, NIC (new donor-acceptor pair):', tuple(predicted_transcript) )
    #         print()   
    #         return   "NIC_novel", transcript          
    # # else:
    # #     print('Did not find NIC', predicted_transcript)

    # # else:
    # print()
    # print('NNC:', tuple(predicted_transcript) )
    # print()       
    # return "NNC", transcript
