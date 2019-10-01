
import sys
import re

# import parasail
import edlib

from collections import namedtuple
from modules import colinear_solver 
from modules import help_functions

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



def edlib_alignment(query, target):
    result = edlib.align(query, target, task="path", mode="HW")
    cigar_string = result["cigar"]
    locations = result['locations']
    query_alignment, target_alignment, cigar_tuples = cigar_to_seq(cigar_string, query, target)
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

    return locations, result['editDistance'], query_alignment, target_alignment



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

def main(solution, refs, parts_to_exons, exon_id_to_choordinates, read_seq, overlap_threshold, is_rc):


    # chained_parts_seq = []
    # chained_parts_ids = []
    prev_ref_stop = -1
    predicted_transcript = []
    predicted_exons = []
    covered_regions = []
    mam = namedtuple('Mam', ['x', 'y', 'c', 'd', 'val', "identity", "exon_id", "ref_chr_id"])
    mam_instance = []
    unique_exon_locations = set()
    for mem in solution:
        ref_chr_id, ref_start, ref_stop =  mem.exon_part_id.split('_')
        ref_start, ref_stop = int(ref_start), int(ref_stop)
        # get all exons associated with the part
        # print(parts_to_exons)
        unique_exon_locations.add((ref_chr_id, ref_start, ref_stop))
    print()
    print(unique_exon_locations)
    print()
    for (ref_chr_id, ref_start, ref_stop) in unique_exon_locations:
        exon_ids = parts_to_exons[ref_chr_id][(ref_start, ref_stop)]
        print(exon_ids)
        for exon_id in exon_ids:
            e_start, e_stop = exon_id_to_choordinates[exon_id]
            if e_stop - e_start > 10:
                # if is_rc:
                #     ref_seq = help_functions.reverse_complement(refs[ref_chr_id][e_start: e_stop])
                # else:
                ref_seq = refs[ref_chr_id][e_start: e_stop]
                # print(exon_id, e_stop - e_start)
                # align them to the read and get the best approxinate match

                locations, edit_distance, read_alignment, ref_alignment = edlib_alignment(ref_seq, read_seq)
                print(exon_id, e_start, e_stop, len(ref_seq), len(read_seq),edit_distance)

                if edit_distance >= 0:
                    score = (len(ref_seq) - edit_distance)/len(ref_seq)

                    if score > 0.7:
                        start, stop = locations[0]
                        covered_regions.append((start,stop, score, exon_id, ref_chr_id))
                        mam_tuple = mam(e_start, e_stop, start, stop, 
                                (stop - start + 1)*score, score,  exon_id, ref_chr_id)
                        mam_instance.append(mam_tuple)
            else:
                print("not aligning exons smaller than 10bp") # TODO: align these and take all locations

            if  e_stop - e_start >= 0.8*len(read_seq): # read is potentially contained within exon 
                print()
                print("aligning read to exon")
                locations, edit_distance, ref_alignment, read_alignment = edlib_alignment(read_seq, ref_seq)
                print(exon_id, e_start, e_stop, len(ref_seq), len(read_seq),edit_distance)
                print()
                if edit_distance >= 0:
                    score = (len(read_seq) - edit_distance)/len(read_seq)

                    if score > 0.7:
                        start, stop = 0, len(read_seq) - 1
                        covered_regions.append((start,stop, score, exon_id, ref_chr_id))
                        mam_tuple = mam(e_start, e_stop, start, stop, 
                                (stop - start + 1)*score, score,  exon_id, ref_chr_id)
                        mam_instance.append(mam_tuple)
                

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
    print("Solutiuon mam")
    print("val", value, "the mam mam_solution:", mam_solution)
    # print()
    # print()
    # print()
    return 'None', -1, value, mam_solution
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