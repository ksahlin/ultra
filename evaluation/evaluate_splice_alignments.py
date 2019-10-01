
from __future__ import print_function
import os,sys
import argparse
import re
import errno
import itertools

import pickle

from collections import defaultdict

import parasail
import pysam
import gffutils

'''
    Below code taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].replace(" ", "_"), [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break



def cigar_to_seq_mm2(read, full_r_seq, full_q_seq):
    # print()
    r_index = 0
    q_index = 0
    r_seq = full_r_seq[read.reference_start: read.reference_end + 1]
    q_seq = full_q_seq[read.query_alignment_start: read.query_alignment_end + 1]
    r_line, m_line, q_line = [], [], []

    cigar_tuples = read.cigartuples
    # print(cigar_tuples)
    # print(full_r_seq)
    # print(full_q_seq)

    # if query is softclipped in beginning or ref seq starts earlier
    if read.query_alignment_start > 0 or read.reference_start > 0:
        r_line.append( "-"*max(0, read.query_alignment_start - read.reference_start ) + full_r_seq[ : read.reference_start ] )
        q_line.append( "-"*max(0, read.reference_start - read.query_alignment_start ) + full_q_seq[ : read.query_alignment_start ]  ) 
        m_line.append("X"*max(read.query_alignment_start, read.reference_start ) )
        if cigar_tuples[0][0] in set([2,3,4]) and read.query_alignment_start == cigar_tuples[0][1]:
            del cigar_tuples[0]

    # if query is softclipped in end or ref seq extends beyond alignment 
    ref_end_offset = len(full_r_seq) - read.reference_end
    q_end_offset = len(full_q_seq) - read.query_alignment_end

    if q_end_offset > 0 or ref_end_offset > 0:
        # print(q_end_offset, ref_end_offset, "lol")
        if ref_end_offset:
            r_end = full_r_seq[ -ref_end_offset : ] + "-"*max(0, q_end_offset - ref_end_offset ) 
        else:
            r_end =  "-" * q_end_offset 

        if q_end_offset:
            q_end = full_q_seq[ -q_end_offset : ] + "-"*max(0, ref_end_offset - q_end_offset ) 
        else:
            q_end = "-" * ref_end_offset 

        m_end = "X"*max(q_end_offset, ref_end_offset )

        if cigar_tuples[-1][0] in set([2,3,4]) and q_end_offset == cigar_tuples[-1][1]:
            # print("HAHAHAHAHAHA")
            del cigar_tuples[-1]
        # del cigar_tuples[-1]
        # print(r_end, m_end, q_end)
    else:
        r_end, m_end, q_end = "", "", ""
    # print(cigar_tuples)

    for (op_type, op_len) in cigar_tuples:
        # op_len = int(op_len)
        if op_type == 0:
            ref_piece = r_seq[r_index: r_index + op_len]
            query_peace = q_seq[q_index: q_index + op_len]
            # print(ref_piece)
            # print(query_peace)

            r_line.append(ref_piece)
            q_line.append(query_peace)
            match_seq = ''.join(['|' if r_base.upper() == q_base.upper() else '*' for (r_base, q_base) in zip(ref_piece, query_peace)]) # faster with "".join([list of str]) instead of +=
                
            m_line.append(match_seq)
            r_index += op_len
            q_index += op_len

        elif op_type == 1:
            # insertion into reference
            r_line.append('-' * op_len)
            m_line.append(' ' * op_len)
            q_line.append(q_seq[q_index: q_index + op_len])
            #  only query index change
            q_index += op_len
        elif op_type == 2 or op_type == 3 or op_type == 4:
            # deletion from reference
            r_line.append(r_seq[r_index: r_index + op_len])
            m_line.append(' ' * op_len)
            q_line.append('-' * op_len)
            #  only ref index change
            r_index += op_len
            # print("here", r_index)

    r_line.append(r_end)
    m_line.append(m_end)
    q_line.append(q_end)    

    return "".join([s for s in r_line]), "".join([s for s in m_line]), "".join([s for s in q_line])

def cigar_to_seq_mm2_local(read, full_r_seq, full_q_seq):

    """
        Changed to parsing = and X instead of simple M.
    """
    # print()
    r_index = 0
    q_index = 0
    r_seq = full_r_seq[read.reference_start: read.reference_end + 1]
    q_seq = full_q_seq[read.query_alignment_start: read.query_alignment_end + 1]
    r_line, m_line, q_line = [], [], []

    cigar_tuples = read.cigartuples
    r_end, m_end, q_end = "", "", ""
    ins, del_, subs, matches = 0, 0, 0, 0
    # print(cigar_tuples)
    for (op_type, op_len) in cigar_tuples:
        # print((op_type, op_len))
        # op_len = int(op_len)
        if op_type == 7:
            ref_piece = r_seq[r_index: r_index + op_len]
            query_peace = q_seq[q_index: q_index + op_len]
            # print(ref_piece)
            # print(query_peace)

            r_line.append(ref_piece)
            q_line.append(query_peace)
            match_seq = ''.join(['|' if r_base.upper() == q_base.upper() else '*' for (r_base, q_base) in zip(ref_piece, query_peace)]) # faster with "".join([list of str]) instead of +=
            matches += op_len 
            m_line.append(match_seq)
            r_index += op_len
            q_index += op_len

        elif op_type == 8:
            ref_piece = r_seq[r_index: r_index + op_len]
            query_peace = q_seq[q_index: q_index + op_len]
            r_line.append(ref_piece)
            q_line.append(query_peace)

            subs += op_len 
            m_line.append('*' * op_len)
            r_index += op_len
            q_index += op_len


        elif op_type == 1:
            # insertion into reference
            r_line.append('-' * op_len)
            m_line.append(' ' * op_len)
            q_line.append(q_seq[q_index: q_index + op_len])
            #  only query index change
            q_index += op_len
            ins += op_len
        elif op_type == 2 or op_type == 3:
            # deletion from reference
            r_line.append(r_seq[r_index: r_index + op_len])
            m_line.append(' ' * op_len)
            q_line.append('-' * op_len)
            #  only ref index change
            r_index += op_len
            # print("here", r_index)
            if op_type == 2:
                del_ += op_len
        elif op_type == 4:
            # softclip
            pass
        elif op_type == 0:
            sys.exit("You have to use =/X cigar operators in sam file (M was detected). In minimap2 you can to provide --eqx flag to fix this.")


    r_line.append(r_end)
    m_line.append(m_end)
    q_line.append(q_end)    

    return "".join([s for s in r_line]), "".join([s for s in m_line]), "".join([s for s in q_line]), ins, del_, subs, matches


def decide_primary_locations(sam_file, args): # maybe this function is not needed if only one primary alignment from minimap2
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    reads_primary = {}
    reads_multiple_primary = set()
    reads_tmp = {}

    for read in SAM_file.fetch(until_eof=True):
        if read.flag == 0 or read.flag == 16:
            ins = sum([length for type_, length in read.cigartuples if type_ == 1])
            del_ = sum([length for type_, length in read.cigartuples if type_ == 2 and length < args.min_intron ])
            subs = sum([length for type_, length in read.cigartuples if type_ == 8])
            matches = sum([length for type_, length in read.cigartuples if type_ == 7])
            tot_align = ins + del_ + subs + matches
            identity = matches/float(tot_align)

            # has_large_del = [length for type_, length in read.cigartuples if type_ == 2 and length >=  args.min_intron ]
            # if has_large_del:
            #     print(has_large_del)
            if read.query_name in reads_primary:
                reads_multiple_primary.add(read.query_name)
                if identity >= reads_tmp[read.query_name][0] and  matches >= reads_tmp[read.query_name][1]:
                    reads_primary[read.query_name] = read
                    reads_tmp[read.query_name] = (identity, matches)
                elif identity <= reads_tmp[read.query_name][0] and  matches <= reads_tmp[read.query_name][1]:
                    continue
                else:
                    if identity * matches > reads_tmp[read.query_name][0] * reads_tmp[read.query_name][1]:
                        reads_primary[read.query_name] = read
                        reads_tmp[read.query_name] = (identity, matches)
                    else: 
                        continue

                        
                    # print( "Ambiguous, prefferred:", round(reads_tmp[read.query_name][0],2), reads_tmp[read.query_name][1], "over", round(identity,2), matches)

                #     if tot_align >= reads_tmp[read.query_name][1]:
                #         reads_primary[read.query_name] = read
                #     else:
                #         continue
                # elif tot_align >= reads_tmp[read.query_name][1]:
                #     if identity >= reads_tmp[read.query_name][0]:
                #         reads_primary[read.query_name] = read
                #     else:
                #         continue

            else:
                reads_primary[read.query_name] = read
                reads_tmp[read.query_name] = (identity, matches)
    print("TOTAL READS FLAGGED WITH MULTIPLE PRIMARY:", len(reads_multiple_primary))
    return reads_primary     


def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)



def get_deletion_sites(read):
    deletion_sites = []
    ref_pos = read.reference_start
    
    for i, (t,l) in enumerate(read.cigartuples):
        if t == 2:
            deletion_start = ref_pos
            ref_pos += l
            deletion_stop = ref_pos
            deletion_sites.append( (deletion_start, deletion_stop, l) )

        elif t == 7 or t== 0 or t == 8:
            ref_pos += l
        elif t == 3:
            # splice_start = ref_pos
            ref_pos += l
            # splice_stop = ref_pos
            # splice_sites.append( (splice_start, splice_stop) )

        elif t == 1 or t == 4 or t == 5: # insertion or softclip
            ref_pos += 0

        else: # reference skip or soft/hardclip "~", or match =
            print("UNEXPECTED!", t)
            sys.exit()

    return deletion_sites

# def get_deletion_coordiantes(read):
#     q_cigar = read.cigarstring
#     q_start = read.reference_start
#     q_end = read.reference_end
#     read_cigar_tuples = []
#     result = re.split(r'[=DXSMINH]+', q_cigar)
#     i = 0
#     for length in result[:-1]:
#         i += len(length)
#         type_ = q_cigar[i]
#         i += 1
#         read_cigar_tuples.append((int(length), type_ ))  

#     for type_, length in read.cigartuples:
    
#     read_deletion_coordinate_pairs = get_deletion_sites(read_cigar_tuples, q_start)
#     return read_deletion_coordinate_pairs


def get_error_rate_stats_per_read(reads_primary_locations, reads, annotated_splice_coordinates_pairs, args, reference = {}):
    # SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    alignments = {}
    alignments_detailed = {}
    read_index = 0
    all_called_intron_lengths = []
    for acc in reads_primary_locations:
        read = reads_primary_locations[acc]
        if read.flag == 0 or read.flag == 16:
            # print(read.is_reverse)
            read_index += 1

            ins = sum([length for type_, length in read.cigartuples if type_ == 1])
            del_old = sum([length for type_, length in read.cigartuples if type_ == 2 and length <  args.min_intron])
            subs = sum([length for type_, length in read.cigartuples if type_ == 8])
            matches = sum([length for type_, length in read.cigartuples if type_ == 7])
            called_intron_lengths = [length for type_, length in read.cigartuples if type_ == 3]


            # Deletion is treated separately because of potential short introns -- if in annotation
            if read.reference_name in annotated_splice_coordinates_pairs:
                del_new = 0
                chr_coordinate_pairs = annotated_splice_coordinates_pairs[read.reference_name]
                read_deletion_coordinate_pairs = get_deletion_sites(read)
                for start_del, stop_del, del_length in read_deletion_coordinate_pairs:
                    if (start_del, stop_del) not in chr_coordinate_pairs:
                        del_new += del_length
                    else:
                        print("True intron masked as deletion, length:", del_length, read.reference_name, start_del, stop_del)

            else:
                del_new = sum([length for type_, length in read.cigartuples if type_ == 2 and length <  args.min_intron])

            all_called_intron_lengths.append(called_intron_lengths)

            assert read.query_name not in alignments
            alignments[read.query_name] = (ins, del_new, subs, matches, read.reference_name, read.reference_start, read.reference_end + 1, read.flag, read_index)

        else:
            pass
            # print("secondary", read.flag, read.reference_name) 
    # SAM_file.close()

    # multiple_alingment_counter = 0
    # for read in alignments:
    #     if len(alignments[read]) > 1:
    #         # print(alignments[read])
    #         multiple_alingment_counter += 1
    # print("Nr reads with multiple primary alignments:", multiple_alingment_counter)
    print( '100 smallest introns as N from the aligner:', sorted([l for lst in all_called_intron_lengths for l in lst])[:100] )
    # sys.exit()
    return alignments, alignments_detailed

def get_summary_stats(reads, quantile):
    tot_ins, tot_del, tot_subs, tot_match = 0, 0, 0, 0
    sorted_reads = sorted(reads.items(), key = lambda x: sum(x[1][0:3])/float(sum(x[1][0:4])) )
    for acc, (ins, del_, subs, matches, chr_id, reference_start, reference_end, sam_flag, read_index) in sorted_reads[ : int(len(sorted_reads)*quantile)]:
        # (ins, del_, subs, matches) = reads[acc]

        tot_ins += ins
        tot_del += del_
        tot_subs += subs
        tot_match += matches

    sum_aln_bases = tot_ins + tot_del + tot_subs + tot_match

    return tot_ins, tot_del, tot_subs, tot_match, sum_aln_bases


def print_detailed_values_to_file(alignments_dict, annotations_dict, reads, outfile, reads_unaligned_in_other_method, read_type):
    # read_calss is FSM, NIC, NNC, ISM
    # donwnload human gtf file to compare against, check how sqanti does it.
    # also sent isONclust tsv file to this script to get cluster size 
    alignments_sorted = sorted(alignments_dict.items(), key = lambda x: sum(x[1][0:3])/float(sum(x[1][0:4])) )

    for (acc, (ins, del_, subs, matches, chr_id, reference_start, reference_end, sam_flag, read_index)) in alignments_sorted:
        error_rate = round( 100* (ins + del_ + subs) /float( (ins + del_ + subs + matches) ), 4 )
        read_class = annotations_dict[acc] #"NA" # annotations_dict[acc]
        read_length = len(reads[acc])
        is_unaligned_in_other_method = 1 if acc in reads_unaligned_in_other_method else 0
        # is_missing_from_clustering_or_correction = 1 if acc in reads_missing_from_clustering_correction_output else 0

        info_tuple = (acc, read_type, ins, del_, subs, matches, error_rate, read_length, is_unaligned_in_other_method, *read_class, chr_id, reference_start, reference_end + 1, sam_flag) # 'tot_splices', 'read_sm_junctions', 'read_nic_junctions', 'fsm', 'nic', 'ism', 'nnc', 'no_splices'  )
        # print(*info_tuple)
        # outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}\n".format(*info_tuple))
        outfile.write( ",".join( [str(item) for item in info_tuple] ) + "\n")


def print_quantile_values(alignments_dict):

    alignments_sorted = sorted(alignments_dict.items(), key = lambda x: sum(x[1][0:3])/float(sum(x[1][0:4])) )
    print(alignments_sorted[-5:])

    sorted_error_rates = [ sum(tup[0:3])/float(sum(tup[0:4])) for acc, tup in alignments_sorted]

    insertions = [ tup[0]/float(sum(tup[0:4])) for acc, tup in alignments_sorted]
    deletions = [ tup[1]/float(sum(tup[0:4])) for acc, tup in alignments_sorted]
    substitutions = [ tup[2]/float(sum(tup[0:4])) for acc, tup in alignments_sorted]

    n = len(sorted_error_rates)
    quantiles = [0, 1.0/20, 1.0/10, 1.0/4, 1.0/2, 3.0/4, 9.0/10, 19.0/20, 1 ]
    quantile_errors = []
    quantile_insertions = []
    quantile_deletions = []
    quantile_substitutions = []

    for q in quantiles:
        if q == 1:
            quantile_errors.append(sorted_error_rates[int(q*n)-1])
            quantile_insertions.append(insertions[int(q*n)-1])
            quantile_deletions.append(deletions[int(q*n)-1])
            quantile_substitutions.append(substitutions[int(q*n)-1])
        else:
            quantile_errors.append(sorted_error_rates[int(q*n)])
            quantile_insertions.append(insertions[int(q*n)])
            quantile_deletions.append(deletions[int(q*n)])
            quantile_substitutions.append(substitutions[int(q*n)])

    # quantile_errors = [sorted_error_rates[0], sorted_error_rates[int(n/20)], sorted_error_rates[int(n/10)], 
    #             sorted_error_rates[int(n/4)], sorted_error_rates[int(n/2)], 
    #             sorted_error_rates[int(3*n/4)], sorted_error_rates[int(9*n/10)], sorted_error_rates[int(19*n/20)], 
    #             sorted_error_rates[-1]]

    return quantile_errors, quantile_insertions, quantile_deletions, quantile_substitutions


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

    return  "".join([s for s in q_aln]), "".join([s for s in r_aln])


def parasail_alignment(read, reference, x_acc = "", y_acc = "", match_score = 2, mismatch_penalty = -2, opening_penalty = 2, gap_ext = 1, ends_discrepancy_threshold = 0):
    user_matrix = parasail.matrix_create("ACGT", match_score, mismatch_penalty)
    result = parasail.sg_trace_scan_16(read, reference, opening_penalty, gap_ext, user_matrix)
    if result.saturated:
        print("SATURATED!")
        result = parasail.sg_trace_scan_32(read, reference, opening_penalty, gap_ext, user_matrix)
    if sys.version_info[0] < 3:
        cigar_string = str(result.cigar.decode).decode('utf-8')
    else:
        cigar_string = str(result.cigar.decode, 'utf-8')
    
    read_alignment, ref_alignment = cigar_to_seq(cigar_string, read, reference)
    return read_alignment, ref_alignment

# def get_cluster_sizes(cluster_file, reads):
#     cluster_sizes = defaultdict(int)
#     tmp_reads = {}
#     for line in open(args.cluster_file, "r"):
#         cl_id, acc = line.split() 
#         if acc in reads:
#             cluster_sizes[cl_id] += 1  
#             tmp_reads[acc] = cl_id
#         else:
#             cluster_sizes[cl_id] += 1  
#             tmp_reads[acc.split("_")[0]] = cl_id

#     reads_to_cluster_size = {}
#     for acc, cl_id in tmp_reads.items():
#         reads_to_cluster_size[acc] = cluster_sizes[cl_id]

#     return reads_to_cluster_size


def get_splice_sites(cigar_tuples, first_exon_start, minimum_annotated_intron, annotated_chr_coordinate_pairs):
    splice_sites = []
    ref_pos = first_exon_start
    
    for i, (l,t) in enumerate(cigar_tuples):
        if t == "D":
            # if l >= minimum_annotated_intron -1:
            #     # print("long del", l)
            #     splice_start = ref_pos
            #     ref_pos += l
            #     splice_stop = ref_pos
            #     splice_sites.append( (splice_start, splice_stop) )
                
            # else:
                if (ref_pos, ref_pos + l) in annotated_chr_coordinate_pairs:
                    splice_sites.append( (ref_pos, ref_pos + l) )
                    print("HEERE")
                ref_pos += l

                # if l > 15:
                #     print("Large deletion!!!", l)


        elif t == "=" or t== "M" or t == "X":
            ref_pos += l
        elif t == "N":
            splice_start = ref_pos
            ref_pos += l
            splice_stop = ref_pos
            splice_sites.append( (splice_start, splice_stop) )

        elif t == "I" or t == "S" or t == "H": # insertion or softclip
            ref_pos += 0

        else: # reference skip or soft/hardclip "~", or match =
            print("UNEXPECTED!", t)
            sys.exit()

    return splice_sites

def get_read_candidate_splice_sites(reads_primary_locations, minimum_annotated_intron, annotated_splice_coordinates_pairs):
    read_splice_sites = {}
    for acc in reads_primary_locations:
        read = reads_primary_locations[acc]
        if read.flag == 0 or read.flag == 16:
            # compare cs tag at intron sites
            q_cigar = read.cigarstring
            q_start = read.reference_start
            q_end = read.reference_end
            read_cigar_tuples = []
            result = re.split(r'[=DXSMINH]+', q_cigar)
            i = 0
            for length in result[:-1]:
                i += len(length)
                type_ = q_cigar[i]
                i += 1
                read_cigar_tuples.append((int(length), type_ ))  

            if read.reference_name in annotated_splice_coordinates_pairs:
                annotated_chr_coordinate_pairs = annotated_splice_coordinates_pairs[read.reference_name]
            else:
                annotated_chr_coordinate_pairs = set()

            read_splice_sites[read.query_name] = {}  
            read_splice_sites[read.query_name][read.reference_name] = get_splice_sites(read_cigar_tuples, q_start, minimum_annotated_intron, annotated_chr_coordinate_pairs)
            # if acc == '261:647|8b7c0adb-b728-4191-9257-97bdc4e73ba0_runid=8c239806e6f576cd17d6b7d532976b1fe830f9c6_sampleid=pcs109_sirv_mix2_LC_read=117640_ch=474_start_time=2019-04-13T18:27:44Z_strand=-_strand=+':
            #     print(get_splice_sites(read_cigar_tuples, q_start, minimum_annotated_intron, annotated_chr_coordinate_pairs))
            #     sys.exit()
            # read_deletion_coordinate_pairs = get_deletion_sites(read)
            # read_splice_sites[read.query_name][read.reference_name] = read_deletion_coordinate_pairs
        

    return read_splice_sites


def get_annotated_splicesites(ref_gff_file, infer_genes, outfolder):
    db_name = os.path.join(outfolder, 'database.db')
    if infer_genes:
        fn = gffutils.example_filename(ref_gff_file)
        db = gffutils.create_db(fn, dbfn=db_name, force=True, keep_order=True, merge_strategy='merge', 
                                sort_attribute_values=True)
        db = gffutils.FeatureDB(db_name, keep_order=True)
    else:
        fn = gffutils.example_filename(ref_gff_file)
        db = gffutils.create_db(fn, dbfn=db_name, force=True, keep_order=True, merge_strategy='merge', 
                                sort_attribute_values=True, disable_infer_genes=True, disable_infer_transcripts=True)
        db = gffutils.FeatureDB(db_name, keep_order=True)

    
    # for tr in db.children(gene, featuretype='transcript', order_by='start'):
    #     # print(tr.id, dir(tr)) 
    #     splice_sites = []
    #     for j in db.children(tr, featuretype='exon', order_by='start'):
    #         # print(j, j.start, j.end)
    #         splice_sites.append(j.start -1)
    #         splice_sites.append(j.end)
    #     splice_sites_tmp = splice_sites[1:-1]
    #     splice_sites = []
    #     for i in range(0, len(splice_sites_tmp),2):
    #         splice_sites.append( (splice_sites_tmp[i], splice_sites_tmp[i+1]) )
    #     # splice_sites = [item for item in zip(splice_sites[:-1], splice_sites[1:])]
    #     ref_isoforms[tr.id] = splice_sites


    splice_coordinates = {} # to calc individual fraction of correct sites and NIC
    splice_coordinates_pairs = {} 
    ref_isoforms = {} # To calculate Full splice matches

    # gene_graphs = {} # gene_id : { (exon_start, exon_stop) : set() }
    # collapsed_exon_to_transcript = {}
    minimum_annotated_intron = 1000000000
    for gene in db.features_of_type('gene'):
        chromosome = str(gene.seqid)
        if chromosome not in ref_isoforms:
            ref_isoforms[chromosome] = {}
        if chromosome not in splice_coordinates:
            splice_coordinates[chromosome] = set()
            splice_coordinates_pairs[chromosome] = set()


        
        # #add splice sites
        # splice_sites_tmp = []
        # for exon in db.children(gene, featuretype='exon', order_by='start'):
        #     exon_start, exon_end = exon.start - 1, exon.stop # double check if gff is 1-indexed here!
        #     splice_sites_tmp.append(exon_start)
        #     splice_sites_tmp.append(exon_end)

        # splice_sites = []
        # for i in range(0, len(splice_sites_tmp[1:-1]),2):
        #     splice_start, splice_stop = splice_sites_tmp[i], splice_sites_tmp[i+1]
        #     splice_sites.append( (splice_start, splice_stop) )
        #     splice_coordinates[chromosome].add(splice_start)
        #     splice_coordinates[chromosome].add(splice_stop)
        #     splice_coordinates_pairs[chromosome].add( (splice_start, splice_stop) )
            
            # collapsed_exon_to_transcript[gene.id][ (exon.start, exon.stop) ].update([ transcript_tmp for transcript_tmp in  exon.attributes['transcript_id']])
            # if (exon.start, exon.stop) in already_parsed_exons: 
            # gene_graph.add_node( (exon.start, exon.stop), weight=1  )
            # print(gene_graph.nodes[(exon.start, exon.stop)])

        #add annotated transcripts
        for transcript in db.children(gene, featuretype='transcript', order_by='start'):
            # print(dir(transcript))
            consecutive_exons = [exon for exon in db.children(transcript, featuretype='exon', order_by='start')]
            # print('transcript', transcript.id, transcript.start, transcript.stop, [ (exon.start, exon.stop) for exon in db.children(transcript, featuretype='exon', order_by='start')])
            tmp_splice_sites = []
            for e1,e2 in zip(consecutive_exons[:-1], consecutive_exons[1:]):
                tmp_splice_sites.append( (e1.stop, e2.start -1 ))           
                splice_coordinates[chromosome].add(e1.stop)
                splice_coordinates[chromosome].add(e2.start -1 )
                splice_coordinates_pairs[chromosome].add( (e1.stop, e2.start -1 ) )

                if e2.start -1 - e1.stop < minimum_annotated_intron:
                    minimum_annotated_intron = e2.start -1 - e1.stop
                # print('exon', exon.id, exon.start, exon.stop)
            
            ref_isoforms[chromosome][tuple(tmp_splice_sites)] = transcript.id

      
    return ref_isoforms, splice_coordinates, splice_coordinates_pairs, minimum_annotated_intron

from collections import namedtuple
def get_splice_classifications(annotated_ref_isoforms, annotated_splice_coordinates, annotated_splice_coordinates_pairs,  all_reads_splice_sites, ref_seqs, reads_primary_locations):

    total_reads = 0
    
    #individual splice sites
    total_individual_true = 0
    total_individual_in_data = 0

    #pairs of splice sites
    total_pairs_in_data = 0
    tot_nic_pairs = 0
    tot_sm_pairs = 0

    #whole structure
    total_transcript_fsm = 0
    total_transcript_nic = 0
    total_transcript_ism = 0
    total_transcript_nnc = 0
    total_transcript_no_splices = 0

    read_annotations = {}
    alignments_not_matching_annotated_sites = 0
    # print(annotated_ref_isoforms["SIRV5"])
    # for tr in annotated_ref_isoforms["SIRV5"]:
    #     print(tr)
    canonical_splice = 0
    all_splice = 0
    for read_acc in all_reads_splice_sites:
        read_annotations[read_acc] = {}
        total_reads += 1
        assert len(all_reads_splice_sites[read_acc]) == 1
        chr_id, read_splice_sites = list(all_reads_splice_sites[read_acc].items())[0]

        # print(reads_primary_locations[read_acc].flag, read_splice_sites)
        # print(chr_id)
        # print(annotated_splice_coordinates)
        if chr_id not in annotated_splice_coordinates:
            alignments_not_matching_annotated_sites += 1
            annotated_sites = set()
            annotated_pairs = set()
            annotated_isoforms = set()
        else:
            annotated_sites = annotated_splice_coordinates[chr_id]
            annotated_pairs = annotated_splice_coordinates_pairs[chr_id]
            annotated_isoforms = annotated_ref_isoforms[chr_id]

        # print(annotated_pairs)
        # print( chr_id, read_splice_sites)
        # print(annotated_ref_isoforms[chr_id])
        read_sm_junctions = 0
        read_nic_junctions = 0
        read_splice_letters = []
        read_splice_choords = []
        total_individual_in_data += len(read_splice_sites)*2
        total_pairs_in_data += len(read_splice_sites)

        for splice_site in read_splice_sites:
            start_sp, stop_sp = splice_site
            if reads_primary_locations[read_acc].flag == 0:
                donor = ref_seqs[chr_id][start_sp: start_sp + 2] 
                acceptor = ref_seqs[chr_id][stop_sp - 2: stop_sp]
            else:
                acceptor = reverse_complement(ref_seqs[chr_id][start_sp: start_sp + 2])
                donor = reverse_complement(ref_seqs[chr_id][stop_sp - 2: stop_sp])

            if donor == "GT" and acceptor == "AG":
                canonical_splice += 1
            all_splice += 1

            read_splice_letters.append( donor + str("-") + acceptor )
            read_splice_choords.append( str(start_sp) + str("-") + str(stop_sp) )
            # print(splice_site)
            if (start_sp, stop_sp) in annotated_pairs:
                tot_sm_pairs += 1 
                total_individual_true += 2
                read_sm_junctions += 1

            elif start_sp in annotated_sites and stop_sp in annotated_sites:
                # print((start_sp, stop_sp), annotated_pairs )
                tot_nic_pairs += 1
                total_individual_true += 2
                read_nic_junctions += 1

            elif start_sp in annotated_sites:
                total_individual_true += 1  

            elif stop_sp in annotated_sites:
                total_individual_true += 1  

            # for sp in splice_site:
            #     # print(sp)
            #     # print(annotated_sites)
            #     total += 1
            #     if sp in annotated_sites:
            #         total_true += 1

        # check set intersection between read splice sites and annotated splice sites
        read_nic, read_ism, read_nnc, read_no_splices = 0,0,0,0
        read_type = ''
        transcript_fsm_id = "NA"
        # print(tuple(read_splice_sites))
        # print(annotated_ref_isoforms[chr_id])
        # print(annotated_ref_isoforms)

        if len(read_splice_sites) > 0:

            if tuple(read_splice_sites) in annotated_isoforms:
                total_transcript_fsm += 1
                read_type = 'FSM'
                transcript_fsm_id = annotated_isoforms[ tuple(read_splice_sites) ]


            elif len(read_splice_sites) == read_sm_junctions + read_nic_junctions:
                if read_nic_junctions >= 1:
                    total_transcript_nic += 1
                    read_type = 'NIC'
                    print('NIC', read_acc)

                else:
                    total_transcript_ism += 1
                    read_type = 'ISM'

                    # print(read_acc)
                    # print(tuple(read_splice_sites))
                    # for an in annotated_ref_isoforms[chr_id]:
                    #     print(an)
                    # print()
                    # print(annotated_ref_isoforms[chr_id])
            else:
                total_transcript_nnc += 1
                read_nnc = 1
                read_type = 'NNC'
        else:
            total_transcript_no_splices += 1                
            read_type = 'NO_SPLICE'

        read_annotation = namedtuple('Annotation', ['tot_splices', 'read_sm_junctions', 'read_nic_junctions', 'annotation', "donor_acceptors", "donor_acceptors_choords", "transcript_fsm_id" ])
        if read_splice_letters:
            donor_acceptors = ":".join([str(item) for item in read_splice_letters])
            donor_acceptors_choords = ":".join([str(item) for item in read_splice_choords])

        else: 
            donor_acceptors = "NA"
            donor_acceptors_choords = "NA"
        read_annotations[read_acc] = read_annotation( len(read_splice_sites), read_sm_junctions, read_nic_junctions, read_type, donor_acceptors, donor_acceptors_choords, transcript_fsm_id )
                # print("FSM!!")
    # print(annotated_ref_isoforms[chr_id])
    # print( tuple(read_splice_sites))
    print("Total splice sizes found in cigar in reads (individual, pairs):", total_individual_in_data, total_pairs_in_data, "total matching annotations (individual):", total_individual_true,
             "total annotated junctions (splice match pairs):", tot_sm_pairs,  "total NIC junctions in (pairs):", tot_nic_pairs, "total reads aligned:", total_reads)
    print("total transcripts FSM:", total_transcript_fsm)
    print("total transcripts NIC:", total_transcript_nic)
    print("total transcripts ISM:", total_transcript_ism)
    print("total transcripts NNC:", total_transcript_nnc)
    print("total transcripts no splice sites:", total_transcript_no_splices)
    print("total splice sites:", all_splice)
    print("GT-AG splice sites:", canonical_splice)

    return read_annotations

def pickle_dump(data, filename):
    with open(os.path.join(args.outfolder,filename), 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)


def pickle_load(filename):
    with open(filename, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        data = pickle.load(f)
    return data

# def print_exact_alignments(reads, mm2_reads):
#     tot = 0
#     for acc in reads:
#         if acc not in mm2_reads:
#             print(acc, "Not in corrected reads")
#         else:
#             orig_seq = reads[acc]
#             mm2_seq = corr_reads[acc]
#             read_alignment, ref_alignment = parasail_alignment(orig_seq, corr_seq)
#             orig_dlen = [d_len for ch, d_len in  [(ch, len(list(_))) for ch, _ in itertools.groupby(read_alignment)] if ch == "-" ]
#             corr_dlen = [d_len for ch, d_len in  [(ch, len(list(_))) for ch, _ in itertools.groupby(ref_alignment)] if ch == "-" ] 
#             if corr_dlen:
#                 max_orig_dlen = max( orig_dlen )
#             else:
#                 max_orig_dlen = 0
#             if corr_dlen:
#                 max_corr_dlen = max( corr_dlen )
#             else:
#                 max_corr_dlen = 0

#             if max_corr_dlen > 10 or max_orig_dlen > 10:
#                 print("orig", acc, read_alignment)
#                 print(orig[acc])
#                 print("corr", acc, ref_alignment)
#                 print(corr[acc])
#                 print(corr_detailed[acc][0])
#                 print(corr_detailed[acc][1])
#                 print(corr_detailed[acc][2])

#                 tot += 1
#                 print()
#     print("TOT structural diffs:", tot)

def main(args):
    if args.load_database:
        print()
        print("LOADING FROM PRECOMPUTED DATABASE")
        print()
        annotated_ref_isoforms = pickle_load(os.path.join( args.outfolder, 'annotated_ref_isoforms.pickle') )
        annotated_splice_coordinates = pickle_load(os.path.join( args.outfolder, 'annotated_splice_coordinates.pickle') )
        annotated_splice_coordinates_pairs = pickle_load(os.path.join( args.outfolder, 'annotated_splice_coordinates_pairs.pickle') )
        minimum_annotated_intron = pickle_load(os.path.join( args.outfolder, 'minimum_annotated_intron.pickle') )
    else:
        annotated_ref_isoforms, annotated_splice_coordinates, annotated_splice_coordinates_pairs, minimum_annotated_intron = get_annotated_splicesites(args.gff_file, args.infer_genes, args.outfolder)
        pickle_dump(annotated_ref_isoforms, os.path.join( args.outfolder, 'annotated_ref_isoforms.pickle') )
        pickle_dump(annotated_splice_coordinates, os.path.join( args.outfolder, 'annotated_splice_coordinates.pickle') )
        pickle_dump(annotated_splice_coordinates_pairs, os.path.join( args.outfolder, 'annotated_splice_coordinates_pairs.pickle') )
        pickle_dump(minimum_annotated_intron, os.path.join( args.outfolder, 'minimum_annotated_intron.pickle') )

    reads = { acc.split()[0] : seq for i, (acc, (seq, qual)) in enumerate(readfq(open(args.reads, 'r')))}
    torkel_primary_locations = decide_primary_locations(args.torkel_sam, args)
    mm2_primary_locations = decide_primary_locations(args.mm2_sam, args)

    torkel, torkel_detailed = get_error_rate_stats_per_read(torkel_primary_locations, reads, annotated_splice_coordinates_pairs, args)
    mm2, mm2_detailed = get_error_rate_stats_per_read(mm2_primary_locations, reads, annotated_splice_coordinates_pairs, args)

    print( "Reads successfully aligned:", len(mm2),len(torkel))

    # if args.align:
    #     print_exact_alignments(reads, corr_reads)


    quantile_tot_torkel, quantile_insertions_torkel, quantile_deletions_torkel, quantile_substitutions_torkel = print_quantile_values(torkel)
    quantile_tot_mm2, quantile_insertions_mm2, quantile_deletions_mm2, quantile_substitutions_mm2 = print_quantile_values(mm2)

    torkel_stats = get_summary_stats(torkel, 1.0)
    mm2_stats = get_summary_stats(mm2, 1.0)
    
    print("Distribution of error rates (Percent)")
    print("Reads, Best, top 5%, top 10%, top 25%, Median, top 75%, top 90%, top 95%, Worst")
    print("torkel,{0},{1},{2},{3},{4},{5},{6},{7},{8}".format( *[round(100*round(x,3), 2) for x in quantile_tot_torkel ] ))
    print("mm2,{0},{1},{2},{3},{4},{5},{6},{7},{8}".format( *[round(100*round(x,3), 2) for x in quantile_tot_mm2 ] ))

    outfile = open(os.path.join(args.outfolder, "results.csv"), "w")
    outfile.write("torkel,tot,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n".format( *[round(100*round(x,3), 2) for x in quantile_tot_torkel ], *torkel_stats, len(torkel)))
    outfile.write("mm2,tot,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n".format( *[round(100*round(x,3), 2) for x in quantile_tot_mm2 ], *mm2_stats, len(mm2)))
    
    outfile.write("torkel,ins,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n".format( *[round(100*round(x,3), 2) for x in quantile_insertions_torkel ], *torkel_stats, len(torkel)))
    outfile.write("mm2,ins,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n".format( *[round(100*round(x,3), 2) for x in quantile_insertions_mm2 ], *mm2_stats, len(mm2)))

    outfile.write("torkel,del,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n".format( *[round(100*round(x,3), 2) for x in quantile_deletions_torkel ], *torkel_stats, len(torkel)))
    outfile.write("mm2,del,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n".format( *[round(100*round(x,3), 2) for x in quantile_deletions_mm2 ], *mm2_stats, len(mm2)))

    outfile.write("torkel,subs,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n".format( *[round(100*round(x,3), 2) for x in quantile_substitutions_torkel ], *torkel_stats, len(torkel)))
    outfile.write("mm2,subs,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n".format( *[round(100*round(x,3), 2) for x in quantile_substitutions_mm2 ], *mm2_stats, len(mm2)))

    outfile.close()    

    ############# Splice site analysis ########################################
    ###########################################################################

    # print(annotated_ref_isoforms)
    # print(annotated_splice_coordinates)
    print("SHORTEST INTRON:", minimum_annotated_intron)
    minimum_annotated_intron = max(minimum_annotated_intron,  args.min_intron)
    mm2_splice_sites = get_read_candidate_splice_sites(mm2_primary_locations, minimum_annotated_intron, annotated_splice_coordinates_pairs)
    torkel_splice_sites = get_read_candidate_splice_sites(torkel_primary_locations, minimum_annotated_intron, annotated_splice_coordinates_pairs)

    refs = { acc.split()[0] : seq for i, (acc, (seq, _)) in enumerate(readfq(open(args.refs, 'r')))}

    print('MINIMAP2')
    mm2_splice_results = get_splice_classifications(annotated_ref_isoforms, annotated_splice_coordinates, annotated_splice_coordinates_pairs, mm2_splice_sites, refs, mm2_primary_locations)
    print('uLTRA')
    torkel_splice_results = get_splice_classifications(annotated_ref_isoforms, annotated_splice_coordinates, annotated_splice_coordinates_pairs, torkel_splice_sites, refs, torkel_primary_locations)
    # reads_to_cluster_size = get_cluster_sizes(args.cluster_file, reads)

    reads_unaligned_in_torkel = set(reads.keys()) - set(torkel_primary_locations.keys())
    reads_unaligned_in_mm2 = set(reads.keys()) - set(mm2_primary_locations.keys()) 

    detailed_results_outfile = open(os.path.join(args.outfolder, "results_per_read.csv"), "w")
    detailed_results_outfile.write("acc,read_type,ins,del,subs,matches,error_rate,read_length,cluster_size, is_unaligned_in_other_method,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag\n")

    print_detailed_values_to_file(mm2, mm2_splice_results, reads, detailed_results_outfile, reads_unaligned_in_torkel, "minimap2")    

    print_detailed_values_to_file(torkel, torkel_splice_results, reads, detailed_results_outfile, reads_unaligned_in_mm2, "uLTRA")
    detailed_results_outfile.close()

    print()
    print("Reads successfully aligned (original/corrected):", len(torkel),len(mm2))
    print("Total reads", len(reads))
    print("READS UNALIGNED (ORIGINAL/CORRECTED):", len(reads_unaligned_in_torkel), len(reads_unaligned_in_mm2) )

    ###########################################################################
    ###########################################################################




    # print(torkel_sorted)
    # print(",".join([s for s in torkel_stats]))
    # print(",".join([s for s in mm2_stats]))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('torkel_sam', type=str, help='Path to the original read file')
    parser.add_argument('mm2_sam', type=str, help='Path to the corrected read file')
    parser.add_argument('reads', type=str, help='Path to the read file')
    parser.add_argument('refs', type=str, help='Path to the refs file')
    parser.add_argument('gff_file', type=str, help='Path to the refs file')
    parser.add_argument('outfolder', type=str, help='Output path of results')
    parser.add_argument('--min_intron', type=int, default=15, help='Threchold for what is counted as varation/intron in alignment as opposed to deletion.')
    parser.add_argument('--infer_genes', action= "store_true", help='Include pairwise alignment of original and corrected read.')
    parser.add_argument('--load_database', action= "store_true", help='Load already computed splice junctions and transcript annotations instead of constructing a new database.')
    # parser.add_argument('--align', action= "store_true", help='Include pairwise alignment of original and corrected read.')

    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)

