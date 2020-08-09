
from __future__ import print_function
import os,sys
import argparse
import re
import errno
import itertools

import pickle

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")


import numpy as np
import seaborn as sns
import pandas as pd

from matplotlib_venn import venn3, venn3_circles, venn2

from collections import defaultdict

# import parasail
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
        name, seqs, last = last[1:].split()[0], [], None
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



def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)


def is_overlapping(a_start,a_stop, b_start,b_stop):
    return (int(a_start) <= int(b_start) <= int(a_stop) )  or (int(a_start) <= int(b_stop) <= int(a_stop)) or (int(b_start) <= int(a_start) <= int(a_stop) <= int(b_stop) )


def parse_differing_location_reads(csv_file):
    reads_isonalign = {}
    reads_minimap2 = {}
    reads_desalt = {}
    for line in open(csv_file,'r'):
        #print(line)
        acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_exonic = line.strip().split(",")
        if algorithm == 'uLTRA':
            reads_isonalign[acc] =  (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_exonic) 
        if algorithm == 'minimap2_GTF':
            reads_minimap2[acc] = (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_exonic) 
        if algorithm == 'deSALT_GTF':
            reads_desalt[acc] = (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_exonic)
    return reads_isonalign, reads_minimap2, reads_desalt


def get_FSM_concordance(reads_ultra, reads_minimap2, reads_desalt):
    ds_fsm_distribution = defaultdict(set)
    desalt = set()
    for acc in reads_desalt:
        ds_annot, transcript_fsm_id = reads_desalt[acc][7], reads_desalt[acc][10]
        if ds_annot == 'FSM':
            ds_fsm_distribution[transcript_fsm_id].add(acc)
            desalt.add( (acc, transcript_fsm_id) ) 

    mm2_fsm_distribution = defaultdict(set)
    minimap2 = set()
    for acc in reads_minimap2:
        mm2_annot, transcript_fsm_id = reads_minimap2[acc][7], reads_minimap2[acc][10]
        if mm2_annot == 'FSM':
            mm2_fsm_distribution[transcript_fsm_id].add(acc)
            minimap2.add( (acc, transcript_fsm_id) ) 

    ultra_fsm_distribution = defaultdict(set)
    ultra = set()
    for acc in reads_ultra:
        ultra_annot, transcript_fsm_id = reads_ultra[acc][7], reads_ultra[acc][10]
        if ultra_annot == 'FSM':
            ultra_fsm_distribution[transcript_fsm_id].add(acc)
            ultra.add( (acc, transcript_fsm_id) ) 


    # print("ds_fsm_distribution:", sorted(ds_fsm_distribution.items(), key = lambda x: x[1]))        
    # print("mm2_fsm_distribution:", sorted(mm2_fsm_distribution.items(), key = lambda x: x[1]))        
    # print("ultra_fsm_distribution:", sorted(ultra_fsm_distribution.items(), key = lambda x: x[1]))        
    print("Desalt nr unique isoforms mapped to:", len(ds_fsm_distribution))
    print("minimap2 nr unique isoforms mapped to:", len(mm2_fsm_distribution))
    print("Ultra nr unique isoforms mapped to:", len(ultra_fsm_distribution))

    print(len(ultra), len(desalt), len(minimap2))

    # a = ultra
    # b = desalt
    # c = minimap2
    a_not_b_c = ultra - (desalt | minimap2)
    b_not_a_c = desalt - (ultra | minimap2)
    c_not_a_b = minimap2 - (ultra | desalt)
    a_b_not_c = (ultra & desalt) - minimap2
    a_c_not_b = (ultra & minimap2) - desalt
    b_c_not_a = (desalt & minimap2) - ultra
    a_b_c = ultra & desalt & minimap2

    print("BAD:")
    print("desalt and minimap2:", len(b_c_not_a))
    print("desalt unique:", len(b_not_a_c))
    print("minimap2 unique:", len(c_not_a_b))
    print()
    print("GOOD:")
    print("Ultra and desalt:", len(a_b_not_c))
    print("Ultra and minimap2:", len(a_c_not_b))
    print()
    
    print("NEUTRAL")
    print("ultra unique:", len(a_not_b_c))
    print("In all", len(a_b_c))

    print("ALL FSM READS:", len( (ultra | desalt | minimap2 )) )
    return [ultra, desalt, minimap2], [ultra_fsm_distribution, ds_fsm_distribution, mm2_fsm_distribution]


def venn(data_for_venn, outfolder, filename):
    ultra, desalt, minimap2 = data_for_venn
    if type(ultra) == set:
        total = len((ultra | desalt | minimap2 ))
    else:
        ultra = set(ultra.keys())
        desalt = set(desalt.keys())
        minimap2 = set(minimap2.keys())
        data_for_venn = [ultra, desalt, minimap2]
        total = len((ultra | desalt | minimap2 ))
    r = venn3(data_for_venn, ("uLTRA", "deSALT_GTF", "minimap2_GTF"), subset_label_formatter=lambda x: f"{(x/total):1.1%}")
    plt.savefig(os.path.join(outfolder, filename +".pdf"))
    plt.clf()

def get_success_regions(data_for_success_cases, reads, outfolder):
    ultra_fsm_distribution, ds_fsm_distribution, mm2_fsm_distribution = data_for_success_cases

    ultra_unique = set(ultra_fsm_distribution.keys()) - (set(ds_fsm_distribution.keys()) | set(mm2_fsm_distribution.keys()))
    outfile = open(os.path.join(outfolder, "ultra_unique_FSMs.csv"), "w")
    fa_outfile = open(os.path.join(outfolder, "ultra_unique_FSMs.fa"), "w")

    interesting_success_cases = []
    for tr_id in ultra_unique:
        for acc in ultra_fsm_distribution[tr_id]:
            outfile.write("{0},{1}\n".format(tr_id, acc)) 
        interesting_success_cases.append( (tr_id, len(ultra_fsm_distribution[tr_id])) )
        # nr_fsm_reads = len(ultra_fsm_distribution[tr_id])
    print("TOTAL number of unique interesting_success_cases:", len(interesting_success_cases))
    print("TOTAL number of reads in unique interesting_success_cases:", sum([nr_reads for tr_id, nr_reads in interesting_success_cases]) )
    print("Printing some of the more abundant ones (over 10 reads):")
    data_for_hist = []
    for tr_id, nr_fsm_reads in sorted(interesting_success_cases, key=lambda x: x[1], reverse=True):
        data_for_hist.append(nr_fsm_reads)
        if nr_fsm_reads >= 10:
            print("interesting success case:", tr_id, nr_fsm_reads)
            for acc in ultra_fsm_distribution[tr_id]:
                seq, qual = reads[acc]
                fa_outfile.write(">{0}\n{1}\n".format(acc, seq)) 

    outfile.close()
    fa_outfile.close()
    return data_for_hist


def get_diff_regions(data_for_success_cases, reads, outfolder):
    ultra_fsm_distribution, ds_fsm_distribution, mm2_fsm_distribution = data_for_success_cases

    ultra_missed = (set(ds_fsm_distribution.keys()) & set(mm2_fsm_distribution.keys())) - set(ultra_fsm_distribution.keys())
    outfile = open(os.path.join(outfolder, "ultra_diff_than_mm2_and_desalt.csv"), "w")
    fa_outfile = open(os.path.join(outfolder, "ultra_diff_than_mm2_and_desalt.fa"), "w")
    for tr_id in ultra_missed:
        for acc in ds_fsm_distribution[tr_id]:
            outfile.write("{0},{1}\n".format(tr_id, acc)) 

        nr_fsm_reads = len(ds_fsm_distribution[tr_id])
        if nr_fsm_reads >= 10:
            print("interesting missed case:", tr_id, nr_fsm_reads)
            for acc in ds_fsm_distribution[tr_id]:
                seq,qual = reads[acc]
                fa_outfile.write(">{0}\n{1}\n".format(acc, seq)) 

    outfile.close()
    fa_outfile.close()

def get_overlap_venn(data1, data2, data3, flag1, flag2, flag3):
    '''
        assigns 1 : minimap2, 10: seSALT, 100: uLTRA
        label: 1: mm,  12: mm & ds, 123: mm & ds & ult, ...
    '''
    data1_labelled = set()
    for acc in data1:
        mm2_start, mm2_stop, data1_splice_choords = data1[acc]
        label = flag1
        if acc in data2:
            ds_start, ds_stop, data2_splice_choords = data2[acc]
            if is_overlapping(mm2_start, mm2_stop,  ds_start, ds_stop):
                label += flag2
            else:
                pass
        else:
            pass

        if acc in data3:
            ult_start, ult_stop, data3_splice_choords = data3[acc]
            if is_overlapping(mm2_start, mm2_stop,  ult_start, ult_stop):
                label += flag3
            else:
                pass
        else:
            pass
            
        data1_labelled.add( (acc, label))
    return data1_labelled


def get_splice_sites_match_venn(data1, data2, data3, flag1, flag2, flag3):
    '''
        assigns 1 : minimap2, 10: seSALT, 100: uLTRA
        label: 1: mm,  12: mm & ds, 123: mm & ds & ult, ...
    '''
    data1_labelled = set()
    for acc in data1:
        data1_start, data1_stop, data1_splice_choords = data1[acc]
        label = flag1
        if acc in data2:
            data2_start, data2_stop, data2_splice_choords = data2[acc]
            if data1_splice_choords == data2_splice_choords:
                label += flag2
            else:
                pass
        else:
            pass

        if acc in data3:
            data3_start, data3_stop, data3_splice_choords = data3[acc]
            if data1_splice_choords == data3_splice_choords:
                label += flag3
            else:
                pass
        else:
            pass
            
        data1_labelled.add( (acc, label))
    return data1_labelled


def fill_aligner_data(read_alignments, is_ultra = False):
    positions = {}
    categories = defaultdict(dict)
    aln = {}
    if not is_ultra:
        genomic = set()

    for acc in read_alignments:
        annot, splice_choords, chr_id, start, stop, is_exonic = read_alignments[acc][7], read_alignments[acc][9], read_alignments[acc][11], read_alignments[acc][12], read_alignments[acc][13], read_alignments[acc][15]
        if is_exonic == '0':
            genomic.add(acc)
        elif is_exonic == '1':
            positions[acc] = (start, stop, splice_choords)
            categories[annot][acc] = (start, stop, splice_choords)
        aln[acc] = annot

    if is_ultra:
        return positions, categories, aln
    else:
        return positions, categories, aln, genomic


def fill_ultra_data(read_alignments, is_genomic):
    ult_positions = {}
    ultra_unaligned = set()
    ultra_genomic_reads_categories =  defaultdict(int)
    ultra_exonic_reads_categories =  defaultdict(dict)
    suspicious_fsms = set()
    for acc in reads_isonalign:
        ia_annot, splice_choords, ia_chr, ia_start, ia_stop, ia_is_exonic = reads_isonalign[acc][7], read_alignments[acc][9], reads_isonalign[acc][11], reads_isonalign[acc][12], reads_isonalign[acc][13], reads_isonalign[acc][15]
        if acc in is_genomic:
            ultra_genomic_reads_categories[ia_annot] += 1
            if ia_annot == "FSM":
                suspicious_fsms.add(acc)
        else:
            ultra_exonic_reads_categories[ia_annot][acc] = (ia_start, ia_stop, splice_choords)

        if ia_annot == 'unaligned':
            ultra_unaligned.add(acc)
        if ia_is_exonic == '1':
            ult_positions[acc] = (ia_start, ia_stop, splice_choords)

    return ult_positions, ultra_genomic_reads_categories, ultra_unaligned, suspicious_fsms, ultra_exonic_reads_categories


def get_mapping_location_concordance(reads_isonalign, reads_minimap2, reads_desalt, reads):

    mm_positions, mm_categories, mm_aln, mm_genomic =  fill_aligner_data(reads_minimap2, is_ultra = False)
    ds_positions, ds_categories, ds_aln, ds_genomic =  fill_aligner_data(reads_desalt, is_ultra = False)

    # mm_positions = {}
    # mm_categories = defaultdict(dict)
    # mm_aln = {}
    # mm_genomic = set()
    # for acc in reads_minimap2:
    #     mm2_annot, mm_splice_choords, mm2_chr, mm2_start, mm2_stop, mm_is_exonic = reads_minimap2[acc][7], reads_minimap2[acc][9], reads_minimap2[acc][11], reads_minimap2[acc][12], reads_minimap2[acc][13], reads_minimap2[acc][15]
    #     if mm_is_exonic == '0':
    #         mm_genomic.add(acc)
    #     elif mm_is_exonic == '1':
    #         mm_positions[acc] = (mm2_start, mm2_stop)
    #         mm_categories[mm2_annot][acc] = (mm2_start, mm2_stop, mm_splice_choords)
    #     mm_aln[acc] = mm2_annot
    
    # ds_positions = {}
    # ds_categories = defaultdict(dict)
    # ds_aln = {}
    # ds_genomic = set()
    # for acc in reads_desalt:
    #     ds_annot, ds_chr, ds_start, ds_stop, ds_is_exonic = reads_desalt[acc][7], reads_desalt[acc][11], reads_desalt[acc][12], reads_desalt[acc][13], reads_desalt[acc][15]
    #     if ds_is_exonic == '0':
    #         ds_genomic.add(acc)
    #     elif ds_is_exonic == '1':
    #         ds_positions[acc] = (ds_start, ds_stop)

    #     ds_aln[acc] = ds_annot

    is_genomic = ds_genomic & mm_genomic

    # get the uLTRA categories for likely genomic reads
    ult_positions, ultra_genomic_reads_categories, ultra_unaligned, suspicious_fsms, ult_categories = fill_ultra_data(reads_isonalign, is_genomic)

    # ult_positions = {}
    # ultra_unaligned = set()
    # ultra_genomic_reads_categories =  defaultdict(int)
    # suspicious_fsms = set()
    # for acc in reads_isonalign:
    #     ia_annot, ia_chr, ia_start, ia_stop, ia_is_exonic = reads_isonalign[acc][7], reads_isonalign[acc][11], reads_isonalign[acc][12], reads_isonalign[acc][13], reads_isonalign[acc][15]
    #     if acc in is_genomic:
    #         ultra_genomic_reads_categories[ia_annot] += 1
    #         if ia_annot == "FSM":
    #             suspicious_fsms.add(acc)
    #     if ia_annot == 'unaligned':
    #         ultra_unaligned.add(acc)
    #     if ia_is_exonic == '1':
    #         ult_positions[acc] = (ia_start, ia_stop)


    data_for_venn = {}

    # All reads

    mm_ovl_venn = get_overlap_venn(mm_positions, ds_positions, ult_positions, 1,10,100)
    ds_ovl_venn = get_overlap_venn(ds_positions, mm_positions, ult_positions, 10,1,100)
    ult_ovl_venn = get_overlap_venn(ult_positions, mm_positions, ds_positions, 100,1,10)
    tot = set(set(mm_positions.keys()) | set(ds_positions.keys()) | set(ult_positions.keys()))
    print("total exonic reads evaluated for concordance:", len(tot))
    tot2 = set(mm_ovl_venn | ds_ovl_venn | ult_ovl_venn)
    print("total union of unique entries in the concordance venn diagram:", len(tot2))
    tot3 = set(mm_ovl_venn & ds_ovl_venn & ult_ovl_venn)
    print("total shared between all unique entries:", len(tot3))
    data_for_venn["ALN"] = [ult_ovl_venn, ds_ovl_venn, mm_ovl_venn]

    # Separated into category 
    label_types = ["NNC", "NO_SPLICE", "NIC", "ISM", "FSM"]
    for label_type in label_types:
        mm_reads = mm_categories[label_type]
        ds_reads = ds_categories[label_type]
        ult_reads = ult_categories[label_type]
        if label_type in set(["NIC", "ISM", "FSM"]):
            mm_ovl_venn = get_splice_sites_match_venn(mm_reads, ds_reads, ult_reads, 1, 10, 100)
            ds_ovl_venn = get_splice_sites_match_venn(ds_reads, mm_reads, ult_reads, 10,1,100)
            ult_ovl_venn = get_splice_sites_match_venn(ult_reads, mm_reads, ds_reads, 100,1,10)
        else:
            mm_ovl_venn = get_overlap_venn(mm_reads, ds_reads, ult_reads, 1,10,100)
            ds_ovl_venn = get_overlap_venn(ds_reads, mm_reads, ult_reads, 10,1,100)
            ult_ovl_venn = get_overlap_venn(ult_reads, mm_reads, ds_reads, 100,1,10)
        
        data_for_venn[label_type] = [ult_ovl_venn, ds_ovl_venn, mm_ovl_venn]
        print()
        print(label_type)
        tot = set(set(mm_reads.keys()) | set(ds_reads.keys()) | set(ult_reads.keys()))
        print("total exonic reads evaluated for concordance:", len(tot))
        tot2 = set(mm_ovl_venn | ds_ovl_venn | ult_ovl_venn)
        print("total union of unique entries in the concordance venn diagram:", len(tot2))
        tot3 = set(mm_ovl_venn & ds_ovl_venn & ult_ovl_venn)
        print("total shared between all unique entries:", len(tot3))

    # label_types = ["NIC", "ISM", "FSM"]
    # for label_type in label_types:
    #     mm_reads = mm_categories[label_type]
    #     ds_reads = ds_categories[label_type]
    #     ult_reads = ult_categories[label_type]

    #     mm_ovl_venn = get_splice_sites_match_venn(mm_reads, ds_reads, ult_reads, 1, 10, 100)
    #     ds_ovl_venn = get_splice_sites_match_venn(ds_reads, mm_reads, ult_reads, 10,1,100)
    #     ult_ovl_venn = get_splice_sites_match_venn(ult_reads, mm_reads, ds_reads, 100,1,10)

    #     data_for_venn[label_type] = [ult_ovl_venn, ds_ovl_venn, mm_ovl_venn]

    print("suspicious_fsm reads", len(suspicious_fsms), "20 first:", list(suspicious_fsms)[:20])
    print("ULTRA categories of likely genomic reads:", ultra_genomic_reads_categories)

    mm_categories = defaultdict(int)
    ds_categories = defaultdict(int)
    unaligned_fsms = set()
    for acc in ultra_unaligned:
        mm_categories[ mm_aln[acc] ] += 1
        ds_categories[ ds_aln[acc] ] += 1

        if mm_aln[acc] == ds_aln[acc] == 'FSM':
            unaligned_fsms.add(acc)

    print("Minimap2 categories of ultra unaligned reads:", mm_categories)
    print("Desalt categories of ultra unaligned reads:", ds_categories)
    print("potential missed fsm reads", len(unaligned_fsms), "20 first:", list(unaligned_fsms)[:20])


    outfile = open(os.path.join(outfolder, "unaligned_FSMs.csv"), "w")
    fa_outfile = open(os.path.join(outfolder, "unaligned_FSMs.fa"), "w")
    for acc in unaligned_fsms:
        seq,qual = reads[acc]
        fa_outfile.write(">{0}\n{1}\n".format(acc, seq)) 

    outfile.close()
    fa_outfile.close()

    return data_for_venn



def get_ultra_categories_of_missed_likely_fsm_reads(data_for_venn, reads_isonalign, reads):
    ultra, desalt, minimap2 = data_for_venn
    missed_fsm = (desalt & minimap2) - ultra
    missed_fsm_read_acc = {acc : transcript_fsm_id for acc, transcript_fsm_id in missed_fsm}

    outfile = open(os.path.join(outfolder, "all_missed_fsm_reads.csv"), "w")
    fa_outfile = open(os.path.join(outfolder, "all_missed_fsm_reads.fa"), "w")
    ultra_categories =  defaultdict(int)
    for acc in reads_isonalign:
        ia_annot, ia_chr, ia_start, ia_stop, ia_is_exonic = reads_isonalign[acc][7], reads_isonalign[acc][11], reads_isonalign[acc][12], reads_isonalign[acc][13], reads_isonalign[acc][15]
        if acc in missed_fsm_read_acc:
            tr_id = missed_fsm_read_acc[acc]
            ultra_categories[ia_annot] += 1
            seq, qual = reads[acc]
            fa_outfile.write(">{0}\n{1}\n".format(acc + "_" + tr_id, seq)) 
            outfile.write("{0}\t{1}\n".format(acc,tr_id))
    outfile.close()
    fa_outfile.close()

    print("ULTRA categories of likely FMS reads (predicted by both mm2 and deSALT):", ultra_categories)


def get_unique_NIC(reads_isonalign, reads_minimap2, reads_desalt, reads, outfolder):
    ultra_NIC = defaultdict(set)
    for acc in reads_isonalign:
        # reads_isonalign[acc] =  (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_exonic) 
        ia_annot, splice_choords, ia_chr  = reads_isonalign[acc][7], reads_isonalign[acc][9], reads_isonalign[acc][11]
        if ia_annot == 'NIC':
            ultra_NIC[(splice_choords, ia_chr)].add(acc)

    minimap2_NIC = defaultdict(set)
    for acc in reads_minimap2:
        mm2_annot, splice_choords, ia_chr  = reads_minimap2[acc][7], reads_minimap2[acc][9], reads_minimap2[acc][11]
        if mm2_annot == 'NIC':
            minimap2_NIC[(splice_choords, ia_chr)].add(acc)

    desalt_NIC = defaultdict(set)
    for acc in reads_desalt:
        ds_annot, splice_choords, ia_chr  = reads_desalt[acc][7], reads_desalt[acc][9], reads_desalt[acc][11]
        if ds_annot == 'NIC':
            desalt_NIC[(splice_choords, ia_chr)].add(acc)

    # find the NICs unique to uLTRA
    ultra_unique_NICs = set(ultra_NIC.keys()) - (set(minimap2_NIC.keys()) | set(desalt_NIC.keys()))


    interesting_cases = []
    for nic_id in ultra_unique_NICs:
        if len(ultra_NIC[nic_id]) >= 5:
            interesting_cases.append( (nic_id,len(ultra_NIC[nic_id])) )

    outfile = open(os.path.join(outfolder, "ultra_unique_NICs.csv"), "w")
    fa_outfile = open(os.path.join(outfolder, "ultra_unique_NICs.fa"), "w")
    print("TOTAL uLTRA distinct  NICS:", len(ultra_NIC))
    print("TOTAL uLTRA NIC reads:", sum([len(reads) for nic_id, reads in  ultra_NIC.items()] ))

    for (nic_id, nr_reads) in sorted(interesting_cases, key=lambda x: x[1], reverse=True):
        print("interesting unique NIC case:", nic_id, len(ultra_NIC[nic_id]))

        exons = nic_id[0].split("-")
        exons = exons[1:-1]
        contains_small_exon = False
        smallest_exon = 0
        if exons:
            for e in exons:
                e1,e2 = e.split(":")
                e1, e2 = int(e1), int(e2)
                if e2-e1 < 20:
                    contains_small_exon = True
                    smallest_exon = e2-e1
        if contains_small_exon:
            print("This case contains small exon, extra interesting!", "smallest_exon:", smallest_exon)

        mm2_preds = []
        for acc in ultra_NIC[nic_id]:
            mm2_preds.append(reads_minimap2[acc][7])
            seq, qual = reads[acc]
            fa_outfile.write(">{0}\n{1}\n".format(acc, seq)) 
            outfile.write("{0},{1}\n".format(nic_id, acc)) 
        print(mm2_preds)

    outfile.close()
    fa_outfile.close()
    return ultra_NIC



def plot_hist(x, outfolder, filename):
    ax = sns.distplot(x, bins=50, kde=False, rug=True)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel("Count")
    ax.set_xlabel("Number of FSM reads")
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder, filename +".pdf"))
    plt.savefig(os.path.join(outfolder, filename +".eps"))
    plt.close()


def main(args):

    reads_isonalign, reads_minimap2, reads_desalt = parse_differing_location_reads(args.csvfile)
    fsm_data_for_venn, data_for_success_cases = get_FSM_concordance( reads_isonalign, reads_minimap2, reads_desalt)
    venn(fsm_data_for_venn, args.outfolder, "reads_to_FSM_concordance")
    venn(data_for_success_cases, args.outfolder, "unique_FSM_concordance")

    reads = { acc.split()[0] : (seq, qual) for i, (acc, (seq, qual)) in enumerate(readfq(open(args.reads, 'r')))}
    print("Total reads", len(reads))

    
    # OVERLAP
    overlap_data_for_venn = get_mapping_location_concordance(reads_isonalign, reads_minimap2, reads_desalt, reads)

    for category in overlap_data_for_venn
        data = overlap_data_for_venn[category]
        venn(data, args.outfolder, "reads_{0}_concordance".format(category))

    # FSM
    get_ultra_categories_of_missed_likely_fsm_reads(fsm_data_for_venn, reads_isonalign, reads)
    data_for_hist = get_success_regions(data_for_success_cases, reads, args.outfolder)
    get_diff_regions(data_for_success_cases, reads, args.outfolder)
    get_unique_NIC(reads_isonalign, reads_minimap2, reads_desalt, reads, args.outfolder)
    plot_hist(data_for_hist, args.outfolder, "ultra_unique_found_FSM_distribution")

    # fq_outfile = open(os.path.join(args.outfolder, "diff_mapped.fq"), "w")
    # info_outfile = open(os.path.join(args.outfolder, "diff_mapped.csv"), "w")
    # for acc in diff_mapped:
    #     info = diff_mapped[acc]
    #     info_outfile.write(acc + "," + ",".join([str(i) for i in  info]) + "\n") 
    #     (seq, qual) = reads[acc]   
    #     fq_outfile.write("@{0}\n{1}\n{2}\n{3}\n".format(acc, seq, "+", qual))    

    # fq_outfile.close()
    # info_outfile.close()






if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('reads', type=str, help='Path to the read file')
    parser.add_argument('csvfile', type=str, help='Path to the csv file')
    parser.add_argument('outfolder', type=str, help='Output path of results')

    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)

