
from __future__ import print_function
import os,sys
import argparse
import re
import errno
import itertools

import pickle

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
        acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_genomic = line.strip().split(",")
        if algorithm == 'uLTRA':
            reads_isonalign[acc] =  (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_genomic) 
        if algorithm == 'minimap2':
            reads_minimap2[acc] = (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_genomic) 
        if algorithm == 'deSALT':
            reads_desalt[acc] = (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_genomic)
    return reads_isonalign, reads_minimap2, reads_desalt


def get_FSM_concordance(reads_isonalign, reads_minimap2, reads_desalt):
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


def venn(data_for_venn, outfolder):
    ultra, desalt, minimap2 = data_for_venn
    total = len((ultra | desalt | minimap2 ))
    r = venn3(data_for_venn, ("uLTRA", "deSALT", "minimap2"), subset_label_formatter=lambda x: f"{(x/total):1.1%}")
    plt.savefig(os.path.join(outfolder, "fsm_concordance.pdf"))
    plt.clf()

def get_success_regions(data_for_success_cases, outfolder):
    ultra_fsm_distribution, ds_fsm_distribution, mm2_fsm_distribution = data_for_success_cases

    ultra_unique = set(ultra_fsm_distribution.keys()) - (set(ds_fsm_distribution.keys()) | set(mm2_fsm_distribution.keys()))
    outfile = open(os.path.join(outfolder, "success.csv"), "w")
    for tr_id in ultra_unique:
        for acc in ultra_fsm_distribution[tr_id]:
            outfile.write("{0},{1}\n".format(tr_id, acc)) 
    outfile.close()

def get_mapping_location_concordance(reads_isonalign, reads_minimap2, reads_desalt):

    differing_reads = defaultdict(set)
    mm_aln = defaultdict(set)
    mm_genomic = set()
    for acc in reads_minimap2:
        mm2_annot, mm2_chr, mm2_start, mm2_stop, mm_is_genomic = reads_minimap2[acc][7], reads_minimap2[acc][11], reads_minimap2[acc][12], reads_minimap2[acc][13], reads_minimap2[acc][16]
        if mm_is_genomic == '1':
            mm_genomic.add(acc)

    ds_genomic = set()
    for acc in reads_desalt:
        ds_annot, ds_chr, ds_start, ds_stop ds_is_genomic = reads_desalt[acc][7], reads_desalt[acc][11], reads_desalt[acc][12], reads_desalt[acc][13], reads_desalt[acc][16]
        if ds_is_genomic == '1':
            ds_genomic.add(acc)

    is_genomic = ds_genomic & mm_genomic
    ultra_aln = defaultdict(set)
    for acc in reads_isonalign:
        ia_annot, ia_chr, ia_start, ia_stop, ia_is_genomic = reads_isonalign[acc][7], reads_isonalign[acc][11], reads_isonalign[acc][12], reads_isonalign[acc][13], reads_isonalign[acc][16]

        if acc in is_genomic:
            ultra_aln[ia_annot] += 1

    print("Categories of likely genomic reads:", ultra_aln)
    # categories:
    #  genomic/exonic
    # What are the venn diagrams in overlapping locations for all the exonic reads

    # What is the category type for uLTRA for reads where minimap2 and deSALT agrees its a genomic read


    # interesting to infer: how many unaligned is genomic

        # if mm2_chr != 'unaligned' and  ds_chr != 'unaligned':
        #     if ia_chr == 'unaligned':
        #         differing_reads[acc].add( ( "unaligned_ultra",mm2_chr, ds_annot, ds_start, ds_stop, mm2_annot, mm2_start, mm2_stop,  reads_isonalign[acc]) )
        #     elif mm2_chr == ds_chr and is_overlapping(ds_start, ds_stop, mm2_start, mm2_stop): # they are overlapping
        #         if ia_chr != mm2_chr or not is_overlapping(ia_start, ia_stop, mm2_start, mm2_stop): # not overlapping with isONalign
        #             differing_reads[acc].add( ( "differing_pos", mm2_chr, ds_annot, ds_start, ds_stop, mm2_annot, mm2_start, mm2_stop,  reads_isonalign[acc]) )
        
        # # More detailed analysis for later
        # else:
        #     if mm2_chr == 'unaligned' and  ds_chr == 'unaligned':
        #         if ia_chr != 'unaligned': 
        #             differing_reads[acc].add( ("unaligned_mm2_ds", "-",  "-",  "-",  "-",  "-", reads_isonalign[acc]) )
        #     elif mm2_chr == 'unaligned':
        #         if ia_chr != 'unaligned': 
        #             ds_chr, ds_start, ds_stop = reads_desalt[acc][11], reads_desalt[acc][12], reads_desalt[acc][13]
        #             if not is_overlapping(ia_start, ia_stop, ds_start, ds_stop):
        #                 differing_reads[acc].add( ("unaligned_mm2_diff_ds",ds_chr, ds_annot, ds_start, ds_stop, "-", "-", "-", reads_isonalign[acc]) )
        #     elif ds_chr == 'unaligned':
        #         if ia_chr != 'unaligned': 
        #             mm2_chr, mm2_start, mm2_stop = reads_minimap2[acc][11], reads_minimap2[acc][12], reads_minimap2[acc][13]
        #             if not is_overlapping(ia_start, ia_stop, mm2_start, mm2_stop):
        #                 differing_reads[acc].add( ("unaligned_ds_diff_mm2", mm2_chr, '-', "-", "-", mm2_annot, mm2_start, mm2_stop, reads_isonalign[acc]) )
    # return differing_reads


def main(args):

    reads_isonalign, reads_minimap2, reads_desalt = parse_differing_location_reads(args.csvfile)
    data_for_venn, data_for_success_cases = get_FSM_concordance( reads_isonalign, reads_minimap2, reads_desalt)
    venn(data_for_venn, args.outfolder)
    get_success_regions(data_for_success_cases, args.outfolder)
    get_mapping_location_concordance(reads_isonalign, reads_minimap2, reads_desalt)


    reads = { acc.split()[0] : (seq, qual) for i, (acc, (seq, qual)) in enumerate(readfq(open(args.reads, 'r')))}
    print("Total reads", len(reads))


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

