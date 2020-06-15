
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


def parse_differing_splicing_reads(csv_file):
    reads_isonalign = {}
    # reads_minimap2 = {}
    reads_desalt = {}
    fsm_unique_for_ultra = 0
    for line in open(csv_file,'r'):
        #print(line)
        acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag, is_genomic = line.strip().split(",")
        if algorithm == 'uLTRA':
            reads_isonalign[acc] =  (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag, is_genomic) 
        if algorithm == 'deSALT_GTF':
            reads_desalt[acc] = (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag, is_genomic)

    differing_reads = defaultdict(set)
    desalt_FSM = 0
    ultra_FSM = 0
    for acc in reads_desalt:
        ds_annot, ds_chr, ds_start, ds_stop = reads_desalt[acc][7], reads_desalt[acc][11], reads_desalt[acc][12], reads_desalt[acc][13]
        ia_annot, ia_chr, ia_start, ia_stop = reads_isonalign[acc][7], reads_isonalign[acc][11], reads_isonalign[acc][12], reads_isonalign[acc][13]
        if ds_annot == 'FSM' and  ia_annot != 'FSM':
            differing_reads[acc].add( (ds_annot, ds_start, ds_stop, reads_isonalign[acc]) )

        elif ds_annot != 'FSM' and ia_annot == 'FSM':
            fsm_unique_for_ultra += 1  

        if ds_annot == 'FSM':
            desalt_FSM += 1
        if ia_annot == 'FSM':
            ultra_FSM += 1

    print("fsm_unique_for_ultra:", fsm_unique_for_ultra)        
    print("fsm_unique_for_desalt_GTF:", len(differing_reads))        

    print("fsm_tot_ultra:", ultra_FSM)        
    print("fsm_tot_desalt_GTF:", desalt_FSM)    

    return differing_reads


def main(args):

    diff_spliced = parse_differing_splicing_reads(args.csvfile)
    reads = { acc.split()[0] : (seq, qual) for i, (acc, (seq, qual)) in enumerate(readfq(open(args.reads, 'r')))}
    print("Total reads", len(reads))


    fq_outfile = open(os.path.join(args.outfolder, "diff_spliced.fq"), "w")
    info_outfile = open(os.path.join(args.outfolder, "diff_spliced.csv"), "w")
    for acc in diff_spliced:
        info = diff_spliced[acc]
        info_outfile.write(acc + "," + ",".join([str(i) for i in  info]) + "\n") 
        (seq, qual) = reads[acc]   
        fq_outfile.write("@{0}\n{1}\n{2}\n{3}\n".format(acc, seq, "+", qual))    

    fq_outfile.close()
    info_outfile.close()
    # detailed_results_outfile = open(os.path.join(args.outfolder, "results_per_read.csv"), "w")
    # detailed_results_outfile.write("acc,read_type,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag\n")






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

