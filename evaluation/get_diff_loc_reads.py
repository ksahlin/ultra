
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


def parse_differing_location_reads(csv_file):
    reads_isonalign = {}
    reads_minimap2 = {}
    reads_desalt = {}
    for line in open(csv_file,'r'):
        #print(line)
        acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag = line.strip().split(",")
        if algorithm == 'uLTRA':
            reads_isonalign[acc] =  (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag) 
        if algorithm == 'minimap2':
            reads_minimap2[acc] = (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag) 
        if algorithm == 'deSALT':
            reads_desalt[acc] = (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag)

    differing_reads = defaultdict(set)
    for acc in reads_isonalign:
        if acc in reads_minimap2 and acc in reads_desalt:
            mm2_chr, mm2_start, mm2_stop = reads_minimap2[acc][11], reads_minimap2[acc][12], reads_minimap2[acc][13]
            ds_chr, ds_start, ds_stop = reads_desalt[acc][11], reads_desalt[acc][12], reads_desalt[acc][13]
            ia_chr, ia_start, ia_stop = reads_isonalign[acc][11], reads_isonalign[acc][12], reads_isonalign[acc][13]
            if mm2_chr == ds_chr and ( ( int(ds_start) <= int(mm2_start) <= int(ds_stop) )  or ( int(ds_start) <= int(mm2_stop) <= int(ds_stop) ) ): # they are overlapping
                if ia_chr != mm2_chr and not ( ( int(ds_start) <= int(ia_start) <= int(ds_stop) )  or ( int(ds_start) <= int(ia_stop) <= int(ds_stop) ) ): # not overlapping with isONalign
                    differing_reads[acc].add( ( "differing", mm2_chr, ds_start, ds_stop, mm2_start, mm2_stop,  reads_isonalign[acc]) )
        else:
            if acc not in reads_minimap2 and acc not in reads_desalt: 
                differing_reads[acc].add( ("unaligned_both", "-",  "-",  "-",  "-",  "-", reads_isonalign[acc]) )
            elif acc not in reads_minimap2:
                ds_chr, ds_start, ds_stop = reads_desalt[acc][11], reads_desalt[acc][12], reads_desalt[acc][13]
                differing_reads[acc].add( ("unaligned_mm2",mm2_chr, ds_start, ds_stop, "-", "-",  reads_isonalign[acc]) )
            elif acc not in reads_minimap2:
                mm2_chr, mm2_start, mm2_stop = reads_minimap2[acc][11], reads_minimap2[acc][12], reads_minimap2[acc][13]
                differing_reads[acc].add( ("unaligned_ds",mm2_chr,  "-", "-", mm2_start, mm2_stop,  reads_isonalign[acc]) )
    return differing_reads


def main(args):

    diff_mapped = parse_differing_location_reads(args.csvfile)
    reads = { acc.split()[0] : (seq, qual) for i, (acc, (seq, qual)) in enumerate(readfq(open(args.reads, 'r')))}
    print("Total reads", len(reads))


    fq_outfile = open(os.path.join(args.outfolder, "diff_mapped.fq"), "w")
    info_outfile = open(os.path.join(args.outfolder, "diff_mapped.csv"), "w")
    for acc in differing_reads:
        info = differing_reads[acc]
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

