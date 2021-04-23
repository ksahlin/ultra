
from __future__ import print_function
import os,sys
import argparse
import re
import errno
import itertools

import pickle

from collections import defaultdict


try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

from matplotlib_venn import venn3, venn3_circles, venn2

import numpy as np
import seaborn as sns
import pandas as pd

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


def parse_differing_splicing_reads(csv_file, outfolder):
    reads_ultra = {}
    reads_ultra_mm2 = {}
    reads_minimap2 = {}
    reads_minimap2_gtf = {}
    reads_desalt = {}
    reads_desalt_gtf = {}
    fsm_unique_for_ultra = 0
    for line in open(csv_file,'r'):
        #print(line)
        acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_genomic = line.strip().split(",")
        if algorithm == 'uLTRA':
            reads_ultra[acc] =  (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_genomic) 
        if algorithm == 'uLTRA_mm2':
            reads_ultra_mm2[acc] =  (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_genomic) 

        if algorithm == 'minimap2':
            reads_minimap2[acc] = (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_genomic) 
        if algorithm == 'minimap2_GTF':
            reads_minimap2_gtf[acc] = (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_genomic) 

        if algorithm == 'deSALT':
            reads_desalt[acc] = (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_genomic)
        if algorithm == 'deSALT_GTF':
            reads_desalt_gtf[acc] = (acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_genomic)

    desalt_fsm_distribution = defaultdict(int)
    desalt = set()
    for acc in reads_desalt:
        ds_annot, transcript_fsm_id = reads_desalt[acc][7], reads_desalt[acc][10]
        if ds_annot == 'FSM':
            desalt_fsm_distribution[transcript_fsm_id] += 1
            desalt.add( (acc, transcript_fsm_id) ) 

    desalt_gtf_fsm_distribution = defaultdict(int)
    desalt_gtf = set()
    for acc in reads_desalt_gtf:
        ds_annot, transcript_fsm_id = reads_desalt_gtf[acc][7], reads_desalt_gtf[acc][10]
        if ds_annot == 'FSM':
            desalt_gtf_fsm_distribution[transcript_fsm_id] += 1
            desalt_gtf.add( (acc, transcript_fsm_id) ) 


    mm2_fsm_distribution = defaultdict(int)
    minimap2 = set()
    for acc in reads_minimap2:
        mm2_annot, transcript_fsm_id = reads_minimap2[acc][7], reads_minimap2[acc][10]
        if mm2_annot == 'FSM':
            mm2_fsm_distribution[transcript_fsm_id] += 1
            minimap2.add( (acc, transcript_fsm_id) ) 

    mm2_gtf_fsm_distribution = defaultdict(int)
    minimap2_gtf = set()
    for acc in reads_minimap2_gtf:
        mm2_annot, transcript_fsm_id = reads_minimap2_gtf[acc][7], reads_minimap2_gtf[acc][10]
        if mm2_annot == 'FSM':
            mm2_gtf_fsm_distribution[transcript_fsm_id] += 1
            minimap2_gtf.add( (acc, transcript_fsm_id) ) 

    ultra_fsm_distribution = defaultdict(int)
    ultra = set()
    for acc in reads_ultra:
        ultra_annot, transcript_fsm_id = reads_ultra[acc][7], reads_ultra[acc][10]
        if ultra_annot == 'FSM':
            ultra_fsm_distribution[transcript_fsm_id] += 1
            ultra.add( (acc, transcript_fsm_id) ) 

    ultra_mm2_fsm_distribution = defaultdict(int)
    ultra_mm2 = set()
    for acc in reads_ultra_mm2:
        ultra_annot, transcript_fsm_id = reads_ultra_mm2[acc][7], reads_ultra_mm2[acc][10]
        if ultra_annot == 'FSM':
            ultra_mm2_fsm_distribution[transcript_fsm_id] += 1
            ultra_mm2.add( (acc, transcript_fsm_id) ) 



    print("desalt_gtf_fsm_distribution:", sorted(desalt_gtf_fsm_distribution.items(), key = lambda x: x[1]))        
    print("mm2_gtf_fsm_distribution:", sorted(mm2_gtf_fsm_distribution.items(), key = lambda x: x[1]))        
    print("ultra_fsm_distribution:", sorted(ultra_fsm_distribution.items(), key = lambda x: x[1]))        

    print("Ultra nr unique isoforms mapped to:", len(ultra_fsm_distribution))
    print("Ultra_mm2 nr unique isoforms mapped to:", len(ultra_mm2_fsm_distribution))

    print("minimap2 nr unique isoforms mapped to:", len(mm2_fsm_distribution))
    print("minimap2_gtf nr unique isoforms mapped to:", len(mm2_gtf_fsm_distribution))

    print("Desalt nr unique isoforms mapped to:", len(desalt_fsm_distribution))
    print("Desalt_gtf nr unique isoforms mapped to:", len(desalt_gtf_fsm_distribution))

    print(len(ultra), len(desalt_gtf), len(minimap2_gtf))

    # a = ultra
    # b = desalt_gtf
    # c = minimap2_gtf
    a_not_b_c = ultra - (desalt_gtf | minimap2_gtf)
    b_not_a_c = desalt_gtf - (ultra | minimap2_gtf)
    c_not_a_b = minimap2_gtf - (ultra | desalt_gtf)
    a_b_not_c = (ultra & desalt_gtf) - minimap2_gtf
    a_c_not_b = (ultra & minimap2_gtf) - desalt_gtf
    b_c_not_a = (desalt_gtf & minimap2_gtf) - ultra
    a_b_c = ultra & desalt_gtf & minimap2_gtf

    print("BAD:")
    print("desalt_gtf and minimap2_gtf:", len(b_c_not_a))
    print("desalt_gtf unique:", len(b_not_a_c))
    print("minimap2_gtf unique:", len(c_not_a_b))
    print()
    print("GOOD:")
    print("Ultra and desalt_gtf:", len(a_b_not_c))
    print("Ultra and minimap2_gtf:", len(a_c_not_b))
    print()
    
    print("NEUTRAL")
    print("ultra unique:", len(a_not_b_c))
    print("In all", len(a_b_c))

    print("ALL FSM READS:", len( (ultra | desalt_gtf | minimap2_gtf )) )
    # return differing_reads

    return b_c_not_a, [ultra, desalt_gtf, minimap2_gtf], [ultra_fsm_distribution, ultra_mm2_fsm_distribution, \
                        desalt_fsm_distribution, desalt_gtf_fsm_distribution, \
                        mm2_fsm_distribution, mm2_gtf_fsm_distribution]


def venn(data_for_venn, outfolder):
    ultra, desalt, minimap2 = data_for_venn
    total = len((ultra | desalt | minimap2 ))
    r = venn3(data_for_venn, ("uLTRA_mm2", "deSALT_GTF", "minimap2_GTF"), subset_label_formatter=lambda x: f"{(x/total):1.1%}")
    plt.savefig(os.path.join(outfolder, "sirv_venn.pdf"))
    plt.clf()


def plot_nr_of_isoforms(data_for_mapping_bias, outfolder ):
    m = {0:'uLTRA', 1:'uLTRA_mm2', 2:'deSALT', 3:'deSALT_GTF', 4: 'minimap2', 5: 'minimap2_GTF'}
    g = {}
    # testdata
    # desalt_fsm_distribution = [('SIRV304', 4), ('SIRV702', 53), ('SIRV503', 85), ('SIRV301', 109), ('SIRV705', 514), ('SIRV613', 521), ('SIRV306', 544), ('SIRV103', 631), ('SIRV101', 666), ('SIRV703', 675), ('SIRV303', 1039), ('SIRV201', 1105), ('SIRV106', 1223), ('SIRV107', 1248), ('SIRV302', 1489), ('SIRV202', 1722), ('SIRV405', 1756), ('SIRV409', 1866), ('SIRV204', 1876), ('SIRV307', 2584), ('SIRV404', 4384), ('SIRV610', 4452), ('SIRV509', 4995), ('SIRV203', 6433), ('SIRV604', 6529), ('SIRV305', 6902), ('SIRV601', 7136), ('SIRV612', 7205), ('SIRV611', 8389), ('SIRV608', 11188), ('SIRV105', 13224), ('SIRV605', 15352), ('SIRV309', 16423), ('SIRV102', 16624), ('SIRV406', 16899), ('SIRV616', 17395), ('SIRV607', 24212), ('SIRV410', 24415), ('SIRV308', 25880), ('SIRV310', 29637), ('SIRV606', 31014), ('SIRV507', 49303), ('SIRV615', 67949), ('SIRV602', 70150), ('SIRV109', 71720), ('SIRV614', 77764), ('SIRV609', 91320)]
    # mm2_fsm_distribution = [('SIRV702', 71), ('SIRV503', 113), ('SIRV301', 114), ('SIRV613', 490), ('SIRV306', 580), ('SIRV103', 639), ('SIRV101', 678), ('SIRV705', 701), ('SIRV703', 706), ('SIRV201', 1188), ('SIRV303', 1212), ('SIRV107', 1263), ('SIRV106', 1274), ('SIRV302', 1725), ('SIRV405', 1811), ('SIRV204', 1839), ('SIRV202', 1851), ('SIRV409', 1940), ('SIRV307', 2945), ('SIRV610', 3750), ('SIRV404', 3974), ('SIRV604', 5984), ('SIRV203', 6155), ('SIRV601', 6431), ('SIRV612', 6599), ('SIRV305', 7082), ('SIRV611', 8403), ('SIRV608', 9985), ('SIRV605', 13966), ('SIRV105', 14933), ('SIRV616', 15240), ('SIRV309', 16587), ('SIRV102', 17024), ('SIRV406', 17090), ('SIRV607', 24139), ('SIRV410', 24661), ('SIRV308', 25823), ('SIRV606', 27173), ('SIRV310', 29674), ('SIRV507', 57275), ('SIRV602', 64702), ('SIRV614', 66636), ('SIRV615', 70215), ('SIRV109', 71483), ('SIRV609', 89736)]
    # ultra_fsm_distribution = [('SIRV702', 115), ('SIRV301', 167), ('SIRV501', 498), ('SIRV306', 667), ('SIRV613', 722), ('SIRV103', 813), ('SIRV101', 933), ('SIRV502', 1144), ('SIRV705', 1159), ('SIRV703', 1204), ('SIRV106', 1569), ('SIRV107', 1631), ('SIRV303', 1676), ('SIRV405', 1698), ('SIRV201', 1843), ('SIRV302', 2123), ('SIRV704', 2331), ('SIRV204', 2403), ('SIRV304', 2766), ('SIRV409', 2820), ('SIRV202', 2994), ('SIRV708', 3786), ('SIRV706', 3978), ('SIRV510', 4402), ('SIRV307', 4552), ('SIRV509', 4853), ('SIRV610', 5045), ('SIRV503', 6029), ('SIRV404', 6226), ('SIRV508', 7002), ('SIRV203', 8167), ('SIRV305', 8724), ('SIRV601', 9258), ('SIRV611', 10158), ('SIRV612', 11435), ('SIRV604', 11755), ('SIRV608', 12465), ('SIRV408', 13204), ('SIRV309', 16725), ('SIRV102', 17289), ('SIRV406', 17532), ('SIRV616', 20300), ('SIRV605', 25051), ('SIRV308', 26816), ('SIRV410', 27740), ('SIRV105', 28098), ('SIRV310', 30045), ('SIRV505', 30359), ('SIRV506', 30700), ('SIRV606', 34199), ('SIRV511', 34728), ('SIRV607', 36813), ('SIRV403', 37059), ('SIRV507', 71281), ('SIRV615', 74600), ('SIRV109', 80321), ('SIRV614', 83221), ('SIRV602', 98305), ('SIRV609', 123301)]
    data = pd.DataFrame(index=range(0), columns=['Algorithm', "Loci",  'SIRV_id', 'FSM_count'])
    i = 0
    loci = set()
    for k, method_dict in enumerate(data_for_mapping_bias):
        algorithm = m[k]
        for j, (sirv_id, fsm_count) in enumerate(method_dict.items()):
            gene = sirv_id[:5]
            loci.add(gene)
            data.loc[i] = [algorithm, gene, sirv_id, fsm_count]
            i += 1
    # print(data)
    for l in loci:
        d = data.loc[data['Loci'] == l]
        ax = sns.pointplot(x="SIRV_id", y="FSM_count", hue="Algorithm", data=d)
        plt.setp(ax.collections, alpha=.5) #for the markers
        plt.setp(ax.lines, alpha=.5)       #for the lines
        ax.set(ylabel="FSM count")
        ax.set(xlabel="SIRV ID")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=50)
        plt.tight_layout()
        plt.yscale('log')
        plt.savefig(os.path.join(outfolder, "sirv_counts_{0}.pdf".format(l)))
        plt.clf()      

    # g = sns.catplot(x="SIRV_id", y="FSM_count", hue="Algorithm", kind="point", col="Loci", col_wrap=3, data=data)
    # g.set_ylabels("FSM count")
    # g.set_xlabels("SIRV ID")
    # g.set_titles("{col_name} {col_var}")
    # g.set_xticklabels(rotation=90)
    # plt.yscale('log')
    # plt.savefig(os.path.join(outfolder, "sirv_counts.pdf"))


def main(args):

    desalt_and_minimap_unique, data_for_venn, data_for_mapping_bias = parse_differing_splicing_reads(args.csvfile, args.outfolder)
    reads = { acc.split()[0] : (seq, qual) for i, (acc, (seq, qual)) in enumerate(readfq(open(args.reads, 'r')))}
    # print("Total reads", len(reads))
    venn(data_for_venn, args.outfolder)
    plot_nr_of_isoforms(data_for_mapping_bias, args.outfolder )

    fq_outfile = open(os.path.join(args.outfolder, "desalt_and_minimap_unique.fq"), "w")
    info_outfile = open(os.path.join(args.outfolder, "desalt_and_minimap_unique.csv"), "w")
    for (acc, tr_id) in desalt_and_minimap_unique:
        info_outfile.write(acc + "," + tr_id + "\n") 
        (seq, qual) = reads[acc]   
        fq_outfile.write("@{0}\n{1}\n{2}\n{3}\n".format(acc, seq, "+", qual))    

    fq_outfile.close()
    info_outfile.close()






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

