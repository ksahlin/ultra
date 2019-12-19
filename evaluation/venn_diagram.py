
## Various plots from large table

import sys
import argparse
import os
import random
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")
# import matplotlib.pyplot as plt
# import matplotlib

import numpy as np
import seaborn as sns
import pandas as pd


def total_error_rate(input_csv, outfolder):

    indata = pd.read_csv(input_csv)
    # print(len(df))
    # indata = df.loc[df['q_acc'] == df['r_acc']]
    # print(len(indata))

    g = sns.catplot(x="read_type", y="error_rate", #col="Depth",
                data=indata,  #hue="read_type", hue_order= ["corrected", "original"],
                kind="violin", aspect=1)

    g.set(ylim=(0,15))
    g.set_ylabels("Error rate (%)")
    g.set_xlabels("Method")

    # ax = sns.boxplot(x="p", y=y, hue = "type", data=indata)
    # ax.set_ylim(0,15)
    # ax.set_ylabel("Error rate %")

    plt.savefig(os.path.join(outfolder, "total_error_rate.eps"))
    plt.savefig(os.path.join(outfolder, "total_error_rate.pdf"))
    plt.close()


def label_transcript(row):
   if row['fsm'] == 1 :
      return 'FSM'
   if row['nic'] == 1 :
      return 'NIC'
   if row['ism'] == 1:
      return 'ISM'
   if row['nnc']  == 1:
      return 'NNC'
   if row['no_splices']  == 1:
      return 'NO_SPLICE'

def splice_site_classification_plot(input_csv, outfolder):

    indata = pd.read_csv(input_csv)
    # indata['transcript_type'] = indata.apply (lambda row: label_transcript(row), axis=1)
    # print(len(df))
    # indata = df.loc[df['q_acc'] == df['r_acc']]
    # print(len(indata))

    g = sns.catplot(x="annotation", #col="Depth",
                data=indata,  hue="read_type", hue_order= ["uLTRA", "minimap2"],
                order= ["FSM", "ISM", "NIC", "NNC", 'NO_SPLICE'], kind="count", aspect=1)

    # g.set(ylim=(0,15))
    g.set_ylabels("Count")
    g.set_xlabels("Transcript type")

    # ax = sns.boxplot(x="p", y=y, hue = "type", data=indata)
    # ax.set_ylim(0,15)
    # ax.set_ylabel("Error rate %")

    plt.savefig(os.path.join(outfolder, "splice_site_classification.eps"))
    plt.savefig(os.path.join(outfolder, "splice_site_classification.pdf"))
    plt.close()

def number_splices_fsm(input_csv, outfolder):

    indata = pd.read_csv(input_csv)
    # indata['annotation'] = indata.apply (lambda row: label_transcript(row), axis=1)
    indata = indata[indata['annotation']=='FSM']
    g = sns.catplot(x="tot_splices", #col="Depth",
                data=indata,  hue="read_type", hue_order= ["uLTRA", "minimap2"], kind="count", aspect=1)
    # axes = g.axes
    g.set_ylabels("Count")
    g.set_xlabels("Number of splice sites")
    # axes.set_xticks(np.arange(0, 70, step=5) )
    # axes.set_xlim(xlim=(0, 70))
    # g.set_xlim(0,70)
    # g.set_xticks(np.arange(0, 70, step=5))
    # ax.set_ylabel("Error rate %")

    plt.savefig(os.path.join(outfolder, "nr_splices.eps"))
    plt.savefig(os.path.join(outfolder, "nr_splices.pdf"))
    plt.close()


from collections import defaultdict

def transpose_count(dct):
    d = defaultdict(int)
    for key1, value in dct.items():
        d[value]  += 1
    return d

def venn_diagram(input_csv, outfolder):

    indata = pd.read_csv(input_csv)
    correct_per_method = {}
    for method in ['uLTRA', "minimap2", "deSALT", "deSALT_GTF", "Graphmap2", "Graphgmap2_GTF"]:
        correct_per_method[method] = {}
        method_alnmnts = indata[indata['alignment_method']== method ]

        # method_alnmnts = indata[indata['alignment_method']== method ]
        
    # orig = indata[indata['read_type']=='original']
    # corr = indata[indata['read_type']=='corrected']

    # print('orig:', orig['transcript_fsm_id'].nunique())
    # print('corr', corr['transcript_fsm_id'].nunique())
    # print(set(orig['transcript_fsm_id'].unique()))
    # print(set(corr['transcript_fsm_id'].unique()))

    orig_reads = pd.Series(orig.transcript_fsm_id.values,index=orig.acc).to_dict()
    corr_reads = pd.Series(corr.transcript_fsm_id.values,index=corr.acc).to_dict()

    all_fsm_orig = set(orig['transcript_fsm_id'].unique())
    all_fsm_corr = set(corr['transcript_fsm_id'].unique())
    print('In orig but not in corr', all_fsm_orig - all_fsm_corr)
    print('In corr but not in orig', all_fsm_corr - all_fsm_orig)

    all_fsm_orig = {x for x in all_fsm_orig if x==x}
    all_fsm_corr = {x for x in all_fsm_corr if x==x}
    print('Total unique in orig', len(all_fsm_orig))
    print('Total unique in corr', len(all_fsm_corr))

    # Categories: (0. not present/aligned in Corr, 0'. not present/aligned in Orig)  1. Both nAn, 2. Both FSM same, 2. Both FSM different, 3. Corr_FSM, Orig nAn, 3. Orig_FSM, Corr nAn,
    # print(orig_reads)
    not_present_corr = set(orig_reads) - set(corr_reads) 
    not_present_orig = set(corr_reads) - set(orig_reads) 
    in_both = set(corr_reads) & set(orig_reads)
    both_nan = [ acc for acc in in_both if orig_reads[acc] not in all_fsm_orig and corr_reads[acc] not in all_fsm_corr ]
    both_fsm_same = [ acc for acc in in_both if orig_reads[acc] == corr_reads[acc] and  corr_reads[acc] in all_fsm_corr ]
    both_fsm_different = [ acc for acc in in_both if orig_reads[acc] != corr_reads[acc] and  corr_reads[acc] in all_fsm_corr and  orig_reads[acc] in all_fsm_orig ]
    corr_fsm_orig_nan = [ acc for acc in in_both if orig_reads[acc] != corr_reads[acc] and  corr_reads[acc] in all_fsm_corr and  orig_reads[acc] not in all_fsm_orig ]
    orig_fsm_corr_nan = [ acc for acc in in_both if orig_reads[acc] != corr_reads[acc] and  corr_reads[acc] not in all_fsm_corr and  orig_reads[acc] in all_fsm_orig ]

    print()
    print('not_present_corr', len(not_present_corr))
    print('not_present_orig', len(not_present_orig))
    print('in_both', len(in_both))
    print('both_nan', len(both_nan))
    print('both_fsm_same', len(both_fsm_same))
    print('both_fsm_different', len(both_fsm_different))
    print('corr_fsm_orig_nan', len (corr_fsm_orig_nan))
    print('orig_fsm_corr_nan', len(orig_fsm_corr_nan))
    print()
    reads_per_transcript_orig = transpose_count(orig_reads)
    reads_per_transcript_corr = transpose_count(corr_reads)
    print("transcript_id,FSM_orig,FSM_corr")
    for transcript_id in set(all_fsm_orig | all_fsm_corr):
        print(transcript_id,reads_per_transcript_orig[transcript_id], reads_per_transcript_corr[transcript_id])

    print("sum",sum([reads_per_transcript_orig[tr_id] for tr_id in set(all_fsm_orig | all_fsm_corr)]), sum([reads_per_transcript_corr[tr_id] for tr_id in set(all_fsm_orig | all_fsm_corr)]) )



def main(args):
    
    sns.set(style="whitegrid")
    flatui = ["#2ecc71", "#e74c3c"] # https://chrisalbon.com/python/data_visualization/seaborn_color_palettes/
    sns.set_palette(flatui)    # total_error_rate(args.input_csv, args.outfolder)

    venn_diagram(args.input_csv, args.outfolder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('input_csv', type=str, help='Path to all stats file')
    parser.add_argument('outfolder', type=str, help='Path to all stats file')
    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)

