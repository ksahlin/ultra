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

def correctness_per_exon_size(input_csv, outfolder):
    indata = pd.read_csv(input_csv)
    nr_reads = float(len(indata.loc[indata['alignment_algorithm'] == 'uLTRA'])) # count the reads in a given dataset
    print(nr_reads)

    indata = pd.read_csv(input_csv)
    ax = sns.lineplot(x="exon_size", y="fraction_correct", hue="alignment_algorithm", estimator=None, lw=1, data=indata)
    # g = sns.catplot(x="alignment_classification", #col="Depth",
    #             data=indata,  hue="alignment_algorithm", hue_order= ["uLTRA", "minimap2", "deSALT", "deSALT_GTF", "Graphmap2", "Graphmap2_GTF"],
    #             order= ["correct", "site_diff", "diff_exon_count", "diff_location", 'unaligned'], kind="count", aspect=1)

    # g.set(ylim=(0,15))
    # g.set_ylabels("Count")
    # g.set_xlabels("Alignment type")
    # g.set_xticklabels(rotation=20)
    # plt.yscale('log')
    # ax = g.ax
    # for p in ax.patches:
    #     # ax.annotate('%{:.1f}'.format(p.get_height()), (p.get_x()+0.1, p.get_height()+50))
    #     ax.annotate('{:.2f}%'.format(100*p.get_height()/1000000.0), (p.get_x()+0.01, p.get_height()+1000), rotation=90, fontsize = 'x-small' )
    # # ax = sns.boxplot(x="p", y=y, hue = "type", data=indata)
    # # ax.set_ylim(0,15)
    ax.set_ylabel("Fraction correct")
    ax.set_xlabel("Exon size")
    # (g.set_axis_labels("Alignment type", "Count").set_xticklabels(["Correct", "Inexact", "Exon diff", "Incorrect", "Unaligned"]))
    plt.savefig(os.path.join(outfolder, "correct_per_exon_size.eps"))
    plt.savefig(os.path.join(outfolder, "correct_per_exon_size.pdf"))
    plt.close()



def main(args):
    
    sns.set(style="whitegrid")
    flatui = ["#2ecc71", "#e74c3c"] # https://chrisalbon.com/python/data_visualization/seaborn_color_palettes/
    sns.set_palette(flatui)    # total_error_rate(args.input_csv, args.outfolder)
    correctness_per_exon_size(args.input_csv, args.outfolder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('input_csv', type=str, help='Path to all stats file')
    parser.add_argument('outfolder', type=str, help='Path to all stats file')

    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)