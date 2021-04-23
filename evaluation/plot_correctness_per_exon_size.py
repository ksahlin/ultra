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
    ax = sns.lineplot(x="exon_size", y="fraction_correct", hue="alignment_algorithm",
                          hue_order= ["uLTRA", "uLTRA_mm2", "minimap2", "minimap2_GTF", "deSALT", "deSALT_GTF"],
                          # hue_order= ["uLTRA", "uLTRA_mm2", "minimap2", "minimap2_GTF", "deSALT", "deSALT_GTF", "Graphmap2", "Graphmap2_GTF"],
                           estimator=None, lw=1, data=indata)
    # g = sns.catplot(x="alignment_classification", #col="Depth",
    #             data=indata,  hue="alignment_algorithm", hue_order= ["uLTRA", "minimap2", "deSALT", "deSALT_GTF", "Graphmap2", "Graphmap2_GTF"],
    #             order= ["correct", "site_diff", "diff_exon_count", "diff_location", 'unaligned'], kind="count", aspect=1)

    # g.set(ylim=(0,15))
    # g.set_ylabels("Count")
    # g.set_xlabels("Alignment type")
    # g.set_xticklabels(rotation=20)
    plt.xscale('log')
    plt.legend(loc='lower right')
    # ax = g.ax
    # for p in ax.patches:
    #     # ax.annotate('%{:.1f}'.format(p.get_height()), (p.get_x()+0.1, p.get_height()+50))
    #     ax.annotate('{:.2f}%'.format(100*p.get_height()/1000000.0), (p.get_x()+0.01, p.get_height()+1000), rotation=90, fontsize = 'x-small' )
    # # ax = sns.boxplot(x="p", y=y, hue = "type", data=indata)
    # # ax.set_ylim(0,15)
    ax.set_ylabel("Fraction correct")
    ax.set_xlabel("Exon size")
    # (g.set_axis_labels("Alignment type", "Count").set_xticklabels(["Correct", "Inexact", "Exon diff", "Incorrect", "Unaligned"]))
    plt.savefig(os.path.join(outfolder, "correctness_per_exon_size.eps"))
    plt.savefig(os.path.join(outfolder, "correctness_per_exon_size.pdf"))
    plt.clf()
    plt.close()

def correctness_per_exon_size_binned(input_csv, outfolder):
    # sns.plt.clf()
    data = pd.read_csv(input_csv)
    with sns.plotting_context("paper", font_scale=1.2):
        plt.rc('text', usetex=False)
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize=10)
        plt.rc('ytick', labelsize=12)
        exon_sizes = [0,10,20,50,100,200,500,1000,5000,25000]
        bins = [i for i in exon_sizes ]
        labels = [ '{0}-{1}'.format(i+1,j) for i,j in zip(exon_sizes[:-1], exon_sizes[1:]) ] 
        # print(bins)
        # print(labels)
        exon_sizes_bins = pd.cut(data['exon_size'], bins, labels = labels )
        data["exon_bins"] = pd.Series(exon_sizes_bins, index=data.index)
        print(data)

        # data_grouped = data.groupby('exon_bins', 'alignment_algorithm')['nr_total', 'nr_corr'].sum()
        data_grouped = data.groupby(['exon_bins', 'alignment_algorithm'], as_index=False).agg({'nr_total':'sum', 'nr_corr' : 'sum' })
        print(data_grouped)
        # print(data.groupby('exon_bins')['nr_total', 'nr_corr','alignment_algorithm'].sum())
        data_grouped['fraction_correct'] = data_grouped['nr_corr']/data_grouped['nr_total']

        print(data_grouped)

        # data_grouped = data.groupby(["exon_bins"])
        # data_grouped.describe()
        # print(data_grouped.groupby(['exon_bins']).sum())
        # for key, item in data_grouped:
        #     print(data_grouped.get_group(key), "\n\n")
        # print(data.groupby(["exon_bins"]))
        # print(data.groupby(["exon_bins"])['nr_corr'].value_counts(normalize=True).mul(100))
        # data_tmp = data.groupby(["exon_bins"])['nr_corr'].value_counts(normalize=True).mul(100)
        # full_supp_data = [ {"exon_bins" : index[0], "percentage" : float(val)} for index, val in data_tmp.iteritems() if index[1] == "yes"]
        # full_supp = pd.DataFrame(full_supp_data)

        g = sns.barplot(x="exon_bins", y="fraction_correct", hue = 'alignment_algorithm', 
                        hue_order= ["uLTRA", "uLTRA_mm2", "minimap2", "minimap2_GTF", "deSALT", "deSALT_GTF"], data=data_grouped)

        # plt.xlabel('exon_size', fontsize=14)
        # plt.ylabel('Fraction correct',fontsize=16)


        plt.tick_params(rotation=20)
        plt.ylim(0, 1)
        g.legend(loc=4)
        g.set_ylabel("Fraction correct")
        g.set_xlabel("Exon size")
        plt.savefig(os.path.join(outfolder, "correctness_per_exon_size_binned.eps"))
        plt.savefig(os.path.join(outfolder, "correctness_per_exon_size_binned.pdf"))

        plt.clf()
        plt.close()




def main(args):
    
    sns.set(style="whitegrid")
    flatui = ["#2ecc71", "#e74c3c"] # https://chrisalbon.com/python/data_visualization/seaborn_color_palettes/
    sns.set_palette(flatui)    # total_error_rate(args.input_csv, args.outfolder)
    # correctness_per_exon_size(args.input_csv, args.outfolder)
    correctness_per_exon_size_binned(args.input_csv, args.outfolder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('input_csv', type=str, help='Path to all stats file')
    parser.add_argument('outfolder', type=str, help='Path to all stats file')

    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)