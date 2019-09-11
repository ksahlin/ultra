
from __future__ import print_function
import os,sys
import argparse
import re
import errno
import itertools
import random

from collections import defaultdict

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


def get_transcript_annotations(gff_file): 
    fn = gffutils.example_filename(gff_file)
    db = gffutils.create_db(fn, dbfn='test.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    db = gffutils.FeatureDB('test.db', keep_order=True)
    annotated_transcripts = defaultdict(lambda: defaultdict(set))
    for gene in db.features_of_type('gene'):
        for transcript in db.children(gene, featuretype='transcript', order_by='start'):
            # print(dir(transcript))
            transcript_seq = [(exon.start-1, exon.stop) for exon in db.children(transcript, featuretype='exon', order_by='start')]
            annotated_transcripts[gene.seqid][transcript.id] = transcript_seq

    return annotated_transcripts


def main(args):
    refs = { acc : seq for acc, (seq, qual) in readfq(open(args.refs, 'r'))}
    annotated_transcripts = get_transcript_annotations(args.gff_file)
    classes = ["FSM", "ISM", "NIC_comb", "NIC_nov" "NNC"]
    outfile = open(args.outfile, "w")
    for i in range(args.nr_reads):
        chr_id = random.choice( list(annotated_transcripts.keys()) )
        print(chr_id)
        transcript_id = random.choice( list(annotated_transcripts[chr_id]) )
        print("sampled", chr_id, transcript_id)
        transcript = annotated_transcripts[chr_id][transcript_id]
        read = "".join([refs[chr_id][start:stop] for start,stop in transcript])
        if args.fastq:
            outfile.write("@{0}_{1}_{2}\n{3}\n+\n{4}\n".format(i, chr_id, transcript_id, read, "I"*len(read)))
        else:
            outfile.write(">{0}_{1}_{2}\n{3}\n".format(i, chr_id, transcript_id, read))

    outfile.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Sample reads.")
    parser.add_argument('refs', type=str, help='Path to the refs file')
    parser.add_argument('gff_file', type=str, help='Path to the refs file')
    parser.add_argument('outfile', type=str, help='Output path of results')
    parser.add_argument('--nr_reads', type=int, default=100, help='Threchold for what is counted as varation/intron in alignment as opposed to deletion.')
    parser.add_argument('--fastq', action='store_true', help='Threchold for what is counted as varation/intron in alignment as opposed to deletion.')

    args = parser.parse_args()

    # outfolder = args.outfolder
    # if not os.path.exists(outfolder):
    #     os.makedirs(outfolder)
    main(args)

