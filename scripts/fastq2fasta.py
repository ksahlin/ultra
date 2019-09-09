
from __future__ import print_function
import os,sys
import argparse



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



def main(args):

    reads = { acc : seq for acc, (seq, qual) in readfq(open(args.fastq, 'r'))}
    outfile = open(args.fasta, 'w')

    for acc, seq in reads.items():
        outfile.write('>{0}\n{1}\n'.format(acc, seq)) 
    outfile.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Sample reads.")
    parser.add_argument('fastq', type=str, help='Path to the refs file')
    parser.add_argument('fasta', type=str, help='Path to the refs file')

    args = parser.parse_args()

    # outfolder = args.outfolder
    # if not os.path.exists(outfolder):
    #     os.makedirs(outfolder)
    main(args)

