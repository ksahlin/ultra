import argparse
import sys

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




def main(args):
    if args.keep_only:
        keep_only = set(args.keep_only)
    else:
        keep_only = set()
    
    if args.remove_refs:
        refs_remove = set(args.remove_refs)
    else: 
        refs_remove = set()

    if args.outfile:
        outfile = open(args.outfile,'w')
    else:
        outpath = args.fasta_file+'_filtered.fa'
        outfile = open(outpath,'w')

 
    for acc,(seq,_) in readfq(open(args.fasta_file, 'r')):
        if keep_only:
            if acc in keep_only:
                outfile.write(">{0}\n{1}\n".format(acc, seq))
                print("Kept:", acc)
                continue  
            else:
                print("Filtered:", acc)
                continue

        if acc in refs_remove:
            print("Filtered:", acc)
            continue

        if len(seq) >= args.min_size:
            print("Wrote to file:", acc)
            outfile.write(">{0}\n{1}\n".format(acc, seq))
    outfile.close()



if __name__ == '__main__':

    parser = argparse.ArgumentParser("Parses a fasta file and output all sequences longer than min_size to fasta format to /path/to/fasta_file_filtered.fa.")
    parser.add_argument('fasta_file', type=str, help='Fasta file. ')
    parser.add_argument('--min_size', type=int, default =0, help='Min size of sequences.')
    parser.add_argument('--outfile', type=str, default="", help='Path to file to write to. default is [INPUT FILE NAME]_filtered.fa')
    parser.add_argument('--remove_refs', type=str, nargs='+', default =[], help='Provide a list of accessions to be removed, e.g., "chr1, chr2,chr3...".')
    parser.add_argument('--keep_only', type=str, nargs='+', default =[], help='Provide a list of accessions to be kept, e.g., "chr1, chr2,chr3...".')

    args = parser.parse_args()
    print("removing:", args.remove_refs)
    main(args)
