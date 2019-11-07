import os
import argparse
import numpy as np
import random
import math
from collections import defaultdict
import errno

'''
    Below function taken from https://github.com/lh3/readfq/blob/master/readfq.py
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



# log how many reads are simulated for each variant.


def simulate_read(i, transcript_acc, isoform ):
    # error_lvls = [0.8, 0.9, 0.92, 0.94, 0.96,0.98, 0.99, 0.995]
    error_lvls = [0.95]
    read = []
    qual = []
    was_del = False
    del_length = 0
    prev_base_pos = 0
    err = 0
    # print(isoform)
    for l, n in enumerate(isoform):
        if l < prev_base_pos + del_length:
            # print('skipping')
            continue
        p_correct_reading = random.choice(error_lvls)
        p_error = 1.0 - p_correct_reading
        # print(round(10**(-math.log(p_correct_reading))), p_error, -math.log(p_error,10)*10)

        r = random.uniform(0,1)
        if r > p_correct_reading:
            error = True
        else:
            error = False

        if error:
            # print(l)
            r = random.uniform(0,1)
            if r < 0.5: #deletion
                was_del = p_error
                del_length += 1
                r_del = random.uniform(0,1)
                while r_del >= 0.5:
                    del_length += 1
                    r_del = random.uniform(0,1)
                # print(l, del_length)
                err += del_length
            elif 0.5 <= r < 0.8:
                read.append(random.choice("ACGT"))
                qual.append( round(-math.log(p_error,10)*10) )
                prev_base_pos = l
                err += 1
            else:
                ins_len = 1
                read.append(n)
                qual.append( round(-math.log(p_error,10)*10) )
                r_ins = random.uniform(0,1)
                while r_ins >= 0.5:
                    read.append(random.choice("ACGT"))
                    r_ins = random.uniform(0,1)
                    qual.append( round(-math.log(0.7,10)*10) )
                    ins_len += 1
                prev_base_pos = l
                err += ins_len

        else:
            if was_del: # adding uncertainty from prevous deleted base
                read.append(n)
                qual.append( round(-math.log(was_del,10)*10) )
            else:
                read.append(n)
                qual.append( round(-math.log(p_error,10)*10) )
            prev_base_pos = l

            was_del = False
            del_length = 0

    read_seq = "".join([n for n in read])
    qual_seq = "".join([chr(q + 33) for q in qual])
    acc = str(transcript_acc) + "_" +  str(i) #= (read_seq, qual_seq)
    print(err, len(isoform), float(err)/ len(isoform))
    # print(read_seq)
    # print(qual_seq)

    return acc, read_seq, qual_seq


def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def main(args):
    sequence_transcripts = {acc: seq for acc, (seq, _) in readfq(open(args.sequence_material,"r")) }

    # just generate all numbers at once and draw from this 5x should be enough
    ont_reads = {}
    reads_generated_log = defaultdict(int)
    errors = []
    all_transctript_accessions = list(sequence_transcripts.keys())
    for i in range(args.read_count):
        acc = random.choice( all_transctript_accessions)
        transcript = sequence_transcripts[acc]
        read_acc, read, qual = simulate_read(i, acc, transcript)
        ont_reads[read_acc] = (read, qual)
        if i % 5000 == 0:
            print(i, "reads simulated.")


    # for acc, abundance in misc_functions.iteritems(reads_generated_log):
    #     args.logfile.write("{0}\t{1}\n".format(acc, abundance))

    # n = float(len(errors))
    # mu =  sum(errors) / n
    # sigma = (sum(list(map((lambda x: x ** 2 - 2 * x * mu + mu ** 2), errors))) / (n - 1)) ** 0.5
    # min_error = min(errors)
    # max_error = max(errors)
    # errors.sort()
    # if len(errors) %2 == 0:
    #     median_error = (errors[int(len(errors)/2)-1] + errors[int(len(errors)/2)]) / 2.0
    # else:
    #     median_error = errors[int(len(errors)/2)]

    # args.logfile.write("mean error: {0}, sd error:{1}, min_error:{2}, max_error:{3}, median_error:{4}\n".format(mu, sigma, min_error, max_error, median_error))

    outfile = open(args.outfile, "w")
    if args.fasta:
        for acc, (read_seq,qual_seq) in sorted(ont_reads.items(), key = lambda x: len(x[1]), reverse = True):
            outfile.write(">{0}\n{1}\n".format(acc, read_seq))
    else:
        for acc, (read_seq,qual_seq) in sorted(ont_reads.items(), key = lambda x: len(x[1]), reverse = True):
            outfile.write("@{0}\n{1}\n{2}\n{3}\n".format(acc, read_seq, "+", qual_seq))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate pacbio reads from a set of transcripts.")
    parser.add_argument('sequence_material', type=str, help='The fasta file with sequences to be sequenced.')
    parser.add_argument('outfile', type=str, help='Output path to fasta file')
    parser.add_argument('read_count', type=int, help='Number of reads to simulate.')
    parser.add_argument('--fasta', action="store_true", help='Output in fasta format')
    # parser.add_argument('config', type=str, help='config file')


    args = parser.parse_args()
    path_, file_prefix = os.path.split(args.outfile)
    mkdir_p(path_)
    args.logfile = open(os.path.join(path_, file_prefix + ".log"), "w")
    main(args)
