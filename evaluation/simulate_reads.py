import os
import argparse
import numpy as np
import misc_functions
import random
import math
from collections import defaultdict

# log how many reads are simulated for each variant.


def simulate_read(i, transcript_acc, isoform ):
    error_lvls = [0.7, 0.8, 0.9, 0.92, 0.94, 0.96,0.98, 0.99, 0.995]
    read = []
    qual = []
    was_del = False
    del_length = 0
    prev_base_pos = 0
    for l, n in enumerate(isoform):
        # if l <= 15 or l >= len(isoform) - 15: # no errors first and last 15 bases
        #     p_correct_reading = 1.0
        #     p_error = 1.0 - 0.995
        # else:
        if l < prev_base_pos + del_length:
            print('skipping')
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
            r = random.uniform(0,1)
            if r < 0.6: #deletion
                was_del = p_error
                del_length += 1
                r_del = random.uniform(0,1)
                while r_del >= 0.7:
                    del_length += 1
                    r_del = random.uniform(0,1)
                    qual.append( round(-math.log(0.7,10)*10) )
            print(del_length)
            elif 0.6 <= r < 0.9:
                read.append(random.choice("ACGT"))
                qual.append( round(-math.log(p_error,10)*10) )
                prev_base_pos = l

            else:
                read.append(n)
                qual.append( round(-math.log(p_error,10)*10) )

                r_ins = random.uniform(0,1)
                while r_ins >= 0.7:
                    read.append(random.choice("ACGT"))
                    r_ins = random.uniform(0,1)
                    qual.append( round(-math.log(0.7,10)*10) )
                prev_base_pos = l

        else:
            if was_del: # adding uncertainty from prevous deleted base
                read.append(n)
                qual.append( round(-math.log(was_del,10)*10) )
            else:
                read.append(n)
                qual.append( round(-math.log(p_error,10)*10) )
            prev_base_pos = l

            was_del = False
            del_length = 1

    read_seq = "".join([n for n in read])
    qual_seq = "".join([chr(q + 33) for q in qual])
    acc = str(transcript_acc) + "_" +  str(i) #= (read_seq, qual_seq)

    # print(read_seq)
    # print(qual_seq)

    return acc, read_seq, qual_seq


def main(args):
    sequence_transcripts = {}
    sequence_transcripts = dict(misc_functions.read_fasta(open(args.sequence_material,"r")))

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
    for acc, (read_seq,qual_seq) in sorted(ont_reads.items(), key = lambda x: len(x[1]), reverse = True):
        outfile.write("@{0}\n{1}\n{2}\n{3}\n".format(acc, read_seq, "+", qual_seq))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate pacbio reads from a set of transcripts.")
    parser.add_argument('sequence_material', type=str, help='The fasta file with sequences to be sequenced.')
    parser.add_argument('outfile', type=str, help='Output path to fasta file')
    parser.add_argument('read_count', type=int, help='Number of reads to simulate.')
    # parser.add_argument('config', type=str, help='config file')


    args = parser.parse_args()
    path_, file_prefix = os.path.split(args.outfile)
    misc_functions.mkdir_p(path_)
    args.logfile = open(os.path.join(path_, file_prefix + ".log"), "w")
    main(args)
