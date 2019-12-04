import os
import argparse
import numpy as np
import random
import math
from collections import defaultdict
import errno

import gffutils
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
    err_rate = float(err)/ len(isoform)
    

    #Example: RBBP7|ENSG00000102054|ENST00000380087|16852551;16870038;16849244;16862955;16869076;16857594;16853682;16852046;16852749;16858676;16845828;16844341|16852628;16870414;16849301;16863100;16869220;16857709;16853842;16852122;16852875;16858849;16845938;16845103_4_0.08101157308186883
    tmp1 = transcript_acc.split("|")
    shortened_read_acc = "|".join(tmp1[:3]) + "_" +  str(i) + "_" + str(err_rate)
    # transcript_acc = "|".join(new_acc)
    # shortened_read_acc = transcript_acc
    # acc = str(transcript_acc) + "_" +  str(i) + "_" + str(err_rate) #= (read_seq, qual_seq)
    full_read_acc = str(transcript_acc) + "_" +  str(i) + "_" + str(err_rate)
    print(err, len(isoform), float(err)/ len(isoform))
    # print(read_seq)
    # print(qual_seq)

    return shortened_read_acc, full_read_acc, read_seq, qual_seq


def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

from itertools import chain, combinations
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def generate_nics(db, sequence_material):
    refs = {acc : seq for acc, (seq, _) in readfq(open(args.sequence_material,"r"))}

    nic_transcripts = {}
    for gene in db.features_of_type('gene'):
        # genes_to_ref[gene.id] = str(gene.seqid)
        print("here", gene.id,  str(gene.seqid))
        chr_id = gene.seqid
        annotated = set()
        nr_transcripts = 0
        for transcript in db.children(gene, featuretype='transcript', order_by='start'):  
            exons = tuple((exon.seqid, exon.start - 1, exon.stop) for exon in db.children(transcript, featuretype='exon', order_by='start'))
            annotated.add(exons)
            # print(exons)
            nr_transcripts += 1

        gene_exons = [(exon.seqid, exon.start - 1, exon.stop) for exon in db.children(gene, featuretype='exon', order_by='start')]

        # randomly select internal exons with p=0.5 and check whether this is already in annotated, if not add to nic
        print( len(annotated),len(gene_exons))

        if len(gene_exons) > 3:
            nr_fails = 0
            nr_nic = 0
            while nr_nic < len(annotated):
                new_internal_exons = [ e for e in gene_exons[1:-1] if random.uniform(0, 1) > 0.5]
                candidate_nic = tuple([gene_exons[0]] + new_internal_exons +  [gene_exons[-1]])
                print( "l", len(candidate_nic))
                if candidate_nic in annotated:
                    nr_fails +=1
                else:
                    # print(candidate_nic)
                    nic_id = "{0}|{1}|{2}|{3}|{4}".format(str(gene.id), str(gene.seqid), nr_nic, ";".join([ str(start) for (s_id, start, stop) in candidate_nic ]), ";".join([str(stop) for (s_id, start, stop) in candidate_nic ]) )
                    exons_seqs = []
                    for s_id, start,stop in candidate_nic: 
                        seq = refs[chr_id][start : stop] 
                        exons_seqs.append(seq)
                    nic_seq = "".join([s for s in exons_seqs])
                    nic_transcripts[nic_id] = nic_seq
                    nr_nic += 1
                    nr_fails = 0
                if nr_fails > 5:
                    break
    print(len(nic_transcripts))

    return nic_transcripts

def main(args):
    if args.nic:
        database = os.path.join(args.outfolder,'database.db')
        if os.path.isfile(database):
            print("Database found in directory using this one.")
            print("If you want to recreate the database, please remove the file: {0}".format(database))
            print()
            db = gffutils.FeatureDB(database, keep_order=True)
        elif not args.disable_infer:
            fn = gffutils.example_filename(args.gtf)
            db = gffutils.create_db(fn, dbfn=database, force=True, keep_order=True, merge_strategy='merge', 
                                    sort_attribute_values=True)
            db = gffutils.FeatureDB(database, keep_order=True)
        else:
            fn = gffutils.example_filename(args.gtf)
            db = gffutils.create_db(fn, dbfn=database, force=True, keep_order=True, merge_strategy='merge', 
                                    sort_attribute_values=True, disable_infer_genes=True, disable_infer_transcripts=True)
            db = gffutils.FeatureDB(database, keep_order=True)

        sequence_transcripts = generate_nics(db, args.sequence_material)
        # sequence_transcripts = {}

    else:
        sequence_transcripts = {seq : acc for acc, (seq, _) in readfq(open(args.sequence_material,"r")) }
        print(len(sequence_transcripts))
        sequence_transcripts = {acc: seq for seq, acc in sequence_transcripts.items() }

    # just generate all numbers at once and draw from this 5x should be enough
    ont_reads = {}
    reads_generated_log = defaultdict(int)
    errors = []
    all_transctript_accessions = list(sequence_transcripts.keys())
    for i in range(args.read_count):
        acc = random.choice( all_transctript_accessions)
        transcript = sequence_transcripts[acc]
        read_acc, full_read_acc,  read, qual = simulate_read(i, acc, transcript)
        ont_reads[read_acc] = (read, qual, full_read_acc)
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

    outfile_fasta = open(args.outfile_prefix + ".fa", "w")
    accsessions_outfile = open(args.full_acc_file, "w")
    # if args.fasta:
    for read_acc, (read_seq,qual_seq, full_read_acc) in sorted(ont_reads.items(), key = lambda x: len(x[1]), reverse = True):
        outfile_fasta.write(">{0}\n{1}\n".format(read_acc, read_seq))
        accsessions_outfile.write("{0},{1}\n".format(read_acc, full_read_acc))
    # else:
    outfile_fastq = open(args.outfile_prefix + ".fq", "w")
    for read_acc, (read_seq,qual_seq, full_read_acc) in sorted(ont_reads.items(), key = lambda x: len(x[1]), reverse = True):
        outfile_fastq.write("@{0}\n{1}\n{2}\n{3}\n".format(read_acc, read_seq, "+", qual_seq))
        accsessions_outfile.write("{0},{1}\n".format(read_acc, full_read_acc))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate pacbio reads from a set of transcripts.")
    parser.add_argument('sequence_material', type=str, help='The fasta file with sequences to be sequenced.')
    parser.add_argument('outfile_prefix', type=str, help='Output path to fasta file')
    parser.add_argument('read_count', type=int, help='Number of reads to simulate.')
    # parser.add_argument('--fasta', action="store_true", help='Output in fasta format')
    parser.add_argument('--nic', action="store_true", help='Simulate NIC transcripts')
    parser.add_argument('--gtf', type=str, default = '', help='GTF to simulate NIC from.')
    parser.add_argument('--disable_infer', action="store_true", help='GTF to simulate NIC from.')
    # parser.add_argument('config', type=str, help='config file')


    args = parser.parse_args()
    path_, file_prefix = os.path.split(args.outfile_prefix)
    mkdir_p(path_)
    args.logfile = open(os.path.join(path_, file_prefix + ".log"), "w")
    args.outfolder = path_
    args.full_acc_file = os.path.join(path_, "accessions_map.csv")
    main(args)
