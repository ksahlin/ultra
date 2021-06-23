
from __future__ import print_function
import os,sys
import argparse
import re
import errno
import itertools

import pickle
import math

from collections import defaultdict

import parasail
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


def decide_primary_locations(sam_file, args): # maybe this function is not needed if only one primary alignment from minimap2
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    reads_primary = {}
    reads_multiple_primary = set()
    # reads_tmp = {}

    for read in SAM_file.fetch(until_eof=True):
        if read.flag == 0 or read.flag == 16:
            # ins = sum([length for type_, length in read.cigartuples if type_ == 1])
            # del_ = sum([length for type_, length in read.cigartuples if type_ == 2 and length < args.min_intron ])
            # subs = sum([length for type_, length in read.cigartuples if type_ == 8])
            # matches = sum([length for type_, length in read.cigartuples if type_ == 7])
            # tot_align = ins + del_ + subs + matches
            # identity = matches/float(tot_align)

            if read.query_name in reads_primary:
                reads_multiple_primary.add(read.query_name)

                # if identity >= reads_tmp[read.query_name][0] and  matches >= reads_tmp[read.query_name][1]:
                #     reads_primary[read.query_name] = read
                #     reads_tmp[read.query_name] = (identity, matches)
                # elif identity <= reads_tmp[read.query_name][0] and  matches <= reads_tmp[read.query_name][1]:
                #     continue
                # else:
                #     if identity * matches > reads_tmp[read.query_name][0] * reads_tmp[read.query_name][1]:
                #         reads_primary[read.query_name] = read
                #         reads_tmp[read.query_name] = (identity, matches)
                #     else: 
                #         continue


            else:
                reads_primary[read.query_name] = read
                # reads_tmp[read.query_name] = (identity, matches)
    print("TOTAL READS FLAGGED WITH MULTIPLE PRIMARY:", len(reads_multiple_primary))
    return reads_primary     


def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)


def print_detailed_values_to_file(error_rates, alignment_results, reads, outfile, alignment_algorithm):
    for acc in alignment_results:
        aln_class, nr_exons_true, nr_exons_inferred, max_diff = alignment_results[acc]
        err_rate = error_rates[acc]
        # read = read_alignments[acc]
        read_length = len(reads[acc])
        info_tuple = (acc, alignment_algorithm, err_rate, read_length, aln_class, nr_exons_true, nr_exons_inferred, max_diff) 
        outfile.write( ",".join( [str(item) for item in info_tuple] ) + "\n")


def get_exon_sites(cigar_tuples, first_exon_start, last_exon_end, annotated_chr_coordinate_pairs):
    ref_pos = first_exon_start
    exon_sites = [ref_pos]
    
    for i, (l,t) in enumerate(cigar_tuples):
        if t == "D":
                if (ref_pos, ref_pos + l) in annotated_chr_coordinate_pairs:
                    exon_sites.append( ref_pos )
                    exon_sites.append( ref_pos + l )
                    # print("HEERE")
                ref_pos += l

        elif t == "=" or t== "M" or t == "X":
            ref_pos += l
        elif t == "N":
            splice_start = ref_pos
            ref_pos += l
            splice_stop = ref_pos
            exon_sites.append( splice_start)
            exon_sites.append( splice_stop )

        elif t == "I" or t == "S" or t == "H": # insertion or softclip
            ref_pos += 0

        else: # reference skip or soft/hardclip "~", or match =
            print("UNEXPECTED!", t)
            sys.exit()

    exon_sites.append(last_exon_end)
    exon_sites = [(exon_sites[i],exon_sites[i+1]) for i in range(0, len(exon_sites), 2) ]  
    return exon_sites

def get_read_alignment_exon_sites(reads_primary_locations, annotated_splice_coordinates_pairs):
    read_exon_sites = {}
    for acc in reads_primary_locations:
        read = reads_primary_locations[acc]
        assert read.flag == 0 or read.flag == 16
        q_cigar = read.cigarstring
        q_start = read.reference_start
        q_end = read.reference_end
        read_cigar_tuples = []
        result = re.split(r'[=DXSMINH]+', q_cigar)
        i = 0
        for length in result[:-1]:
            i += len(length)
            type_ = q_cigar[i]
            i += 1
            read_cigar_tuples.append((int(length), type_ ))  

        mod_ref = modify_ref_header_for_alignment(read.reference_name) # convert 1,..,22,X,Y,MT to chr1, .., chr22, chrX, chrY, chrM
        if mod_ref in annotated_splice_coordinates_pairs:
            annotated_chr_coordinate_pairs = annotated_splice_coordinates_pairs[mod_ref]
        else:
            annotated_chr_coordinate_pairs = set()

        read_exon_sites[read.query_name] = {}  
        read_exon_sites[read.query_name][mod_ref] = get_exon_sites(read_cigar_tuples, q_start, q_end, annotated_chr_coordinate_pairs)

    return read_exon_sites


def get_annotated_splicesites(ref_gff_file, infer_genes, outfolder):
    db_name = os.path.join(outfolder, 'database.db')
    if infer_genes:
        fn = gffutils.example_filename(ref_gff_file)
        db = gffutils.create_db(fn, dbfn=db_name, force=True, keep_order=True, merge_strategy='merge', 
                                sort_attribute_values=True)
        db = gffutils.FeatureDB(db_name, keep_order=True)
    else:
        fn = gffutils.example_filename(ref_gff_file)
        db = gffutils.create_db(fn, dbfn=db_name, force=True, keep_order=True, merge_strategy='merge', 
                                sort_attribute_values=True, disable_infer_genes=True, disable_infer_transcripts=True)
        db = gffutils.FeatureDB(db_name, keep_order=True)

    splice_coordinates = {} # to calc individual fraction of correct sites and NIC
    splice_coordinates_pairs = {} 
    ref_isoforms = {} # To calculate Full splice matches
    minimum_annotated_intron = 1000000000
    for gene in db.features_of_type('gene'):
        chromosome = str(gene.seqid)
        if chromosome not in ref_isoforms:
            ref_isoforms[chromosome] = {}
        if chromosome not in splice_coordinates:
            splice_coordinates[chromosome] = set()
            splice_coordinates_pairs[chromosome] = set()

        #add annotated transcripts
        for transcript in db.children(gene, featuretype='transcript', order_by='start'):
            # print(dir(transcript))
            consecutive_exons = [exon for exon in db.children(transcript, featuretype='exon', order_by='start')]
            # print('transcript', transcript.id, transcript.start, transcript.stop, [ (exon.start, exon.stop) for exon in db.children(transcript, featuretype='exon', order_by='start')])
            tmp_splice_sites = []
            for e1,e2 in zip(consecutive_exons[:-1], consecutive_exons[1:]):
                tmp_splice_sites.append( (e1.stop, e2.start -1 ))           
                splice_coordinates[chromosome].add(e1.stop)
                splice_coordinates[chromosome].add(e2.start -1 )
                splice_coordinates_pairs[chromosome].add( (e1.stop, e2.start -1 ) )

                if e2.start -1 - e1.stop < minimum_annotated_intron:
                    minimum_annotated_intron = e2.start -1 - e1.stop
                # print('exon', exon.id, exon.start, exon.stop)
            
            ref_isoforms[chromosome][tuple(tmp_splice_sites)] = transcript.id

    return ref_isoforms, splice_coordinates, splice_coordinates_pairs, minimum_annotated_intron


# def mod_long_acc(acc):
#     if len(acc) > 250:
#         acc = acc.split()[0]
#         if len(acc) > 250: # simulated data
#             #RBBP7|ENSG00000102054|ENST00000380087|16852551;16870038;16849244;16862955;16869076;16857594;16853682;16852046;16852749;16858676;16845828;16844341|16852628;16870414;16849301;16863100;16869220;16857709;16853842;16852122;16852875;16858849;16845938;16845103_4_0.08101157308186883
#             print("here:", acc)
#             tmp1 = acc.split("|")
#             tmp2 = acc.split("_")
#             new_acc = tmp1[:3] +  tmp2[-2:]
#             acc = "|".join(new_acc)
#             print(acc)

#     return acc

def get_alignment_classifications(true_exon_sites, aligned_exon_sites, correct_dist_threshold, is_ultra = False):
    read_annotations = {}
    total_count_exon_sizes = defaultdict(int)
    correct_count_exon_sizes = defaultdict(int)

    for acc in true_exon_sites:
        exons_true = true_exon_sites[acc]

        # if is_ultra:
        #     acc_mod = mod_long_acc(acc)
        # else:
        #     acc_mod = acc

        if acc not in aligned_exon_sites:
            read_annotations[acc] = ('unaligned', len(exons_true), 0, 0)
        else:
            chr_, exons_aligned = list(aligned_exon_sites[acc].items())[0]
            # exons_aligned = list(aln_dict.values())[0]
            aln_start, aln_stop = exons_aligned[0][0], exons_aligned[-1][1]
            true_start, true_stop = exons_true[0][0], exons_true[-1][1]
            if aln_start > true_stop or  aln_stop < true_start:
                read_annotations[acc] = ("diff_location", len(exons_true), len(exons_aligned), 0)

            elif len(exons_true) != len(exons_aligned):
                read_annotations[acc] = ("diff_exon_count", len(exons_true), len(exons_aligned), 0)
            else:
                max_diff = 0
                for i in range(len(exons_true)):
                    (true_start, true_stop) = exons_true[i]
                    (aln_start, aln_stop) = exons_aligned[i]
                    print(true_start - aln_start, true_stop - aln_stop)
                    # tmp_diff = max(math.fabs(true_start - (aln_start + 1)), math.fabs(true_stop - aln_stop) ) # (aln_start + 1) to make 1-based coordinate because biomart is 1 based (for annotated datasets ENS and SIM_ANN)
                    tmp_diff = max(math.fabs(true_start - aln_start), math.fabs(true_stop - aln_stop) ) # For SIM_NIC dataset which we simulate ourselves
                    
                    if tmp_diff <= correct_dist_threshold:
                        correct_count_exon_sizes[true_stop-true_start + 1] += 1

                    if tmp_diff > max_diff:
                        max_diff = tmp_diff

                if max_diff > correct_dist_threshold:
                    is_correct = False
                else:
                    is_correct = True

                if is_correct:
                    read_annotations[acc] = ("correct", len(exons_true), len(exons_aligned), max_diff)
                else:
                    read_annotations[acc] = ("site_diff", len(exons_true), len(exons_aligned), max_diff)

        for (true_start, true_stop) in exons_true:
            total_count_exon_sizes[true_stop-true_start + 1] += 1

    return read_annotations, total_count_exon_sizes, correct_count_exon_sizes


def pickle_dump(data, filename):
    with open(os.path.join(args.outfolder,filename), 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)


def pickle_load(filename):
    with open(filename, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        data = pickle.load(f)
    return data


def get_error_rates(reads):
    ## ENSG00000105143|ENST00000598504|14961310|14979693_0_0.09273518580515568
    error_rates = {}
    for acc in reads:
        err_rate = acc.split("_")[-1]
        error_rates[acc] = err_rate
    return error_rates

def modify_ref_header_for_alignment(header):
    if header.isdigit() or header == 'X' or header == 'Y':
        return 'chr'+ header
    elif header == 'MT':
        return 'chrM'
    else:
        return header


def modify_reference_headers(refs):
    modified = False
    for header in list(refs.keys()):
        if header.isdigit() or header == 'X' or header == 'Y':
            chr_id = 'chr'+ header
        elif header == 'MT':
            chr_id = 'chrM'
        else:
            chr_id = header

        # we have modified 
        if chr_id != header:
            modified = True
            seq = refs[header]
            del refs[header]
            refs[chr_id] = seq
    return modified

def get_true_exon_sites(accessions_map):
    true_exon_sites = {}
    for line in open(accessions_map, 'r'):
        acc, full_acc = line.split(',')
        if "_" in acc.split("|")[0]:
            exon_starts = full_acc.split("|")[3].split(";")
            tmp_acc2 = full_acc.split("|")[4]
            exon_stops = tmp_acc2.split("_")[0].split(";")
        else:
            tmp_acc = full_acc.split("_")[0]
            exon_starts, exon_stops = tmp_acc.split("|")[3].split(";"), tmp_acc.split("|")[4].split(";") 
        # print(acc)
        # print(exon_starts, exon_stops)

        e = sorted([int(pos) for pos in  exon_starts + exon_stops])
        true_exon_sites[acc] = [(e[i],e[i+1]) for i in range(0, len(e), 2) ]  
    return true_exon_sites
    # for acc in reads:
    #     e = sorted([int(pos) for pos in  acc.split("|")[3].split(";")])
    #     true_exon_sites[acc] = [(e[i],e[i+1]) for i in range(0, len(e), 2) ]  

def print_correctness_per_exon_size(correctness_per_exon_size_outfile, total_count_exon_sizes, correct_count_exon_sizes, alignment_method):
    for e_size in sorted(list(total_count_exon_sizes.keys())):
        nr_total = total_count_exon_sizes[e_size]
        if e_size in correct_count_exon_sizes:
            nr_corr = correct_count_exon_sizes[e_size]
        else:
            nr_corr = 0
        correctness_per_exon_size_outfile.write("{0},{1},{2},{3},{4}\n".format(e_size, nr_total, nr_corr, float(nr_corr)/nr_total, alignment_method))


def main(args):

    if args.load_database:
        print()
        print("LOADING FROM PRECOMPUTED DATABASE")
        print()
        annotated_ref_isoforms = pickle_load(os.path.join( args.outfolder, 'annotated_ref_isoforms.pickle') )
        annotated_splice_coordinates = pickle_load(os.path.join( args.outfolder, 'annotated_splice_coordinates.pickle') )
        annotated_splice_coordinates_pairs = pickle_load(os.path.join( args.outfolder, 'annotated_splice_coordinates_pairs.pickle') )
        minimum_annotated_intron = pickle_load(os.path.join( args.outfolder, 'minimum_annotated_intron.pickle') )
    else:
        annotated_ref_isoforms, annotated_splice_coordinates, annotated_splice_coordinates_pairs, minimum_annotated_intron = get_annotated_splicesites(args.gff_file, args.infer_genes, args.outfolder)
        pickle_dump(annotated_ref_isoforms, os.path.join( args.outfolder, 'annotated_ref_isoforms.pickle') )
        pickle_dump(annotated_splice_coordinates, os.path.join( args.outfolder, 'annotated_splice_coordinates.pickle') )
        pickle_dump(annotated_splice_coordinates_pairs, os.path.join( args.outfolder, 'annotated_splice_coordinates_pairs.pickle') )
        pickle_dump(minimum_annotated_intron, os.path.join( args.outfolder, 'minimum_annotated_intron.pickle') )


    reads = { acc.split()[0] : seq for i, (acc, (seq, qual)) in enumerate(readfq(open(args.reads, 'r')))}
    print("Total reads", len(reads))

    error_rates = get_error_rates(reads)
    true_exon_sites = get_true_exon_sites(args.accessions_map)
    refs = { acc.split()[0] : seq for i, (acc, (seq, _)) in enumerate(readfq(open(args.refs, 'r')))}
    modify_reference_headers(refs)
    detailed_results_outfile = open(os.path.join(args.outfolder, "results_per_read.csv"), "w")
    detailed_results_outfile.write("acc,alignment_algorithm,error_rate,read_length,alignment_classification,nr_exons_true,nr_exons_inferred,max_diff\n")
    #acc, alignment_algorithm, err_rate, read_length, aln_class, nr_exons_true, nr_exons_inferred, max_diff
    correctness_per_exon_size_outfile = open(os.path.join(args.outfolder, "correctness_per_exon_size.csv"), "w")
    correctness_per_exon_size_outfile.write("exon_size,nr_total,nr_corr,fraction_correct,alignment_algorithm\n")

    if args.ultra_sam:
        torkel_primary_locations = decide_primary_locations(args.ultra_sam, args)
        torkel_exon_sites = get_read_alignment_exon_sites(torkel_primary_locations, annotated_splice_coordinates_pairs)
        print('uLTRA')
        torkel_alignment_results, total_count_exon_sizes, correct_count_exon_sizes = get_alignment_classifications(true_exon_sites, torkel_exon_sites, args.correct_diff)
        print_correctness_per_exon_size(correctness_per_exon_size_outfile, total_count_exon_sizes, correct_count_exon_sizes, "uLTRA")
        reads_unaligned_in_torkel = set(reads.keys()) - set(torkel_primary_locations.keys())
        print_detailed_values_to_file(error_rates, torkel_alignment_results, reads, detailed_results_outfile, "uLTRA")
        print("Reads successfully aligned uLTRA:", len(torkel_primary_locations))
        print("READS UNALIGNED uLTRA:", len(reads_unaligned_in_torkel) )

    if args.ultra_mm2_sam:
        ultra_mm2_primary_locations = decide_primary_locations(args.ultra_mm2_sam, args)
        ultra_mm2_exon_sites = get_read_alignment_exon_sites(ultra_mm2_primary_locations, annotated_splice_coordinates_pairs)
        print('uLTRA')
        torkel_alignment_results, total_count_exon_sizes, correct_count_exon_sizes = get_alignment_classifications(true_exon_sites, ultra_mm2_exon_sites, args.correct_diff)
        print_correctness_per_exon_size(correctness_per_exon_size_outfile, total_count_exon_sizes, correct_count_exon_sizes, "uLTRA_mm2")
        reads_unaligned_in_ultra_mm2 = set(reads.keys()) - set(ultra_mm2_primary_locations.keys())
        print_detailed_values_to_file(error_rates, torkel_alignment_results, reads, detailed_results_outfile, "uLTRA_mm2")
        print("Reads successfully aligned uLTRA:", len(ultra_mm2_primary_locations))
        print("READS UNALIGNED uLTRA:", len(reads_unaligned_in_ultra_mm2) )

    if args.mm2_sam:
        mm2_primary_locations = decide_primary_locations(args.mm2_sam, args)
        mm2_exon_sites = get_read_alignment_exon_sites(mm2_primary_locations, annotated_splice_coordinates_pairs)
        print('MINIMAP2')
        mm2_alignment_results, total_count_exon_sizes, correct_count_exon_sizes = get_alignment_classifications(true_exon_sites, mm2_exon_sites, args.correct_diff)
        print_correctness_per_exon_size(correctness_per_exon_size_outfile, total_count_exon_sizes, correct_count_exon_sizes, "minimap2")
        reads_unaligned_in_mm2 = set(reads.keys()) - set(mm2_primary_locations.keys()) 
        print_detailed_values_to_file(error_rates, mm2_alignment_results, reads, detailed_results_outfile, "minimap2")    
        print("Reads successfully aligned mm2:", len(mm2_primary_locations))
        print("READS UNALIGNED mm2:", len(reads_unaligned_in_mm2) )

    if args.mm2_gtf_sam:
        mm2_gtf_primary_locations = decide_primary_locations(args.mm2_gtf_sam, args)
        mm2_gtf_exon_sites = get_read_alignment_exon_sites(mm2_gtf_primary_locations, annotated_splice_coordinates_pairs)
        print('MINIMAP2')
        mm2_gtf_alignment_results, total_count_exon_sizes, correct_count_exon_sizes = get_alignment_classifications(true_exon_sites, mm2_gtf_exon_sites, args.correct_diff)
        print_correctness_per_exon_size(correctness_per_exon_size_outfile, total_count_exon_sizes, correct_count_exon_sizes, "minimap2_GTF")
        reads_unaligned_in_mm2_gtf = set(reads.keys()) - set(mm2_gtf_primary_locations.keys()) 
        print_detailed_values_to_file(error_rates, mm2_gtf_alignment_results, reads, detailed_results_outfile, "minimap2_GTF")    
        print("Reads successfully aligned mm2:", len(mm2_gtf_primary_locations))
        print("READS UNALIGNED mm2:", len(reads_unaligned_in_mm2_gtf) )

    if args.graphmap2_sam:
        graphmap2_primary_locations = decide_primary_locations(args.graphmap2_sam, args)
        graphmap2_exon_sites = get_read_alignment_exon_sites(graphmap2_primary_locations, annotated_splice_coordinates_pairs)
        print("GraphMap2")
        graphmap2_alignment_results, total_count_exon_sizes, correct_count_exon_sizes = get_alignment_classifications(true_exon_sites, graphmap2_exon_sites, args.correct_diff)
        print_correctness_per_exon_size(correctness_per_exon_size_outfile, total_count_exon_sizes, correct_count_exon_sizes, "Graphmap2")
        reads_unaligned_in_graphmap2 = set(reads.keys()) - set(graphmap2_primary_locations.keys()) 
        print_detailed_values_to_file(error_rates, graphmap2_alignment_results, reads, detailed_results_outfile, "Graphmap2")
        print("Reads successfully aligned graphmap2:", len(graphmap2_primary_locations))
        print("READS UNALIGNED graphmap2:", len(reads_unaligned_in_graphmap2) )

    if args.graphmap2_gtf_sam:
        graphmap2_gtf_primary_locations = decide_primary_locations(args.graphmap2_gtf_sam, args)
        graphmap2_gtf_exon_sites = get_read_alignment_exon_sites(graphmap2_gtf_primary_locations, annotated_splice_coordinates_pairs)
        print("GraphMap2_GTF")
        graphmap2_gtf_alignment_results, total_count_exon_sizes, correct_count_exon_sizes = get_alignment_classifications(true_exon_sites, graphmap2_gtf_exon_sites, args.correct_diff)
        print_correctness_per_exon_size(correctness_per_exon_size_outfile, total_count_exon_sizes, correct_count_exon_sizes, "Graphmap2_GTF")
        reads_unaligned_in_graphmap2_gtf = set(reads.keys()) - set(graphmap2_gtf_primary_locations.keys()) 
        print_detailed_values_to_file(error_rates, graphmap2_gtf_alignment_results, reads, detailed_results_outfile, "Graphmap2_GTF")
        print("Reads successfully aligned graphmap2:", len(graphmap2_gtf_primary_locations))
        print("READS UNALIGNED graphmap2:", len(reads_unaligned_in_graphmap2_gtf) )

    if args.desalt_sam:
        desalt_primary_locations = decide_primary_locations(args.desalt_sam, args)
        desalt_exon_sites = get_read_alignment_exon_sites(desalt_primary_locations, annotated_splice_coordinates_pairs)
        print("deSALT")
        desalt_alignment_results, total_count_exon_sizes, correct_count_exon_sizes = get_alignment_classifications(true_exon_sites, desalt_exon_sites, args.correct_diff)
        print_correctness_per_exon_size(correctness_per_exon_size_outfile, total_count_exon_sizes, correct_count_exon_sizes, "deSALT")
        reads_unaligned_in_desalt = set(reads.keys()) - set(desalt_primary_locations.keys()) 
        print_detailed_values_to_file(error_rates, desalt_alignment_results, reads, detailed_results_outfile, "deSALT")
        print("Reads successfully aligned deSALT:", len(desalt_primary_locations))
        print("READS UNALIGNED deSALT:", len(reads_unaligned_in_desalt) )

    if args.desalt_gtf_sam:
        desalt_gtf_primary_locations = decide_primary_locations(args.desalt_gtf_sam, args)
        desalt_gtf_exon_sites = get_read_alignment_exon_sites(desalt_gtf_primary_locations, annotated_splice_coordinates_pairs)
        print("deSALT_gtf")
        desalt_gtf_alignment_results, total_count_exon_sizes, correct_count_exon_sizes = get_alignment_classifications(true_exon_sites, desalt_gtf_exon_sites, args.correct_diff)
        print_correctness_per_exon_size(correctness_per_exon_size_outfile, total_count_exon_sizes, correct_count_exon_sizes, "deSALT_GTF")
        reads_unaligned_in_desalt_gtf = set(reads.keys()) - set(desalt_gtf_primary_locations.keys()) 
        print_detailed_values_to_file(error_rates, desalt_gtf_alignment_results, reads, detailed_results_outfile, "deSALT_GTF")
        print("Reads successfully aligned deSALT_gtf:", len(desalt_gtf_primary_locations))
        print("READS UNALIGNED deSALT_gtf:", len(reads_unaligned_in_desalt_gtf) )

    detailed_results_outfile.close()
    correctness_per_exon_size_outfile.close()



    ###########################################################################
    ###########################################################################




    # print(torkel_sorted)
    # print(",".join([s for s in torkel_stats]))
    # print(",".join([s for s in mm2_stats]))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('--ultra_sam', type=str, default = '', help='Path to the original read file')
    parser.add_argument('--ultra_mm2_sam', type=str, default = '', help='Path to the original read file')
    parser.add_argument('--mm2_sam', type=str, default = '', help='Path to the corrected read file')
    parser.add_argument('--mm2_gtf_sam', type=str, default = '', help='Path to the corrected read file')
    parser.add_argument('--desalt_sam', type=str, default = '', help='Path to the corrected read file')
    parser.add_argument('--desalt_gtf_sam', type=str, default = '', help='Path to the corrected read file')
    parser.add_argument('--graphmap2_sam', type=str, default = '', help='Path to the corrected read file')
    parser.add_argument('--graphmap2_gtf_sam', type=str, default = '', help='Path to the corrected read file')
    parser.add_argument('reads', type=str, help='Path to the read file')
    parser.add_argument('refs', type=str, help='Path to the refs file')
    parser.add_argument('gff_file', type=str, help='Path to the refs file')
    parser.add_argument('accessions_map', type=str, help='Path to accessions mapping file')
    parser.add_argument('outfolder', type=str, help='Output path of results')
    parser.add_argument('--min_intron', type=int, default=15, help='Threchold for what is counted as varation/intron in alignment as opposed to deletion.')
    parser.add_argument('--infer_genes', action= "store_true", help='Include pairwise alignment of original and corrected read.')
    parser.add_argument('--load_database', action= "store_true", help='Load already computed splice junctions and transcript annotations instead of constructing a new database.')
    parser.add_argument('--simulated', action= "store_true", help='Adds extra analysis that can be done for simulated data since known true locations.')
    parser.add_argument('--correct_diff', type=int, default=1, help='Adds extra analysis that can be done for simulated data since known true locations.')

    # parser.add_argument('--align', action= "store_true", help='Include pairwise alignment of original and corrected read.')

    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)

