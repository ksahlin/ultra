
from __future__ import print_function
import os,sys
import argparse
import re
import errno
import itertools

import pickle

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
    reads_tmp = {}

    for read in SAM_file.fetch(until_eof=True):
        if read.flag == 0 or read.flag == 16:
            ins = sum([length for type_, length in read.cigartuples if type_ == 1])
            del_ = sum([length for type_, length in read.cigartuples if type_ == 2 and length < args.min_intron ])
            subs = sum([length for type_, length in read.cigartuples if type_ == 8])
            matches = sum([length for type_, length in read.cigartuples if type_ == 7])
            tot_align = ins + del_ + subs + matches
            identity = matches/float(tot_align)

            # has_large_del = [length for type_, length in read.cigartuples if type_ == 2 and length >=  args.min_intron ]
            # if has_large_del:
            #     print(has_large_del)
            if read.query_name in reads_primary:
                reads_multiple_primary.add(read.query_name)
                if identity >= reads_tmp[read.query_name][0] and  matches >= reads_tmp[read.query_name][1]:
                    reads_primary[read.query_name] = read
                    reads_tmp[read.query_name] = (identity, matches)
                elif identity <= reads_tmp[read.query_name][0] and  matches <= reads_tmp[read.query_name][1]:
                    continue
                else:
                    if identity * matches > reads_tmp[read.query_name][0] * reads_tmp[read.query_name][1]:
                        reads_primary[read.query_name] = read
                        reads_tmp[read.query_name] = (identity, matches)
                    else: 
                        continue

                        
                    # print( "Ambiguous, prefferred:", round(reads_tmp[read.query_name][0],2), reads_tmp[read.query_name][1], "over", round(identity,2), matches)

                #     if tot_align >= reads_tmp[read.query_name][1]:
                #         reads_primary[read.query_name] = read
                #     else:
                #         continue
                # elif tot_align >= reads_tmp[read.query_name][1]:
                #     if identity >= reads_tmp[read.query_name][0]:
                #         reads_primary[read.query_name] = read
                #     else:
                #         continue

            else:
                reads_primary[read.query_name] = read
                reads_tmp[read.query_name] = (identity, matches)
    print("TOTAL READS FLAGGED WITH MULTIPLE PRIMARY:", len(reads_multiple_primary))
    return reads_primary     


def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)


def print_detailed_values_to_file(error_rates, annotations_dict, reads, outfile, reads_unaligned_in_other_method, read_type, read_alignments):
    for acc in reads:
        if acc in error_rates:
            err_rate = error_rates[acc]
        else:
            err_rate = "-"

        if acc not in read_alignments:
            chr_id = 'unaligned'
            reference_start = '-'  #read.reference_name, read.reference_start, read.reference_end + 1, read.flag,
            reference_end = '-'
            flag = '-'
        else:
            read = reads[acc]
            reference_name, reference_start, reference_end, flag = read.reference_name, read.reference_start, read.reference_end + 1, read.flag

        read_class = annotations_dict[acc] 
        read_length = len(reads[acc])
        is_unaligned_in_other_method = 1 if acc in reads_unaligned_in_other_method else 0
        info_tuple = (acc, read_type, err_rate, read_length, is_unaligned_in_other_method, *read_class, reference_name, reference_start, reference_end, flag) # 'tot_splices', 'read_sm_junctions', 'read_nic_junctions', 'fsm', 'nic', 'ism', 'nnc', 'no_splices'  )
        outfile.write( ",".join( [str(item) for item in info_tuple] ) + "\n")


def get_splice_sites(cigar_tuples, first_exon_start, minimum_annotated_intron, annotated_chr_coordinate_pairs):
    splice_sites = []
    ref_pos = first_exon_start
    
    for i, (l,t) in enumerate(cigar_tuples):
        if t == "D":
            # if l >= minimum_annotated_intron -1:
            #     # print("long del", l)
            #     splice_start = ref_pos
            #     ref_pos += l
            #     splice_stop = ref_pos
            #     splice_sites.append( (splice_start, splice_stop) )
                
            # else:
                if (ref_pos, ref_pos + l) in annotated_chr_coordinate_pairs:
                    splice_sites.append( (ref_pos, ref_pos + l) )
                    print("HEERE")
                ref_pos += l

                # if l > 15:
                #     print("Large deletion!!!", l)


        elif t == "=" or t== "M" or t == "X":
            ref_pos += l
        elif t == "N":
            splice_start = ref_pos
            ref_pos += l
            splice_stop = ref_pos
            splice_sites.append( (splice_start, splice_stop) )

        elif t == "I" or t == "S" or t == "H": # insertion or softclip
            ref_pos += 0

        else: # reference skip or soft/hardclip "~", or match =
            print("UNEXPECTED!", t)
            sys.exit()

    return splice_sites

def get_read_candidate_splice_sites(reads_primary_locations, minimum_annotated_intron, annotated_splice_coordinates_pairs):
    read_splice_sites = {}
    for acc in reads_primary_locations:
        read = reads_primary_locations[acc]
        if read.flag == 0 or read.flag == 16:
            # compare cs tag at intron sites
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

            if read.reference_name in annotated_splice_coordinates_pairs:
                annotated_chr_coordinate_pairs = annotated_splice_coordinates_pairs[read.reference_name]
            else:
                annotated_chr_coordinate_pairs = set()

            read_splice_sites[read.query_name] = {}  
            read_splice_sites[read.query_name][read.reference_name] = get_splice_sites(read_cigar_tuples, q_start, minimum_annotated_intron, annotated_chr_coordinate_pairs)

        

    return read_splice_sites


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


from collections import namedtuple
def get_splice_classifications(annotated_ref_isoforms, annotated_splice_coordinates, annotated_splice_coordinates_pairs,  all_reads_splice_sites, ref_seqs, reads_primary_locations):

    total_reads = 0
    
    #individual splice sites
    total_individual_true = 0
    total_individual_in_data = 0

    #pairs of splice sites
    total_pairs_in_data = 0
    tot_nic_pairs = 0
    tot_sm_pairs = 0

    #whole structure
    total_transcript_fsm = 0
    total_transcript_nic = 0
    total_transcript_ism = 0
    total_transcript_nnc = 0
    total_transcript_no_splices = 0

    read_annotations = {}
    alignments_not_matching_annotated_sites = 0
    # print(annotated_ref_isoforms["SIRV5"])
    # for tr in annotated_ref_isoforms["SIRV5"]:
    #     print(tr)
    canonical_splice = 0
    all_splice = 0
    for read_acc in all_reads_splice_sites:
        read_annotations[read_acc] = {}
        total_reads += 1
        assert len(all_reads_splice_sites[read_acc]) == 1
        chr_id, read_splice_sites = list(all_reads_splice_sites[read_acc].items())[0]

        # print(reads_primary_locations[read_acc].flag, read_splice_sites)
        # print(chr_id)
        # print(annotated_splice_coordinates)
        if chr_id not in annotated_splice_coordinates:
            alignments_not_matching_annotated_sites += 1
            annotated_sites = set()
            annotated_pairs = set()
            annotated_isoforms = set()
        else:
            annotated_sites = annotated_splice_coordinates[chr_id]
            annotated_pairs = annotated_splice_coordinates_pairs[chr_id]
            annotated_isoforms = annotated_ref_isoforms[chr_id]

        # print(annotated_pairs)
        # print( chr_id, read_splice_sites)
        # print(annotated_ref_isoforms[chr_id])
        read_sm_junctions = 0
        read_nic_junctions = 0
        read_splice_letters = []
        read_splice_choords = []
        total_individual_in_data += len(read_splice_sites)*2
        total_pairs_in_data += len(read_splice_sites)

        for splice_site in read_splice_sites:
            start_sp, stop_sp = splice_site
            if reads_primary_locations[read_acc].flag == 0:
                donor = ref_seqs[chr_id][start_sp: start_sp + 2] 
                acceptor = ref_seqs[chr_id][stop_sp - 2: stop_sp]
            else:
                acceptor = reverse_complement(ref_seqs[chr_id][start_sp: start_sp + 2])
                donor = reverse_complement(ref_seqs[chr_id][stop_sp - 2: stop_sp])

            if donor == "GT" and acceptor == "AG":
                canonical_splice += 1
            all_splice += 1

            read_splice_letters.append( donor + str("-") + acceptor )
            read_splice_choords.append( str(start_sp) + str("-") + str(stop_sp) )
            # print(splice_site)
            if (start_sp, stop_sp) in annotated_pairs:
                tot_sm_pairs += 1 
                total_individual_true += 2
                read_sm_junctions += 1

            elif start_sp in annotated_sites and stop_sp in annotated_sites:
                # print((start_sp, stop_sp), annotated_pairs )
                tot_nic_pairs += 1
                total_individual_true += 2
                read_nic_junctions += 1

            elif start_sp in annotated_sites:
                total_individual_true += 1  

            elif stop_sp in annotated_sites:
                total_individual_true += 1  

            # for sp in splice_site:
            #     # print(sp)
            #     # print(annotated_sites)
            #     total += 1
            #     if sp in annotated_sites:
            #         total_true += 1

        # check set intersection between read splice sites and annotated splice sites
        read_nic, read_ism, read_nnc, read_no_splices = 0,0,0,0
        read_type = ''
        transcript_fsm_id = "NA"
        # print(tuple(read_splice_sites))
        # print(annotated_ref_isoforms[chr_id])
        # print(annotated_ref_isoforms)

        if len(read_splice_sites) > 0:

            if tuple(read_splice_sites) in annotated_isoforms:
                total_transcript_fsm += 1
                read_type = 'FSM'
                transcript_fsm_id = annotated_isoforms[ tuple(read_splice_sites) ]


            elif len(read_splice_sites) == read_sm_junctions + read_nic_junctions:
                if read_nic_junctions >= 1:
                    total_transcript_nic += 1
                    read_type = 'NIC'
                    # print('NIC', read_acc)

                else:
                    total_transcript_ism += 1
                    read_type = 'ISM'

                    # print(read_acc)
                    # print(tuple(read_splice_sites))
                    # for an in annotated_ref_isoforms[chr_id]:
                    #     print(an)
                    # print()
                    # print(annotated_ref_isoforms[chr_id])
            else:
                total_transcript_nnc += 1
                read_nnc = 1
                read_type = 'NNC'
        else:
            total_transcript_no_splices += 1                
            read_type = 'NO_SPLICE'

        read_annotation = namedtuple('Annotation', ['tot_splices', 'read_sm_junctions', 'read_nic_junctions', 'annotation', "donor_acceptors", "donor_acceptors_choords", "transcript_fsm_id" ])
        if read_splice_letters:
            donor_acceptors = ":".join([str(item) for item in read_splice_letters])
            donor_acceptors_choords = ":".join([str(item) for item in read_splice_choords])

        else: 
            donor_acceptors = "NA"
            donor_acceptors_choords = "NA"
        read_annotations[read_acc] = read_annotation( len(read_splice_sites), read_sm_junctions, read_nic_junctions, read_type, donor_acceptors, donor_acceptors_choords, transcript_fsm_id )
                # print("FSM!!")
    # print(annotated_ref_isoforms[chr_id])
    # print( tuple(read_splice_sites))
    print("Total splice sizes found in cigar in reads (individual, pairs):", total_individual_in_data, total_pairs_in_data, "total matching annotations (individual):", total_individual_true,
             "total annotated junctions (splice match pairs):", tot_sm_pairs,  "total NIC junctions in (pairs):", tot_nic_pairs, "total reads aligned:", total_reads)
    print("total transcripts FSM:", total_transcript_fsm)
    print("total transcripts NIC:", total_transcript_nic)
    print("total transcripts ISM:", total_transcript_ism)
    print("total transcripts NNC:", total_transcript_nnc)
    print("total transcripts no splice sites:", total_transcript_no_splices)
    print("total splice sites:", all_splice)
    print("GT-AG splice sites:", canonical_splice)

    return read_annotations


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
    torkel_primary_locations = decide_primary_locations(args.torkel_sam, args)
    mm2_primary_locations = decide_primary_locations(args.mm2_sam, args)

    if args.simulated:
        error_rates = get_error_rates(reads)
    else:
        error_rates = {}


    ############# Splice site analysis ########################################
    ###########################################################################

    # print(annotated_ref_isoforms)
    # print(annotated_splice_coordinates)
    print("SHORTEST INTRON:", minimum_annotated_intron)
    minimum_annotated_intron = max(minimum_annotated_intron,  args.min_intron)
    mm2_splice_sites = get_read_candidate_splice_sites(mm2_primary_locations, minimum_annotated_intron, annotated_splice_coordinates_pairs)
    torkel_splice_sites = get_read_candidate_splice_sites(torkel_primary_locations, minimum_annotated_intron, annotated_splice_coordinates_pairs)

    refs = { acc.split()[0] : seq for i, (acc, (seq, _)) in enumerate(readfq(open(args.refs, 'r')))}

    print('MINIMAP2')
    mm2_splice_results = get_splice_classifications(annotated_ref_isoforms, annotated_splice_coordinates, annotated_splice_coordinates_pairs, mm2_splice_sites, refs, mm2_primary_locations)
    print('uLTRA')
    torkel_splice_results = get_splice_classifications(annotated_ref_isoforms, annotated_splice_coordinates, annotated_splice_coordinates_pairs, torkel_splice_sites, refs, torkel_primary_locations)
    # reads_to_cluster_size = get_cluster_sizes(args.cluster_file, reads)

    reads_unaligned_in_torkel = set(reads.keys()) - set(torkel_primary_locations.keys())
    reads_unaligned_in_mm2 = set(reads.keys()) - set(mm2_primary_locations.keys()) 

    detailed_results_outfile = open(os.path.join(args.outfolder, "results_per_read.csv"), "w")

    detailed_results_outfile.write("acc,read_type,error_rate,read_length,is_unaligned_in_other_method,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag\n")
    print_detailed_values_to_file(error_rates, mm2_splice_results, reads, detailed_results_outfile, reads_unaligned_in_torkel, "minimap2", mm2_primary_locations)    
    print_detailed_values_to_file(error_rates, torkel_splice_results, reads, detailed_results_outfile, reads_unaligned_in_mm2, "uLTRA", torkel_primary_locations)

    detailed_results_outfile.close()

    print()
    print("Reads successfully aligned (original/corrected):", len(torkel_primary_locations),len(mm2_primary_locations))
    print("Total reads", len(reads))
    print("READS UNALIGNED (ORIGINAL/CORRECTED):", len(reads_unaligned_in_torkel), len(reads_unaligned_in_mm2) )

    ###########################################################################
    ###########################################################################




    # print(torkel_sorted)
    # print(",".join([s for s in torkel_stats]))
    # print(",".join([s for s in mm2_stats]))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('torkel_sam', type=str, help='Path to the original read file')
    parser.add_argument('mm2_sam', type=str, help='Path to the corrected read file')
    parser.add_argument('reads', type=str, help='Path to the read file')
    parser.add_argument('refs', type=str, help='Path to the refs file')
    parser.add_argument('gff_file', type=str, help='Path to the refs file')
    parser.add_argument('outfolder', type=str, help='Output path of results')
    parser.add_argument('--min_intron', type=int, default=15, help='Threchold for what is counted as varation/intron in alignment as opposed to deletion.')
    parser.add_argument('--infer_genes', action= "store_true", help='Include pairwise alignment of original and corrected read.')
    parser.add_argument('--load_database', action= "store_true", help='Load already computed splice junctions and transcript annotations instead of constructing a new database.')
    parser.add_argument('--simulated', action= "store_true", help='Adds extra analysis that can be done for simulated data since known true locations.')

    # parser.add_argument('--align', action= "store_true", help='Include pairwise alignment of original and corrected read.')

    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)
