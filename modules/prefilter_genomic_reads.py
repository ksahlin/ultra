
import os 
import subprocess
import sys

from struct import *
import intervaltree
import pysam

from collections import defaultdict

from modules import help_functions

def get_ultra_indexed_choordinates(ref_part_sequences, indexfolder, outfolder):
    id_to_chr = help_functions.pickle_load( os.path.join(indexfolder, 'id_to_chr.pickle') )
    indexed_regions = defaultdict(intervaltree.IntervalTree)

    for sequence_id, seq  in ref_part_sequences.items():
        chr_id, start, stop = unpack('LLL',sequence_id)
        ref_seq_name = id_to_chr[chr_id]
        indexed_regions[ref_seq_name].addi(start, stop, None)

    return indexed_regions

def align_with_minimap2(refs_path, read_path, outfolder, nr_cores, k_size):
    minimap2_samfile_path = os.path.join(outfolder, "minimap2.sam")
    with open(minimap2_samfile_path, "w") as output_file:
        # print('Running spoa...', end=' ')
        sys.stdout.flush()
        stderr_file = open(os.path.join(outfolder, "minimap2_errors.1") , 'w')
        # minimap2 --eqx -t 62 -ax splice -k13 -w 5 -G 500k {input.index} {input.fastq}
        subprocess.check_call([ 'minimap2',  '-ax', 'splice', '--eqx' , '-k', str(k_size), '-t', str(nr_cores), refs_path, read_path], stdout=output_file, stderr=stderr_file)
        sys.stdout.flush()
    output_file.close()
    return minimap2_samfile_path

def get_exons_from_cigar(read):
    # https://pysam.readthedocs.io/en/latest/api.html
    # M   BAM_CMATCH  0
    # I   BAM_CINS    1
    # D   BAM_CDEL    2
    # N   BAM_CREF_SKIP   3
    # S   BAM_CSOFT_CLIP  4
    # H   BAM_CHARD_CLIP  5
    # P   BAM_CPAD    6
    # =   BAM_CEQUAL  7
    # X   BAM_CDIFF   8
    # B   BAM_CBACK   9

    aligned_choordinates = set()
    q_start = read.reference_start
    q_end = read.reference_end
    exon_sites = []
    ref_pos = q_start
    exon_start = q_start
    aln_cig_types = {2, 7, 0, 8} # 

    for i, (t, l) in enumerate(read.cigartuples):
        if t in aln_cig_types:
            ref_pos += l
        elif t == 3:
            exon_sites.append( (exon_start, ref_pos) )
            ref_pos += l
            exon_start = ref_pos

        # elif t == "I" or t == "S" or t == "H": # insertion or softclip
        #     ref_pos += 0
        # else: # reference skip or soft/hardclip "~", or match =
        #     print("UNEXPECTED!", t)
        #     sys.exit()

    exon_sites.append( (exon_start, ref_pos) )
    return exon_sites


def overlap_size(a, b, c, d):
    max_start = max(a, c)
    min_stop = min(b, d)
    return min_stop - max_start

def is_overlapping(a_start,a_stop, b_start,b_stop):
    return (int(a_start) <= int(b_start) <= int(a_stop) )  or (int(a_start) <= int(b_stop) <= int(a_stop)) or (int(b_start) <= int(a_start) <= int(a_stop) <= int(b_stop) )


def filter_reads_to_align(minimap2_samfile_path, indexed_regions, outfolder, genomic_frac_cutoff):
    SAM_file = pysam.AlignmentFile(minimap2_samfile_path, "r", check_sq=False)

    reads_to_align = open(os.path.join(outfolder, "reads_after_genomic_filtering.fastq"), "w")
    unindexed_aligned = pysam.AlignmentFile(os.path.join(outfolder, "unindexed.sam"), "w", template=SAM_file)
    indexed_aligned = pysam.AlignmentFile(os.path.join(outfolder, "indexed.sam"), "w", template=SAM_file)

    # reads_unindexed = {}
    # reads_indexed = {}
    nr_reads_unindexed = 0 

    for read in SAM_file.fetch(until_eof=True):
        qual_chars = pysam.qualities_to_qualitystring( read.query_qualities )
        if read.flag == 0 or read.flag == 16:

            aligned_choordinates = get_exons_from_cigar(read)
            total_overlap = 0
            total_aligned_length = 0
            for (a_start, a_stop) in aligned_choordinates:
                overlaps = indexed_regions[read.reference_name].overlap(a_start, a_stop)
                for ovl in overlaps: #continue here to get sizxe of overlap
                    total_overlap += overlap_size(a_start, a_stop, ovl.begin, ovl.end)
                total_aligned_length += a_stop - a_start

            if 1 - (total_overlap/total_aligned_length) > genomic_frac_cutoff:
                read.set_tag('XA', "")
                read.set_tag('XC', "uLTRA_unindexed")
                unindexed_aligned.write(read)
                # reads_unindexed[read.query_name] = read
                nr_reads_unindexed += 1
            else:
                seq = help_functions.reverse_complement(read.query_sequence) if read.is_reverse else read.query_sequence
                reads_to_align.write(">{0}\n{1}\n+\n{2}\n".format(read.query_name, seq, qual_chars))
                indexed_aligned.write(read)
                # reads_indexed[read.query_name] = read
        
        elif read.flag == 4: # unmapped
            reads_to_align.write(">{0}\n{1}\n+\n{2}\n".format(read.query_name, read.query_sequence, qual_chars))


    unindexed_aligned.close()
    indexed_aligned.close()
    reads_to_align.close()

    return nr_reads_unindexed, reads_to_align.name


def main(ref_part_sequences, ref, reads, outfolder, indexfolder, nr_cores, genomic_frac_cutoff, k_size):
    indexed_regions = get_ultra_indexed_choordinates(ref_part_sequences, indexfolder, outfolder)
    minimap2_samfile_path = align_with_minimap2(ref, reads, outfolder, nr_cores, k_size)
    print('minimap2 done.')
    nr_reads_unindexed, path_reads_to_align = filter_reads_to_align(minimap2_samfile_path, indexed_regions, outfolder, genomic_frac_cutoff)
    return nr_reads_unindexed, path_reads_to_align





