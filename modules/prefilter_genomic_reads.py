
import os 
import subprocess
import sys

from struct import *
import intervaltree
import pysam

from collections import defaultdict

from modules import help_functions

def get_ultra_indexed_choordinates(ref_part_sequences):
    indexed_regions = defaultdict(intervaltree.IntervalTree)

    for sequence_id, seq  in ref_part_sequences.items():
        chr_id, start, stop = unpack('LLL',sequence_id)
        indexed_regions[chr_id].addi(start, stop, None)

    return indexed_regions

def align_with_minimap2(refs_path, read_path, outfolder, nr_cores):
    minimap2_samfile_path = os.path.join(outfolder, "minimap2.sam")
    with open(minimap2_samfile_path, "w") as output_file:
        # print('Running spoa...', end=' ')
        sys.stdout.flush()
        stderr_file = open(os.path.join(outfolder, "minimap2_errors.1") , 'w')
        # minimap2 --eqx -t 62 -ax splice -k13 -w 5 -G 500k {input.index} {input.fastq}
        subprocess.check_call([ 'minimap2',  '-ax', 'splice', '--eqx' , '-k13', '-w5', '-G', '500k', '-t', str(nr_cores), refs_path, read_path], stdout=output_file, stderr=stderr_file)
        # print('Done.')
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

        else: # reference skip or soft/hardclip "~", or match =
            print("UNEXPECTED!", t)
            sys.exit()
    exon_sites.append( (exon_start, ref_pos) )
    return exon_sites




def parse_alignments_and_mask(minimap2_samfile_path, indexed_regions):
    SAM_file = pysam.AlignmentFile(minimap2_samfile_path, "r", check_sq=False)
    reads_to_ignore = {}

    for read in SAM_file.fetch(until_eof=True):
        if read.flag == 0 or read.flag == 16:

            aligned_choordinates = get_exons_from_cigar(read)
            total_overlap = 0
            total_aligned_length = 0
            for (a_start, a_stop) in aligned_choordinates:
                overlaps =  indexed_regions[read.reference_name].overlap(a_start, a_stop)
                for ovl in overlaps: continue here to get sizxe of overlap
                    total_overlap += overlaps
                print("overlaps", overlaps)
                total_overlap += overlaps
                total_aligned_length += a_stop - a_start

            # is_exonic = 1 if indexed_regions[reference_name].overlaps(reference_start, reference_end) else 0

            if total_overlap/total_aligned_length < 0.7:
                reads_to_ignore[read.query_name] = read
    return reads_to_ignore


def print_read_categories(reads_to_ignore, reads, outfolder):
    reads_to_align = open(os.path.join(outfolder, "reads_after_genomic_filtering.fasta"), "w")
    genomic_aligned = pysam.AlignmentFile(open(os.path.join(outfolder, "genomic.sam"), "w", template=samfile))

    for acc, (seq, _) in help_functions.readfq(open(reads,"r")):
        if acc in reads_to_ignore:
            genomic_aligned.write(read)
        else:
            reads_to_align.write(">{0}\n{1}\n".format(acc, seq))

    genomic_aligned.close()



def main(ref_part_sequences, ref, reads, outfolder, nr_cores):
    indexed_regions = get_ultra_indexed_choordinates(ref_part_sequences)
    minimap2_samfile_path = align_with_minimap2(ref, reads, outfolder, nr_cores)
    reads_to_ignore = parse_alignments_and_mask(minimap2_samfile_path, indexed_regions)
    path_reads_to_align = print_read_categories(reads_to_ignore, reads, outfolder)
    return set(reads_to_ignore.keys()), path_reads_to_align






