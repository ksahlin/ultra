import pysam


def get_segments(read_aln, ref_aln, predicted_exons):
    segments = []
    ref_seq_break_points = set()
    prev = 0
    print(predicted_exons)
    for p1,p2 in predicted_exons:
        ref_seq_break_points.add( p2 - p1 + prev )
        prev += p2-p1
    
    print('ref_seq_break_points', ref_seq_break_points)

    curr_ref_pos = 0
    ref_aln_break_points = []
    for i, n in enumerate(ref_aln):
        if curr_ref_pos in ref_seq_break_points:
            ref_aln_break_points.append(i)

        if n != '-':
            curr_ref_pos += 1
    if n != '-' and curr_ref_pos in ref_seq_break_points:
        ref_aln_break_points.append(i)

    print('ref_aln_break_points', ref_aln_break_points)

    e_start = 0
    for e_stop in ref_aln_break_points:
        segments.append( (read_aln[ e_start : e_stop ], ref_aln[ e_start : e_stop ]) )
        e_start = e_stop
    return segments

def get_type(n1, n2):
    if n1 == n2:
        return '='
    elif n1 == '-':
        return 'D'
    elif n2 == '-':
        return 'I'
    else:
        return 'X'

def get_cigars(segments):
    segments_cigar = []
    for read, ref in segments:
        c = []
        prev_type = get_type(read[0], ref[0])
        length = 1
        for n1,n2 in zip(read[1:], ref[1:]):
            curr_type = get_type(n1, n2)
            if curr_type == prev_type:
                length += 1
            else:
                c.append(str(length) + prev_type)
                length = 1
                prev_type = curr_type
        

        c.append(str(length) + prev_type)
        segment_cigar = "".join([item for item in c ])
        segments_cigar.append(segment_cigar)

    print('segments_cigar', segments_cigar)
    return segments_cigar #"".join([str(length)+ type_ for length, type_ in c ])

def get_genomic_cigar(read_aln, ref_aln, predicted_exons):

    segments = get_segments(read_aln, ref_aln, predicted_exons)
    cigars = get_cigars(segments)
    genomic_cigar = []
    intron_lengths = [e2[0] - e1[1] for e1, e2 in zip(predicted_exons[:-1], predicted_exons[1:])]
    for i in range(len(cigars)):
        if i <= len(intron_lengths) -1:
            genomic_cigar.append( cigars[i] + '{0}N'.format( intron_lengths[i] ) )
        else:
            genomic_cigar.append( cigars[i]  )

    genomic_cigar = "".join(s for s in genomic_cigar)
    return genomic_cigar



def main(read_id, ref_id, classification, predicted_exons, read_aln, ref_aln, annotated_to_transcript_id, alignment_outfile):
    print(ref_id, classification, predicted_exons, read_aln, ref_aln, alignment_outfile)
    read_sam_entry = pysam.AlignedSegment(alignment_outfile.header)
    if classification != 'unaligned':
        genomic_cigar = get_genomic_cigar(read_aln, ref_aln, predicted_exons)
        print('genomic cigar:', genomic_cigar)
        read_sam_entry.cigarstring = genomic_cigar
        read_sam_entry.reference_start = predicted_exons[0][0]

    else:
        read_sam_entry.cigarstring = '*'
        read_sam_entry.reference_start = -1

    read_sam_entry.query_name = read_id
    read_sam_entry.flag = 0 # TODO: add reverse complements
    read_sam_entry.reference_name = ref_id
    read_sam_entry.mapping_quality = 60 # TODO: calculate mapping quality 
    print(annotated_to_transcript_id)
    read_sam_entry.set_tag('AN', annotated_to_transcript_id)
    read_sam_entry.set_tag('CN', classification)

    # read_sam_entry.reference_star = 





    alignment_outfile.write(read_sam_entry)
    # sys.exit()