
from collections import defaultdict

from modules.help_functions import readfq

def reverse_mapping(d):
    d_dict = defaultdict(list)
    for k,v in sorted(d.items(), key = lambda x: x[0]):
        for i in v:
            d_dict[i].append(k)
    return d_dict


def create_graph_from_exon_parts(db, min_mem): 
    """
        We need to link parts --> exons and exons --> transcripts
    """
    # print(dir(db))
    genes_to_ref = {} # gene_id : { (exon_start, exon_stop) : set() }
    exon_id_to_choordinates = {}
    splices_to_transcripts = defaultdict(dict)
    parts_to_transcript_annotations = defaultdict(lambda: defaultdict(set))
    transcripts_to_parts_annotations = defaultdict(lambda: defaultdict(set))
    all_parts_pairs_annotations = defaultdict(lambda: defaultdict(set))
    all_part_sites_annotations = defaultdict(set)
    # annotated_transcripts = defaultdict(set)

    parts_to_exons = defaultdict(lambda: defaultdict(set))
    for i, exon in enumerate(db.features_of_type('exon', order_by='seqid')):
        # print(exon.id, exon.start, exon.stop, exon.seqid)
        exon_id_to_choordinates[exon.id] = (exon.start - 1, exon.stop)
        # creating the augmentation
        if i == 0: # initialization
            prev_seq_id = exon.seqid
            active_start = exon.start - 1
            active_stop = exon.stop
            active_exons = set() 
            active_exons.add(exon.id)

        if exon.seqid != prev_seq_id: # switching chromosomes
            parts_to_exons[prev_seq_id][(active_start, active_stop)] = active_exons
            prev_seq_id = exon.seqid
            active_start = exon.start - 1
            active_stop = exon.stop
            active_exons = set() 
            active_exons.add(exon.id)           

        if exon.start - 1 > active_stop:
            parts_to_exons[exon.seqid][(active_start, active_stop)] = active_exons
            active_exons = set()
            active_exons.add(exon.id)

            active_start = exon.start - 1
            active_stop = exon.stop

        else:
            active_exons.add(exon.id)    
            active_stop = max(active_stop, exon.stop)

        assert active_start <= exon.start - 1

    parts_to_exons[exon.seqid][(active_start, active_stop)] = active_exons

    for sid in parts_to_exons:
        print(sid, parts_to_exons[sid])
        print()
    
    exons_to_parts = reverse_mapping(parts_to_exons)

    # sys.exit()

    for gene in db.features_of_type('gene'):
        genes_to_ref[gene.id] = str(gene.seqid)

        for transcript in db.children(gene, featuretype='transcript', order_by='start'):
            transcript_parts = []
            transcript_exons = []
            for exon in db.children(transcript, featuretype='exon', order_by='start'):
                transcript_parts +=  exons_to_parts[exon.id]
                transcript_exons.append( (exon.start-1, exon.stop) )

            transcript_splices = [ (e1[1],e2[0]) for e1, e2 in zip(transcript_exons[:-1],transcript_exons[1:])]

            splices_to_transcripts[gene.seqid][ tuple(transcript_splices)] = transcript.id
            parts_to_transcript_annotations[gene.seqid][ tuple(transcript_parts) ].add(  transcript.id )
            transcripts_to_parts_annotations[gene.seqid][ transcript.id ].add( tuple(transcript_parts)  )
            for part_start, part_stop in transcript_parts:
                all_parts_pairs_annotations[str(gene.seqid)][( part_start, part_stop )].add( transcript.id )
                all_part_sites_annotations[str(gene.seqid)].add(part_start)
                all_part_sites_annotations[str(gene.seqid)].add(part_stop)

            # print("transcript_parts", tuple(transcript_parts))

    # print(parts_to_exons)
    # sys.exit()
    # for sid in parts_to_exons:
    #     print(sid, parts_to_exons[sid])
    #     print()
    # sys.exit()
    # print(splices_to_transcripts)
    return  genes_to_ref, parts_to_exons, splices_to_transcripts, parts_to_transcript_annotations, transcripts_to_parts_annotations,  all_parts_pairs_annotations, all_part_sites_annotations, exon_id_to_choordinates 



def get_sequences_from_choordinates(parts_to_exons, genes_to_ref, ref):
    refs = {acc : seq for acc, (seq, _) in readfq(open(ref,"r"))}
    segments = {}
    for gene_id in parts_to_exons:
        parts_instance = parts_to_exons[gene_id]
        chromosome = genes_to_ref[gene_id]
        segments[chromosome] = {}
        for part in parts_instance:
            start,stop = part[0], part[1]
            seq = refs[chromosome][start : stop] # gtf 1 indexed and last coordinate is inclusive

            segments[chromosome][part] = seq
    # print(segments)
    return segments




    