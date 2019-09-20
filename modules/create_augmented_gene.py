
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
    exons_to_ref = {} # gene_id : { (exon_start, exon_stop) : set() }
    exon_id_to_choordinates = {}
    splices_to_transcripts = defaultdict(dict)
    all_splice_pairs_annotations = defaultdict(lambda: defaultdict(set))
    all_splice_sites_annotations = defaultdict(set)
    # annotated_transcripts = defaultdict(set)

    parts_to_exons = defaultdict(lambda: defaultdict(set))
    for i, exon in enumerate(db.features_of_type('exon', order_by='seqid')):
        exons_to_ref[exon.id] = exon.seqid
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

    # print(parts_to_exons["SIRV3"])
    # sys.exit()

    for sid in parts_to_exons:
        print(sid, parts_to_exons[sid])
        print()
    
    # exons_to_parts = reverse_mapping(parts_to_exons)
    # sys.exit()

    for gene in db.features_of_type('gene'):
        genes_to_ref[gene.id] = str(gene.seqid)

        for transcript in db.children(gene, featuretype='transcript', order_by='start'):
            transcript_exons = []
            for exon in db.children(transcript, featuretype='exon', order_by='start'):
                transcript_exons.append( (exon.start-1, exon.stop) )
            transcript_splices = [ (e1[1],e2[0]) for e1, e2 in zip(transcript_exons[:-1],transcript_exons[1:])]

            splices_to_transcripts[gene.seqid][ tuple(transcript_splices)] = transcript.id
            for site1, site2 in transcript_splices:
                all_splice_pairs_annotations[str(gene.seqid)][(site1, site2)].add( transcript.id )
                all_splice_sites_annotations[str(gene.seqid)].add(site1)
                all_splice_sites_annotations[str(gene.seqid)].add(site2)


    transcripts_to_splices = reverse_mapping(splices_to_transcripts)

            # print("transcript_parts", tuple(transcript_parts))
    # print(all_splice_pairs_annotations)
    # print(all_splice_sites_annotations)
    # print(splices_to_transcripts)
    # print(parts_to_exons)
    # sys.exit()
    # for sid in parts_to_exons:
    #     print(sid, parts_to_exons[sid])
    #     print()
    # sys.exit()
    # print(splices_to_transcripts)
    return  genes_to_ref, exons_to_ref, parts_to_exons, splices_to_transcripts, transcripts_to_splices, all_splice_pairs_annotations, all_splice_sites_annotations, exon_id_to_choordinates



def get_part_sequences_from_choordinates(parts_to_exons, genes_to_ref, refs):
    segments = {}
    for gene_id in parts_to_exons:
        parts_instance = parts_to_exons[gene_id]
        chromosome = genes_to_ref[gene_id]
        segments[chromosome] = {}
        for part in parts_instance:
            start,stop = part[0], part[1]
            seq = refs[chromosome][start : stop] 

            segments[chromosome][part] = seq
    # print(segments)
    return segments

def get_exon_sequences_from_choordinates(exon_id_to_choordinates, exons_to_ref, refs):
    exon_sequences = defaultdict(dict)
    for exon_id in exon_id_to_choordinates:
        start,stop = exon_id_to_choordinates[exon_id]
        chromosome = exons_to_ref[exon_id]
        seq = refs[chromosome][start : stop] 
        exon_sequences[chromosome][(start,stop)] = seq
    # print(segments)
    return exon_sequences




    