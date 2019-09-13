
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
    parts_to_exons = {}
    exon_id_to_choordinates = {}
    splices_to_transcripts = defaultdict(dict)
    parts_to_transcript_annotations = defaultdict(lambda: defaultdict(set))
    transcripts_to_parts_annotations = defaultdict(lambda: defaultdict(set))
    all_parts_pairs_annotations = defaultdict(lambda: defaultdict(set))
    all_part_sites_annotations = defaultdict(set)
    # annotated_transcripts = defaultdict(set)
    for gene in db.features_of_type('gene'):
        # print(dir(gene))
        # print(gene.id, gene.seqid, gene.start, gene.stop, gene.attributes)
        genes_to_ref[gene.id] = str(gene.seqid)
        # exon_id_to_choordinates[str(gene.seqid)] = {}
        parts_to_exons[str(gene.seqid)] = defaultdict(set)
        # splices_to_transcripts[gene.seqid] = {}
        #add nodes
        exons = [exon for exon in db.children(gene, featuretype='exon', order_by='start') ]
        chord_to_exon = defaultdict(list)

        parts_to_exons_for_gene = {}        
        active_exons = set() #set(chord_to_exon[all_starts_and_stops[0][0]])

        for i, e in enumerate(exons):
            chord_to_exon[e.start - 1].append(e.id)
            chord_to_exon[e.stop].append(e.id)
            exon_id_to_choordinates[e.id] = (e.start - 1, e.stop)

            # creating the augmentation
            if i == 0:
                active_start = e.start - 1
                active_stop = e.stop
                active_exons.add(e.id)

            if e.start - 1 > active_stop:
                parts_to_exons_for_gene[(active_start, active_stop)] = active_exons
                active_exons = set()
                active_exons.add(e.id)

                active_start = e.start - 1
                active_stop = e.stop

            else:
                active_exons.add(e.id)    
                active_stop = e.stop

        parts_to_exons_for_gene[(active_start, active_stop)] = active_exons

        # print(parts_to_exons_for_gene)

        exon_to_chord = {e.id : (e.start-1, e.stop) for e in exons}
        # print([(e.start-1, e.stop) for e in exons])



        # print("PARTS to exons",  parts_to_exons_for_gene)
        exons_to_parts = reverse_mapping(parts_to_exons_for_gene)
        # print("EXONS to PARTS",  exons_to_parts)

        # extend the parts that are smaller than min_mem to length min_mem + 1
        # parts_to_exons[gene.id] = parts_to_exons_for_gene
        for (start,stop), active_exons in parts_to_exons_for_gene.items():
            if (start,stop) in parts_to_exons[str(gene.seqid)]:
                parts_to_exons[str(gene.seqid)][(start,stop)].update(active_exons)
            else:
                parts_to_exons[str(gene.seqid)][(start,stop)] = active_exons

        # for exon in db.children(gene, featuretype='exon', order_by='start'):
        #     splices_to_transcripts[gene.id][ (exon.start, exon.stop) ].update([ transcript_tmp for transcript_tmp in  exon.attributes['transcript_id']])

        # for transcript in db.children(gene, featuretype='transcript', order_by='start'):
        #     annotated_transcripts[gene.seqid].add( tuple( '_'.join([str(item) for item in (gene.seqid, exon.start, exon.stop)]) for exon in db.children(transcript, featuretype='exon', order_by='start') ) )

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

    # print(splices_to_transcripts)
    return  genes_to_ref, parts_to_exons, splices_to_transcripts, parts_to_transcript_annotations, transcripts_to_parts_annotations,  all_parts_pairs_annotations, all_part_sites_annotations, exon_id_to_choordinates #annotated_transcripts



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




    