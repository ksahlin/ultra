
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
    # genes_to_ref = {} # gene_id : { (exon_start, exon_stop) : set() }
    exons_to_ref = {} # gene_id : { (exon_start, exon_stop) : set() }
    exon_id_to_choordinates = {}
    splices_to_transcripts = defaultdict(dict)
    all_splice_pairs_annotations = defaultdict(lambda: defaultdict(set))
    all_splice_sites_annotations = defaultdict(set)
    # annotated_transcripts = defaultdict(set)

    parts_to_exons = defaultdict(lambda: defaultdict(set))
    for i, exon in enumerate(db.features_of_type('exon', order_by='seqid')):
        if exon.seqid.isdigit() or exon.seqid == 'X' or exon.seqid == 'Y':
            chr_id = 'chr'+ exon.seqid
        elif exon.seqid == 'MT':
            chr_id = 'chrM'
        else:
            chr_id = exon.seqid

        exons_to_ref[exon.id] = chr_id
        exon_id_to_choordinates[exon.id] = (exon.start - 1, exon.stop)
        # creating the augmentation
        if i == 0: # initialization
            prev_seq_id = chr_id
            active_start = exon.start - 1
            active_stop = exon.stop
            active_exons = set() 
            active_exons.add(exon.id)

        if chr_id != prev_seq_id: # switching chromosomes
            parts_to_exons[prev_seq_id][(active_start, active_stop)] = active_exons
            prev_seq_id = chr_id
            active_start = exon.start - 1
            active_stop = exon.stop
            active_exons = set() 
            active_exons.add(exon.id)           

        if exon.start - 1 > active_stop:
            parts_to_exons[chr_id][(active_start, active_stop)] = active_exons
            active_exons = set()
            active_exons.add(exon.id)

            active_start = exon.start - 1
            active_stop = exon.stop

        else:
            active_exons.add(exon.id)    
            active_stop = max(active_stop, exon.stop)

        assert active_start <= exon.start - 1

    parts_to_exons[chr_id][(active_start, active_stop)] = active_exons

    # print(parts_to_exons["SIRV3"])
    # sys.exit()
    print()
    for sid in parts_to_exons:
        print(sid, len(parts_to_exons[sid]))
        print()

    # for transcript in db.features_of_type('transcript'):
    #     genes_to_ref[transcript.id] = str(transcript.seqid)
    #     print("here", transcript.id,  str(transcript.seqid))  
    # for gene in db.features_of_type('gene'):
    #     genes_to_ref[gene.id] = str(gene.seqid)
    #     print("here", gene.id,  str(gene.seqid))  
    # exons_to_parts = reverse_mapping(parts_to_exons)
    # sys.exit()
    # ccc = set()
    # for gene in db.features_of_type('gene'):
    #     genes_to_ref[gene.id] = str(gene.seqid)
    #     print("here", gene.id,  str(gene.seqid))
    for transcript in db.features_of_type('transcript', order_by='seqid'): #db.children(gene, featuretype='transcript', order_by='start'):
        if transcript.seqid.isdigit() or transcript.seqid == 'X' or transcript.seqid == 'Y':
            chr_id = 'chr'+ transcript.seqid
        elif transcript.seqid == 'MT':
            chr_id = 'chrM'
        else:
            chr_id = transcript.seqid

        # print("here", transcript.id,  str(chr_id))  
        transcript_exons = []
        for exon in db.children(transcript, featuretype='exon', order_by='start'):
            transcript_exons.append( (exon.start-1, exon.stop) )
        internal_transcript_splices = [ (e1[1],e2[0]) for e1, e2 in zip(transcript_exons[:-1],transcript_exons[1:])]
        
        # internal transcript splices
        splices_to_transcripts[chr_id][ tuple(internal_transcript_splices)] = transcript.id
        for site1, site2 in internal_transcript_splices:
            all_splice_pairs_annotations[str(chr_id)][(site1, site2)].add( transcript.id )
            # if site2 == 3105:
            #     print('LOOOL',chr_id)
            #     ccc.add(chr_id)

            all_splice_sites_annotations[str(chr_id)].add(site1)
            all_splice_sites_annotations[str(chr_id)].add(site2)
        
        # add start and end splice to all_splice_sites_annotations 
        all_splice_sites_annotations[str(chr_id)].add(transcript_exons[0][0])
        all_splice_sites_annotations[str(chr_id)].add(transcript_exons[-1][-1])

    # if ccc:
    #     print(ccc)
    #     sys.exit()
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
    return  exons_to_ref, parts_to_exons, splices_to_transcripts, transcripts_to_splices, all_splice_pairs_annotations, all_splice_sites_annotations, exon_id_to_choordinates



def get_part_sequences_from_choordinates(parts_to_exons, refs):
    segments = {}
    for chr_id in parts_to_exons:
        parts_instance = parts_to_exons[chr_id]
        # chromosome = genes_to_ref[chr_id]
        segments[chr_id] = {}
        for part in parts_instance:
            start,stop = part[0], part[1]
            seq = refs[chr_id][start : stop] 
            segments[chr_id][part] = seq
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




    