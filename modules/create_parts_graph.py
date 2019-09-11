
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
        We need to link parts --> exons --> transcripts
    """
    # print(dir(db))
    genes_to_ref = {} # gene_id : { (exon_start, exon_stop) : set() }
    parts_to_exons = {}
    exon_choordinates = {}
    exons_to_transcripts = defaultdict(dict)
    parts_to_transcript_annotations = defaultdict(lambda: defaultdict(set))
    transcripts_to_parts_annotations = defaultdict(lambda: defaultdict(set))
    all_parts_pairs_annotations = defaultdict(lambda: defaultdict(set))
    all_part_sites_annotations = defaultdict(set)
    # annotated_transcripts = defaultdict(set)
    for gene in db.features_of_type('gene'):
        # print(dir(gene))
        # print(gene.id, gene.seqid, gene.start, gene.stop, gene.attributes)
        genes_to_ref[gene.id] = str(gene.seqid)
        exon_choordinates[str(gene.seqid)] = {}
        # parts_to_exons[gene.id] = defaultdict(set)
        # exons_to_transcripts[gene.seqid] = {}
        parts_to_exons_for_gene = {}        
        #add nodes
        exons = [exon for exon in db.children(gene, featuretype='exon', order_by='start') ]
        chord_to_exon = defaultdict(list)
        for e in exons:
            chord_to_exon[e.start - 1].append(e.id)
            chord_to_exon[e.stop].append(e.id)
            exon_choordinates[str(gene.seqid)][ (e.start - 1, e.stop)] = e.id


        exon_to_chord = {e.id : (e.start-1, e.stop) for e in exons}
        print([(e.start-1, e.stop) for e in exons])


        all_starts = [(e.start -1, 'start') for e in exons]
        all_stops = [(e.stop, 'stop') for e in exons]
        # all_starts_and_stops = sorted( set(all_starts + all_stops))
        all_starts_and_stops = sorted( set(all_starts + all_stops), key = lambda x: (x[0],len(x[1])) ) # all stops before starts if the same coordinate

        print()
        print(str(gene.seqid))
        print()
        print(all_starts_and_stops)
        active_exons = set() #set(chord_to_exon[all_starts_and_stops[0][0]])
        for p1, p2 in zip(all_starts_and_stops[:-1], all_starts_and_stops[1:]):
            if p1[1] == 'stop':
                active_exons =  active_exons - set(chord_to_exon[p1[0]])
            elif p1[1] == 'start':
                active_exons =  active_exons | set(chord_to_exon[p1[0]])
            
            # if p1[0] == 11122 or p2[0] == 11153:
            #     print(p1, p2)
            #     print(all_starts_and_stops2)
            #     sys.exit()


            if active_exons:  
            # if (p1[1], p2[1]) != ('stop', 'start'):
            #     if p1[1] == 'start' and p2[1] == 'stop':
            #         active_exons = set(chord_to_exon[p2[0]]) | set(chord_to_exon[p1[0]])
            #     elif p1[1] == 'stop' and p2[1] == 'stop':
            #         active_exons = set(chord_to_exon[p2[0]]) - set(chord_to_exon[p1[0]])
            #     elif p1[1] == 'start' and p2[1] == 'start':
            #         active_exons = set(chord_to_exon[p1[0]]) - set(chord_to_exon[p2[0]])

                part_length = int(p2[0]) - int(p1[0]) 
                if part_length < min_mem:
                    if p1[1] == 'start' and p2[1] == 'stop':
                        print('Need to extend over junction because exon smaller tham min mem. Treat this case!')

                        active_start_and_stop =  set(chord_to_exon[p2[0]]) & set(chord_to_exon[p1[0]])
                        active_starts_only =  set(chord_to_exon[p1[0]]) - active_start_and_stop
                        active_stops_only = set(chord_to_exon[p2[0]]) -  active_start_and_stop 
                        active_both_directions = active_exons - (active_start_and_stop - active_starts_only - active_stops_only)
                        print(active_start_and_stop)
                        print(active_starts_only)
                        print(active_stops_only)
                        print(active_both_directions)

                        # sys.exit()
                        if active_starts_only:
                            e_ids = [e_id for e_id in active_starts_only ]
                            tightest_bound_downstream = min([ exon_to_chord[e_id][1] for e_id in e_ids])
                            print(p1,p2, tightest_bound_downstream)                 
                            print("Tackle this!! HERE")

                            if tightest_bound_downstream - p2[0] > min_mem - part_length: # enough room to extend
                                parts_to_exons_for_gene[(p1[0], p2[0] + (min_mem - part_length))] = active_starts_only  | active_both_directions
                                print("MANAGED downstream!")
                                # sys.exit()

                        if active_stops_only:

                            e_ids = [e_id for e_id in active_stops_only ]
                            tightest_bound_upstream = max([ exon_to_chord[e_id][0] for e_id in e_ids])
                            print(p1,p2, tightest_bound_upstream)

                            if p1[0] - tightest_bound_upstream > min_mem - part_length: # enough room to extend
                                parts_to_exons_for_gene[(p1[0] - (min_mem - part_length), p2[0])] = active_stops_only  | active_both_directions
                                print("MANAGED upstream!", (p1[0] - (min_mem - part_length), p2[0]))
                                # sys.exit()
                        
                        # if active_both_directions:
                        #     e_ids = [e_id for e_id in active_both_directions ]
                        #     tightest_bound_downstream = min([ exon_to_chord[e_id][1] for e_id in e_ids])
                        #     tightest_bound_upstream = max([ exon_to_chord[e_id][0] for e_id in e_ids])

                        #     if (tightest_bound_downstream - p2[0]) + (p1[0] - tightest_bound_upstream) > min_mem - part_length: # enough room to extend
                        #         parts_to_exons_for_gene[(p1[0] - (min_mem - part_length), p2[0] + (min_mem - part_length))] =  active_both_directions
                
                        # if (p1[0] - (min_mem - part_length), p2[0] + (min_mem - part_length)) == (11122, 11153):
                        #     print('looool',p1, p2, active_exons)
                        #     print(exon_to_chord)
                        #     sys.exit()

                        if active_start_and_stop:
                            print("DIDN'T MANAGE!", p1[0], p2[0], )
                            print("failed for short exon")
                            # sys.exit()

                    elif p1[1] == 'stop' and p2[1] == 'stop':
                        print('extending smaller part upstream',p1, p2)
                        # exon_id =  chord_to_exon[p1[0]][0]
                        # furthest_upstream = exon_to_chord[exon_id][0] 
                        e_ids = [e_id for e_id in chord_to_exon[p1[0]] ]
                        tightest_upstream = max([ exon_to_chord[e_id][0] for e_id in e_ids])

                        if p1[0] - tightest_upstream > min_mem - part_length: # enough room to extend
                            parts_to_exons_for_gene[(p1[0] - (min_mem - part_length), p2[0])] = active_exons
                            print("extended:", (p1[0] - (min_mem - part_length), p2[0]) , "active_exons:", active_exons)
                        else:
                            print("Not enough room to extend!!")
                            # sys.exit()

                    elif p1[1] == 'start' and p2[1] == 'start':
                        print('extending smaller part to downstream',p1, p2)
                        # exon_id =  chord_to_exon[p2[0]][0]
                        # furthest_downstream = exon_to_chord[exon_id][1]
                        e_ids = [e_id for e_id in chord_to_exon[p2[0]]]
                        tightest_bound_downstream = min([ exon_to_chord[e_id][1] for e_id in e_ids])

                        if tightest_bound_downstream - p2[0] > min_mem - part_length: # enough room to extend
                            parts_to_exons_for_gene[(p1[0], p2[0] + (min_mem - part_length))] = active_exons
                            print("extended:", (p1[0], p2[0] + (min_mem - part_length)), "active_exons:", active_exons )
                        else:
                            print("Not enough room to extend!!")
                            # sys.exit()

                    elif p1[1] == 'stop' and p2[1] == 'start':
                        e_ids = [e_id for e_id in active_exons ]
                        tightest_bound_upstream = max([ exon_to_chord[e_id][0] for e_id in e_ids])
                        tightest_bound_downstream = min([ exon_to_chord[e_id][1] for e_id in e_ids])
                        print(p1,p2,tightest_bound_upstream, tightest_bound_downstream)                 
                        print("Tackle this!! HERE")

                        if p1[0] - tightest_bound_upstream > min_mem - part_length: # enough room to extend
                            parts_to_exons_for_gene[(p1[0] - (min_mem - part_length), p2[0])] = active_exons
                            print("MANAGED upstream!", (p1[0] - (min_mem - part_length), p2[0]))
                            # sys.exit()

                        elif tightest_bound_downstream - p2[0] > min_mem - part_length: # enough room to extend
                            parts_to_exons_for_gene[(p1[0], p2[0] + (min_mem - part_length))] = active_exons
                            print("MANAGED downstream!")
                            # sys.exit()

                        else:
                            print("DIDN'T MANAGE!", p1[0], p2[0], )

                        # sys.exit()
                    else:
                        print("IM HERE BUG")
                        sys.exit()

                else:
                    parts_to_exons_for_gene[(p1[0], p2[0])] = active_exons


            else:
                print("HEERE:", p1, p2 )

            # Update active set
            if p2[1] == 'stop':
                active_exons = active_exons - set(chord_to_exon[p2[0]])



        print("PARTS to exons",  parts_to_exons_for_gene)
        exons_to_parts = reverse_mapping(parts_to_exons_for_gene)
        print("EXONS to PARTS",  exons_to_parts)

        # extend the parts that are smaller than min_mem to length min_mem + 1
        parts_to_exons[gene.id] = parts_to_exons_for_gene

        # for exon in db.children(gene, featuretype='exon', order_by='start'):
        #     exons_to_transcripts[gene.id][ (exon.start, exon.stop) ].update([ transcript_tmp for transcript_tmp in  exon.attributes['transcript_id']])

        # for transcript in db.children(gene, featuretype='transcript', order_by='start'):
        #     annotated_transcripts[gene.seqid].add( tuple( '_'.join([str(item) for item in (gene.seqid, exon.start, exon.stop)]) for exon in db.children(transcript, featuretype='exon', order_by='start') ) )

        for transcript in db.children(gene, featuretype='transcript', order_by='start'):
            transcript_parts = []
            transcript_exons = []
            for exon in db.children(transcript, featuretype='exon', order_by='start'):
                transcript_parts +=  exons_to_parts[exon.id]
                transcript_exons.append( (exon.start-1, exon.stop) )

            exons_to_transcripts[gene.seqid][ tuple(transcript_exons)] = transcript.id
            parts_to_transcript_annotations[gene.seqid][ tuple(transcript_parts) ].add(  transcript.id )
            transcripts_to_parts_annotations[gene.seqid][ transcript.id ].add( tuple(transcript_parts)  )
            for part_start, part_stop in transcript_parts:
                all_parts_pairs_annotations[str(gene.seqid)][( part_start, part_stop )].add( transcript.id )
                all_part_sites_annotations[str(gene.seqid)].add(part_start)
                all_part_sites_annotations[str(gene.seqid)].add(part_stop)

            # print("transcript_parts", tuple(transcript_parts))

    # sys.exit()

    # print(exons_to_transcripts)
    return  genes_to_ref, parts_to_exons, exons_to_transcripts, parts_to_transcript_annotations, transcripts_to_parts_annotations,  all_parts_pairs_annotations, all_part_sites_annotations,exon_choordinates #annotated_transcripts



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




    