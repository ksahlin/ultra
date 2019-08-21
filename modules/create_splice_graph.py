
from collections import defaultdict
import networkx as nx

def parse_tsv(): 
    pass


def create_graph(db): 

    print(dir(db))
    gene_database = {} # gene_id : { (exon_start, exon_stop) : set() }
    topological_sorts = {} # gene_id : { (exon_start, exon_stop) : set() }
    collapsed_exon_to_transcript = {}
    for gene in db.features_of_type('gene'):
        # print(dir(gene))
        print(gene.id, gene.seqid, gene.start, gene.stop, gene.attributes)
        gene_graph = nx.DiGraph()
        collapsed_exon_to_transcript[gene.id] = defaultdict(set)
        already_parsed_exons = set()
        
        #add nodes
        for exon in db.children(gene, featuretype='exon', order_by='start'):
            collapsed_exon_to_transcript[gene.id][ (exon.start, exon.stop) ].update([ transcript_tmp for transcript_tmp in  exon.attributes['transcript_id']])
            if (exon.start, exon.stop) in already_parsed_exons:
                
                gene_graph.add_node( (exon.start, exon.stop) )

        #add edges
        for transcript in db.children(gene, featuretype='transcript', order_by='start'):
            # print(dir(transcript))
            # print('transcript', transcript.id, transcript.start, transcript.stop)
            consecutive_exons = [exon for exon in db.children(transcript, featuretype='exon', order_by='start')]
            for e1,e2 in zip(consecutive_exons[:-1], consecutive_exons[1:]):
                # print('exon', exon.id, exon.start, exon.stop)
                gene_graph.add_edge( (e1.start, e1.stop),  (e2.start, e2.stop)  )
                

        print(gene_graph.edges())
        gene_database[gene.id] = gene_graph
        
        top_sort = list(nx.topological_sort(gene_graph))
        topological_sorts[gene.id] = top_sort

    return gene_database, topological_sorts

def extract_spanning_paths()

def collapse_identical_for_mumer():
    '''
        Creates a mapping from one segment on the genome to several exons containig the segment, this to remove redundancy in database when calling mummer
    '''
    pass





# def parse_gff_for_transcripts(gff_input): 
#     fn = gffutils.example_filename(gff_input)
#     db = gffutils.create_db(fn, dbfn='test.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
#     db = gffutils.FeatureDB('test.db', keep_order=True)
#     # gene = db["PB.1016"]
#     # print(transcript)
#     ref_isoforms = {}
#     for tr in db.children(gene, featuretype='transcript', order_by='start'):
#         # print(tr.id, dir(tr)) 
#         splice_sites = []
#         for j in db.children(tr, featuretype='exon', order_by='start'):
#             # print(j, j.start, j.end)
#             splice_sites.append(j.start -1)
#             splice_sites.append(j.end)
#         splice_sites_tmp = splice_sites[1:-1]
#         splice_sites = []
#         for i in range(0, len(splice_sites_tmp),2):
#             splice_sites.append( (splice_sites_tmp[i], splice_sites_tmp[i+1]) )
#         # splice_sites = [item for item in zip(splice_sites[:-1], splice_sites[1:])]
#         ref_isoforms[tr.id] = splice_sites
#     print(ref_isoforms)
#     return ref_isoforms





            # if (exon.start, exon.stop) in already_parsed_exons:
            #     print(exon.attributes['transcript_id'])
            #     collapsed_exon_to_transcript[gene.id][ (exon.start, exon.stop) ].update([ transcript_tmp for transcript_tmp in  exon.attributes['transcript_id']])
            #     print(collapsed_exon_to_transcript[gene.id])
            # else:
            #     already_parsed_exons.add( (exon.start, exon.stop) )
            #     collapsed_exon_to_transcript[gene.id][ (exon.start, exon.stop) ].update([ transcript_tmp for transcript_tmp in  exon.attributes['transcript_id']])

            # print('exon', exon.id, exon.start, exon.stop)
            # transcript_id
            # print(dir(exon), exon.attributes)
            # exon_to_transcript[exon.id].add(transcript.id)       


        # for transcript in db.children(gene, featuretype='transcript', order_by='start'):
        #     # print(dir(transcript))
        #     print('transcript', transcript.id, transcript.start, transcript.stop)
        #     for exon in db.children(transcript, featuretype='exon', order_by='start'):
        #         print('exon', exon.id, exon.start, exon.stop)
        #         exon_to_transcript[exon.id].add(transcript.id)


        


        # for transcript in db.children(gene, featuretype='transcript', order_by='start'):
        #     print(transcript)
    #     print(gene)
    # for i in db.children(gene, featuretype='exon', order_by='start'):
    # for tr in db.children(gene, featuretype='transcript', order_by='start'):
    #     # print(tr.id, dir(tr)) 
    #     splice_sites = []
    #     for j in db.children(tr, featuretype='exon', order_by='start'):
    #         # print(j, j.start, j.end)
    #         splice_sites.append(j.start -1)
    #         splice_sites.append(j.end)
    #     splice_sites_tmp = splice_sites[1:-1]
    #     splice_sites = []
    #     for i in range(0, len(splice_sites_tmp),2):
    #         splice_sites.append( (splice_sites_tmp[i], splice_sites_tmp[i+1]) )
    #     # splice_sites = [item for item in zip(splice_sites[:-1], splice_sites[1:])]
    #     ref_isoforms[tr.id] = splice_sites
    # print(ref_isoforms)
    
    # print(collapsed_exon_to_transcript)