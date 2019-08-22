
from collections import defaultdict
import networkx as nx

from modules.help_functions import readfq

def parse_tsv(): 
    pass


def create_graph(db): 
    # print(dir(db))
    gene_graphs = {} # gene_id : { (exon_start, exon_stop) : set() }
    collapsed_exon_to_transcript = {}
    for gene in db.features_of_type('gene'):
        # print(dir(gene))
        # print(gene.id, gene.seqid, gene.start, gene.stop, gene.attributes)
        gene_graph = nx.DiGraph(chr=str(gene.seqid))
        print( gene_graph.graph)
        collapsed_exon_to_transcript[gene.id] = defaultdict(set)
        already_parsed_exons = set()
        
        #add nodes
        for exon in db.children(gene, featuretype='exon', order_by='start'):
            collapsed_exon_to_transcript[gene.id][ (exon.start, exon.stop) ].update([ transcript_tmp for transcript_tmp in  exon.attributes['transcript_id']])
            # if (exon.start, exon.stop) in already_parsed_exons:
                
            gene_graph.add_node( (exon.start, exon.stop) )
            # print(gene_graph.nodes[(exon.start, exon.stop)])

        #add edges
        for transcript in db.children(gene, featuretype='transcript', order_by='start'):
            # print(dir(transcript))
            consecutive_exons = [exon for exon in db.children(transcript, featuretype='exon', order_by='start')]
            print('transcript', transcript.id, transcript.start, transcript.stop, [ (exon.start, exon.stop) for exon in db.children(transcript, featuretype='exon', order_by='start')])

            for e1,e2 in zip(consecutive_exons[:-1], consecutive_exons[1:]):
                # print('exon', exon.id, exon.start, exon.stop)
                gene_graph.add_edge( (e1.start, e1.stop),  (e2.start, e2.stop), weight=1  )
                

        # print(gene_graph.edges())
        gene_graphs[gene.id] = gene_graph
        
    # print(collapsed_exon_to_transcript)
    return gene_graphs, collapsed_exon_to_transcript


def get_sequences_from_choordinates(gene_graphs, ref):
    refs = {acc : seq for acc, (seq, _) in readfq(open(ref,"r"))}
    segments = {}
    for gene_id in gene_graphs:
        gene_graph = gene_graphs[gene_id]
        chromosome = gene_graph.graph['chr']
        for node in gene_graph:
            start,stop = node[0], node[1]
            seq = refs[chromosome][start : stop]

            segments[node] = seq
    # print(segments)
    return segments

# def create_global_source_sink(gene_graphs, topological_sorts):
#     for gene_id in gene_graphs:

def create_global_source_sink(gene_graphs):
    topological_sorts = {} 

    for gene_id in gene_graphs:
        gene_graph = gene_graphs[gene_id]
        in_deg = gene_graph.in_degree
        sources = [node for node, in_degree in in_deg if in_degree == 0]
        out_deg = gene_graph.out_degree
        sinks = [node for node, out_degree in out_deg if out_degree == 0]
        
        print("sources", sources)
        print("sinks", sinks)
        source_node =  ("source", "source") # follow the tuple format of other nodes
        gene_graph.add_node( source_node )
        sink_node =  ("sink", "sink") # follow the tuple format of other nodes
        gene_graph.add_node( sink_node )

        # nr_nodes = len(list(gene_graph.nodes()))
        gene_graph.add_edges_from([(source_node, s) for s in sources], weight= 1)
        gene_graph.add_edges_from([(s, sink_node) for s in sinks], weight = 1)
        print(gene_graph.edges(data=True))
        top_sort = list(nx.topological_sort(gene_graph))
        topological_sorts[gene_id] = top_sort
    return topological_sorts


def derive_path_cover(gene_graphs, topological_sorts):
    for gene_id in gene_graphs:
        top_sort = topological_sorts[gene_id]
        gene_graph = gene_graphs[gene_id]
        # print([ [G[n][nbr]["weight"] for nbr in G.neighbors(n)] for n in order ])
        # print(order)
        print("GRAPH", gene_graph.nodes(data=True))

        path_cover = []
        while True:
            longest_path = nx.algorithms.dag.dag_longest_path(gene_graph, weight='weight')
            longest_path_edges = [(n1, n2) for n1,n2 in zip(longest_path[:-1], longest_path[1:])]
            # print(longest_path)
            # print([(n1, n2) for n1,n2 in zip(longest_path[:-1], longest_path[1:])])
            tot_weight = sum( [gene_graph[n1][n2]["weight"] for n1,n2 in longest_path_edges])
            print( "Longest path", longest_path, "tot weight", tot_weight)
            if tot_weight  == 2:
                break
            else:
                path_cover.append(longest_path)
                for n1,n2 in longest_path_edges:
                    if gene_graph.in_degree[n1] != 0 and  gene_graph.out_degree[n2] != 0:
                        gene_graph[n1][n2]["weight"] = 0

        print("Nr paths:", len(path_cover))
    sys.exit()

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