

import networkx as nx


def calc_last2reach(gene_graph, path_cover, topological_sort):
    # compute index
    index = {}
    for i in range(len(path_cover)):
        for pos, v in enumerate(path_cover[i]):
            index[(v,i)] = pos
    print('INDEX:',  index)

    # initialize last2reach
    last2reach =  { (v,i) : -1 for i in range(len(path_cover)) for v in topological_sort}

    for v in topological_sort:
        for i in range(len(path_cover)):
            if (v,i) in index:
               last2reach[(v,i)] =  index[(v,i)]
            else:
                last2reach[(v,i)] = max([last2reach[(u,i)] for u in gene_graph.predecessors(v)] )

    for i in range(len(path_cover)):
        print('PATH:', path_cover[i])
        tmp_nodes = [ ((v,j), last2reach[(v,j)]) for (v,j) in last2reach if i==j]
        print('LAST2REACH for path:', tmp_nodes)
    sys.exit()
    print('LAST2REACH:', last2reach)
    # gene_graph_compl = gene_graph.reverse()
    # topological_sort_compl = list(nx.topological_sort(gene_graph_compl))
    # for v in topological_sort_compl:
    #     last2reach[(v,i)]
    #     print(v)

