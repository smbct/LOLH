#!/usr/bin/env python3

import network

#------------------------------------------------------------------------------
#                Write the gene clusters into a file
#------------------------------------------------------------------------------
def write_clusters(filename, graph, clusters):
    
    file = open(filename, 'w')
    for cluster in clusters:
        for atom_ind in cluster:
            atom = graph.atoms[atom_ind]
            file.write(atom[0]+'_'+str(atom[1])+ ' ')
        file.write('\n')
    file.close()
    
    return

#------------------------------------------------------------------------------
#                perform the clustering of a coexpression network
#------------------------------------------------------------------------------
def network_clustering(input_file, n_cell_min, score_min, louvain_param, output_file):
    
    data = network.GData()
    data.load_from_file(input_file, n_cell_min, score_min)
    
    graph = network.Graph('coexpression network')
    
    # build graph from raw data, exclude mitochondrial and ribosomal genes
    exclude_mt_rp = True
    filter_edges = True
    graph.create_from_gdata(data, exclude_mt_rp, filter_edges)
    
    graph.compute_clusters(louvain_param)
    
    # discard clusters with less than 11 genes
    n_genes_min = 10
    clusters = [graph.clusters[elt] for elt in graph.clusters if len(graph.clusters[elt]) >= n_genes_min]
    
    print('\n')
    print(len(clusters), ' clusters found')
    
    # print(clusters)
    
    print('\n')
    
    ind_cluster = 0
    for cluster in clusters:
        print('cluster ', ind_cluster, ': ', len(cluster), ' genes')
        ind_cluster += 1
    
    write_clusters(output_file, graph, clusters)