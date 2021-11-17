#!/usr/bin/python

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.colors import colorConverter
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

import networkx as nx
import community as community_louvain
from fa2 import ForceAtlas2


# Data structure to store all the information from the network raw file
#------------------------------------------------------------------------------
class GData:

    #------------------------------------------------------------------------------
    def __init__(self):

        self.atoms = None
        self.edges = None
        self.proportions = None
        self.atom_indexes = None # dictionary to recover atom indexes
        self.predecessors = None # list of predecessors for each vertices

        return

    #------------------------------------------------------------------------------
    def _get_atom_index(self, atom):

        if atom in self.atom_indexes:
            return self.atom_indexes[atom]

        else:
            self.atoms.append(atom)
            self.atom_indexes[atom] = len(self.atoms)-1
            self.edges.append([])
            self.predecessors.append([])
            self.proportions.append((0,0))

            return len(self.atoms)-1

    #------------------------------------------------------------------------------
    def n_atoms(self):
        return len(self.atoms)

    #------------------------------------------------------------------------------
    def load_from_file(self, filename, ncell_min, score_min):

        # read file content
        file = open(filename)
        content = file.read().splitlines()
        file.close()

        self.atoms = []
        self.edges = [] # for each atom, list of connected atoms with: index, score (relative distance), individual pos score, individual negative score
        self.predecessors = []
        self.proportions = []
        self.atom_indexes = {} # dictionary to recover atom indexes

        # read the graph data
        for line in content:

            tokens = line.split(' ')

            if np.min([int(tokens[2]), int(tokens[3])]) >= ncell_min:
                gene_name = tokens[0]
                gene_val = int(tokens[1])

                atom = (gene_name, gene_val)
                atom_ind = self._get_atom_index(atom)

                self.proportions[atom_ind] = (int(tokens[2]), int(tokens[3]))

                token_index = 4
                while token_index < len(tokens):

                    atom_connected = (tokens[token_index], int(tokens[token_index+1]))
                    atom_connected_index = self._get_atom_index(atom_connected)
                    atom_connected_score = float(tokens[token_index+2])
                    atom_connected_pos_score = int(tokens[token_index+3])
                    atom_connected_neg_score = int(tokens[token_index+4])

                    if atom_connected_score >= score_min:
                        #self.edges[atom_ind].append((atom_connected_index, atom_connected_score, atom_connected_pos_score, atom_connected_neg_score))
                        self.edges[atom_connected_index].append((atom_ind, atom_connected_score, atom_connected_pos_score, atom_connected_neg_score))
                        self.predecessors[atom_ind].append(atom_connected_index)
                    token_index += 5

        return


# Class graph to represent a co-expression graph with 2d positions and clusters
#------------------------------------------------------------------------------
class Graph:

    #------------------------------------------------------------------------------
    def __init__(self, name = 'unnamed'):

        self.atoms = None
        self.atom_indexes = None
        self.edges = None
        self.weights = None

        self.positions = None
        self.clusters = None
        self.atom_cluster = None # association of atom indexes to cluster indexes

        self.dendogram = None

        self.nx_graph = None

        self.color_palette = ['#000000','#004949','#009292','#ff6db6','#ffb6db', '#490092','#006ddb','#b66dff','#6db6ff','#b6dbff', '#920000','#924900','#db6d00','#24ff24','#ffff6d']
        #self.color_palette = [elt for elt in mcolors.TABLEAU_COLORS] + ['goldenrod', 'indigo', 'steelblue']

        self.name = name

        np.random.seed(42)

        return

    #------------------------------------------------------------------------------
    def create_from_gdata(self, data, exclude_mt_rp, filter_edges = False):

        # filter_edges: additional filter: exclude isolated edges

        # select only atoms that appear in at least one edge
        selected_atoms = []
        for ind in range(data.n_atoms()):
            selected = False
            gene_name = data.atoms[ind][0]
            if not exclude_mt_rp or ('RPL' not in gene_name and 'RPS' not in gene_name and 'MT' not in gene_name):
                for edge1 in data.edges[ind]:
                    # look for the reversed edge
                    for edge2 in data.edges[edge1[0]]:
                        if edge2[0] == ind:
                            if not filter_edges or (len(data.edges[ind]) > 1 or len(data.edges[edge1[0]]) > 1):
                                selected = True
                                break
                    if selected:
                        break
            if selected:
                selected_atoms.append(ind)

        self.atoms = [data.atoms[ind] for ind in selected_atoms]
        self.atom_indexes = {self.atoms[ind]:ind for ind in range(len(self.atoms))}

        # dictionaries to access between original data and reduced graph indexes for vertices
        gdata_to_graph_index = [-1 for _ in range(data.n_atoms())]
        graph_to_gdata_index = [-1 for _ in range(len(self.atoms))]
        graph_ind = 0
        for gdata_ind in selected_atoms:
            gdata_to_graph_index[gdata_ind] = graph_ind
            graph_to_gdata_index[graph_ind] = gdata_ind
            graph_ind += 1

        # compute the min and max weight to normalize weights
        min_weight = -1
        max_weight = -1
        for gdata_ind in selected_atoms:
            for edge in data.edges[gdata_ind]:
                if min_weight < 0 or edge[1] < min_weight:
                    min_weight = edge[1]
                if max_weight < 0 or edge[1] > max_weight:
                    max_weight = edge[1]

        # create a list of weighted edges of the graph: undirected -> mean weight
        self.edges = []
        self.weights = []

        added_edges = [] # contains all edges added to the graph yet
        nb = 0

        for gdata_ind in selected_atoms:
            if nb % 100 == 0:
                print(nb, ' vertices processed over ', len(selected_atoms))
            nb += 1
            for edge in data.edges[gdata_ind]:
                gdata_ind2 = edge[0]

                # make sure both vertices have been selected
                if gdata_ind2 in selected_atoms:

                    pair = None
                    if gdata_ind < gdata_ind2:
                        pair = [gdata_ind, gdata_ind2]
                    else:
                        pair = [gdata_ind2, gdata_ind]

                    # make sure this edge has not been created yet
                    if not pair in added_edges:

                        added_edges.append(pair)

                        # look for a reversed edge
                        rweight = None
                        for redge in data.edges[gdata_ind2]:
                            if redge[0] == gdata_ind:
                                rweight = redge[1]
                                break
                        weight = None
                        if rweight == None:
                            weight = edge[1]
                        else:
                            weight = (edge[1]+rweight)/2.

                        weight = (weight-min_weight)/(max_weight-min_weight)

                        # make sure the other vertex is in the graph as well
                        graph_ind2 = gdata_to_graph_index[gdata_ind2]
                        if graph_ind2 >= len(selected_atoms):
                            print('error: ', graph_ind2)

                        self.edges.append((gdata_to_graph_index[gdata_ind], graph_ind2))
                        self.weights.append(weight)

        return


    #------------------------------------------------------------------------------
    def create_from_gdata_regul(self, data, exclude_mt_rp, filter_edges = False):

        # filter_edges: additional filter: exclude isolated edges

        candidates = [True for _ in range(data.n_atoms())]
        stop = False

        pending = []
        pending2 = [ind for ind in range(data.n_atoms())]

        while len(pending2) > 0: # look for a fix point

            pending = pending2
            pending2 = []

            for ind in pending:

                excluded = False

                if candidates[ind]:

                    gene_name = data.atoms[ind][0]

                    # exclude ribosomal and mitochondrial genes
                    if exclude_mt_rp and ('RPL' in gene_name or 'RPS' in gene_name or 'MT' in gene_name):
                        excluded = True

                    # exclude genes if they are connected to only one other gene
                    if filter_edges and not excluded:
                        nb_candidates_edges = 0
                        last_candidates_edge = -1
                        for edge in data.edges[ind]:
                            ind2 = edge[0]
                            if candidates[ind2]:
                                nb_candidates_edges += 1
                                last_candidates_edge = ind2
                            if nb_candidates_edges > 1:
                                break

                        if nb_candidates_edges == 1:

                            nb_candidates_edges_bis = 0
                            last_candidates_edge_bis = -1

                            for edge in data.edges[last_candidates_edge]:
                                ind3 = edge[0]
                                if candidates[ind3]:
                                    nb_candidates_edges_bis += 1
                                    last_candidates_edge_bis = ind3
                                if nb_candidates_edges_bis > 1:
                                    break

                            if nb_candidates_edges_bis == 0:
                                excluded = True

                            elif nb_candidates_edges_bis == 1 and last_candidates_edge_bis == ind:
                                excluded = True

                    # exclude vertices if the neighbours are not connected to enough vertices
                    if not excluded:
                        nb_successors = 0
                        for edge in data.edges[ind]:
                            ind2 = edge[0]
                            if candidates[ind2]:
                                nb_successors += 1
                        if nb_successors < 2:
                            # look for the number of predecessor
                            nb_predecessors = 0
                            for ind2 in data.predecessors[ind]:
                                if candidates[ind2]:
                                    nb_predecessors += 1

                            if nb_predecessors < 2 or nb_predecessors >= 30:
                                excluded = True


                    if excluded:
                        candidates[ind] = False
                        pending2 += [edge[0] for edge in data.edges[ind] if candidates[edge[0]]]
                        pending2 += [ind2 for ind2 in data.predecessors[ind] if candidates[ind2]]

            pending2 = np.unique(pending2)

        selected_atoms = [ind for ind in range(data.n_atoms()) if candidates[ind]]
        print('N selected atoms: ', len(selected_atoms))

        #print('selected atoms: ', selected_atoms)

        self.atoms = [data.atoms[ind] for ind in selected_atoms]
        self.atom_indexes = {self.atoms[ind]:ind for ind in range(len(self.atoms))}

        # dictionaries to access between original data and reduced graph indexes for vertices
        gdata_to_graph_index = [-1 for _ in range(data.n_atoms())]
        graph_to_gdata_index = [-1 for _ in range(len(self.atoms))]
        graph_ind = 0
        for gdata_ind in selected_atoms:
            gdata_to_graph_index[gdata_ind] = graph_ind
            graph_to_gdata_index[graph_ind] = gdata_ind
            graph_ind += 1

        # compute the min and max weight to normalize weights
        min_weight = -1
        max_weight = -1
        for gdata_ind in selected_atoms:
            for edge in data.edges[gdata_ind]:
                if min_weight < 0 or edge[1] < min_weight:
                    min_weight = edge[1]
                if max_weight < 0 or edge[1] > max_weight:
                    max_weight = edge[1]

        # create a list of weighted edges of the graph: undirected -> mean weight
        self.edges = []
        self.weights = []

        nb = 0

        for gdata_ind in selected_atoms:
            if nb % 100 == 0:
                print(nb, ' over ', len(selected_atoms))
            nb += 1
            for edge in data.edges[gdata_ind]:
                gdata_ind2 = edge[0]

                # make sure both vertices have been selected
                if gdata_ind2 in selected_atoms:

                    weight = edge[1]
                    weight = (weight-min_weight)/(max_weight-min_weight)

                    # make sure the other vertex is in the graph as well
                    graph_ind2 = gdata_to_graph_index[gdata_ind2]
                    if graph_ind2 >= len(selected_atoms):
                        print('error: ', graph_ind2)

                    # self.edges.append((gdata_to_graph_index[gdata_ind], graph_ind2))
                    # inverted edges
                    self.edges.append((graph_ind2, gdata_to_graph_index[gdata_ind]))
                    self.weights.append(weight)

        return



    #------------------------------------------------------------------------------
    def plot(self, ax, col_option, arrows = False, cluster_size_limit = 20):

        # arrows: draw arrows instead of edges

        # plot the graph from force atlas 2 2d positions

        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        # plot the clusters
        for ind_cluster in self.clusters:

            if len(self.clusters[ind_cluster]) >= cluster_size_limit:

                pos_cluster = None
                for ind_atom in self.clusters[ind_cluster]:
                    if pos_cluster == None:
                        pos_cluster = self.positions[ind_atom]
                    else:
                        pos_cluster = (pos_cluster[0]+self.positions[ind_atom][0], pos_cluster[1]+self.positions[ind_atom][1])

                pos_cluster = (pos_cluster[0]/len(self.clusters[ind_cluster]), pos_cluster[1]/len(self.clusters[ind_cluster]))

                text = ax.text(pos_cluster[0], pos_cluster[1], 'c' + str(ind_cluster), ha="center", va="center", color='tan', fontsize=12, fontweight='bold', zorder=2, bbox=dict(facecolor='blue', alpha=0.5, boxstyle='round')).set_clip_on(True)

        if col_option == 'clustering_colors':
            # create a legend for the clusters
            legend_elements = [ mpatches.Patch(facecolor=self.color_palette[ind_cluster%len(self.color_palette)], label='c'+str(ind_cluster)) for ind_cluster in self.clusters if len(self.clusters[ind_cluster]) >= cluster_size_limit]
            ax.legend(handles=legend_elements, loc='upper left')

        elif col_option == '01_colors': # legend for 0/1 value
            # create for expressed/unexpressed
            legend_elements = [ mpatches.Patch(facecolor='red', label='unexpressed'), mpatches.Patch(facecolor='forestgreen', label='expressed')]
            ax.legend(handles=legend_elements, loc='upper left')

        # plot the vertices
        for ind in range(len(self.atoms)):

            atom = self.atoms[ind]

            if col_option == 'clustering_colors':
                if len(self.clusters[self.atom_cluster[ind]]) >= cluster_size_limit:
                    atom_color = self.color_palette[self.atom_cluster[ind]%len(self.color_palette)]
                else:
                    atom_color = 'grey'

            elif col_option == '01_colors': # color based on discrete value

                if atom[1] == 0:
                    atom_color = 'red'
                else:
                    atom_color = 'forestgreen'

            else:
                atom_color = 'black'

            # elif col_option == 'TF_01_colors': # color based on transcription factor
            # if atom[0] in TF_names:
            #    if atom[1] == 0:
            #        atom_color = 'red'
            #    else:
            #        atom_color = 'green'


            text = ax.text(self.positions[ind][0], self.positions[ind][1], atom[0]+'_'+str(atom[1]), ha="center", va="center", color=atom_color, fontsize=8, fontweight='bold', zorder=1).set_clip_on(True)



        if arrows:
            X = [self.positions[edge[0]][0] for edge in self.edges]
            Y = [self.positions[edge[0]][1] for edge in self.edges]
            U = [self.positions[edge[1]][0] - self.positions[edge[0]][0] for edge in self.edges]
            V = [self.positions[edge[1]][1] - self.positions[edge[0]][1] for edge in self.edges]
            ax.quiver(X, Y, U, V, scale_units='xy', angles='xy', scale=1, color=[colorConverter.to_rgba('black', alpha=self.weights[ind]) for ind in range(len(self.edges))])

        else:

            drawLines = []
            for edge in self.edges:
                p1 = self.positions[edge[0]]
                p2 = self.positions[edge[1]]
                drawLines.append([p1[0], p2[0]])
                drawLines.append([p1[1], p2[1]])
            ax.plot(*drawLines, zorder=0)

            for ind in range(len(self.edges)):
                #ax.lines[ind].set_linewidth( (nx_edges[ind][2]/9000.)*1. )
                ax.lines[ind].set_color(colorConverter.to_rgba('black', alpha=self.weights[ind]))
                ax.lines[ind].set_linewidth( self.weights[ind] )

        return

    #------------------------------------------------------------------------------
    def _create_nx_graph(self):

        # create the nx graph

        self.nx_graph = nx.Graph()
        self.nx_graph.add_nodes_from([ind for ind in range(len(self.atoms))])

        nx_edges = [(self.edges[ind][0], self.edges[ind][1], self.weights[ind]) for ind in range(len(self.edges))]
        self.nx_graph.add_weighted_edges_from(nx_edges)

        return


    #------------------------------------------------------------------------------
    def compute_positions(self):

        # compute 2d positions with Force Atlas 2 algorithm

        if self.nx_graph == None:
            self._create_nx_graph()

        forceatlas2 = ForceAtlas2(
                                # Behavior alternatives
                                outboundAttractionDistribution=True,  # Dissuade hubs
                                linLogMode=False,  # NOT IMPLEMENTED
                                adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
                                edgeWeightInfluence=1.,

                                # Performance
                                jitterTolerance=1.0,  # Tolerance
                                barnesHutOptimize=True,
                                barnesHutTheta=1.2,
                                multiThreaded=False,  # NOT IMPLEMENTED

                                # Tuning
                                scalingRatio=2.0,
                                strongGravityMode=False,
                                gravity=1.0,

                                # Log
                                verbose=True)

        self.positions = forceatlas2.forceatlas2_networkx_layout(self.nx_graph, pos=None, iterations=2000)

        return

    #------------------------------------------------------------------------------
    def compute_clusters(self, res):

        if self.nx_graph == None:
            self._create_nx_graph()

        partition = community_louvain.best_partition(self.nx_graph, resolution=res)

        # for each atom index, an index of the cluster
        self.clusters = {}

        # association of the cluster index to each atom index
        self.atom_cluster = [-1 for ind in range(len(self.atoms))]

        for atom_index in partition:

            ind_cluster = partition[atom_index]

            if ind_cluster in self.clusters:
                self.clusters[ind_cluster].append(atom_index)
            else:
                self.clusters[ind_cluster] = [atom_index]

            self.atom_cluster[atom_index] = ind_cluster

        return

    #------------------------------------------------------------------------------
    def save(self, filename):

        file = open(filename, 'w')

        # atoms
        for ind_atom in range(len(self.atoms)):
            atom = self.atoms[ind_atom]
            file.write(atom[0] + '_' + str(atom[1]))
            if ind_atom < len(self.atoms)-1:
                file.write(' ')
        file.write('\n')

        # write the positions
        for ind_atom in range(len(self.atoms)):
            pos = self.positions[ind_atom]
            file.write(str(pos[0]) + '_' + str(pos[1]))
            if ind_atom < len(self.atoms)-1:
                file.write(' ')
        file.write('\n')

        # edges (with weights)
        for ind_edge in range(len(self.edges)):
            edge = self.edges[ind_edge]
            file.write(str(edge[0]) + '_' + str(edge[1]) + '_' + str(self.weights[ind_edge]))
            if ind_edge < len(self.edges)-1:
                file.write(' ')
        file.write('\n')

        # clusters
        file.write(str(len(self.clusters)))
        for cluster_index in self.clusters:
            cluster = self.clusters[cluster_index]
            file.write('\n')
            for ind in range(len(cluster)):
                ind_atom = cluster[ind]
                file.write(str(ind_atom))
                if ind < len(cluster)-1:
                    file.write(' ')

        file.close()

        return

    #------------------------------------------------------------------------------
    def load_from_file(self, filename):

        file = open(filename, 'r')
        content = file.read().splitlines()
        file.close()

        # read the atoms
        self.atoms = []
        tokens = content[0].split(' ')
        for elt in tokens:
            elts = elt.split('_')
            self.atoms.append((elts[0], int(elts[1])))
        self.atom_indexes = {self.atoms[ind]:ind for ind in range(len(self.atoms))}

        # read the positions
        self.positions = []
        tokens = content[1].split(' ')
        for elt in tokens:
            elts = elt.split('_')
            self.positions.append((float(elts[0]), float(elts[1])))

        # read the edges
        self.edges = []
        self.weights = []
        tokens = content[2].split(' ')
        for elt in tokens:
            edge_elts = elt.split('_')
            self.edges.append((int(edge_elts[0]), int(edge_elts[1])))
            self.weights.append(float(edge_elts[2]))

        # read the clusters
        n_clusters = int(content[3])
        self.clusters = {ind:[] for ind in range(n_clusters)}
        self.atom_cluster = [-1 for _ in range(len(self.atoms))]
        for cluster_index in range(n_clusters):
            self.clusters[cluster_index] = [int(elt) for elt in content[4+cluster_index].split(' ')]
            for atom_index in self.clusters[cluster_index]:
                self.atom_cluster[atom_index] = cluster_index
        return


    #---------------------------------------------------------------------------
    def print_clusters(self):

        print('clusters from ' + self.name + ' graph:')

        for cluster in self.clusters:
            atoms = [self.atoms[ind] for ind in self.clusters[cluster]]
            disp = 'c' + str(cluster) + ' (' + str(len(atoms)) + ') : '
            for atom in atoms:
                disp += atom[0] + '_' + str(atom[1]) + ', '
            print(disp + '\n')

        return

    #------------------------------------------------------------------------------
    def save_clusters(self, filename):

        file = open(filename, 'w')

        for cluster in self.clusters:

            file.write(str(cluster) + ' ')
            for atom_ind in self.clusters[cluster]:
                atom = self.atoms[atom_ind]
                file.write(atom[0] + '_' + str(atom[1]) + ' ')
            file.write('\n')

        file.close()

        return


    #---------------------------------------------------------------------------
    @staticmethod
    def compare_clusters(graph_1, graph_2, threshold = 0):

        # compare the clusters of two different graph: count the number of atoms in each pair of clusters

        # build a comparison matrix
        shared_atoms = [ [ [ ] for _ in graph_2.clusters ] for _ in graph_1.clusters ]
        for cluster_g1 in graph_1.clusters:
            for cluster_g2 in graph_2.clusters:
                atoms_g1 = [graph_1.atoms[ind] for ind in graph_1.clusters[cluster_g1]]
                atoms_g2 = [graph_2.atoms[ind] for ind in graph_2.clusters[cluster_g2]]
                shared_atoms[cluster_g1][cluster_g2] = [atom for atom in atoms_g1 if atom in atoms_g2]

        # display atoms shared in two clusterss
        for cluster_g1 in graph_1.clusters:
            for cluster_g2 in graph_2.clusters:

                if len(shared_atoms[cluster_g1][cluster_g2]) > threshold:
                    disp = graph_1.name + ' c' + str(cluster_g1) + ' #' + str(len(graph_1.clusters[cluster_g1])) + ' vs ' + graph_2.name + ' c' + str(cluster_g2) + ' '
                    disp += ' #' + str(len(graph_2.clusters[cluster_g2])) + ' '
                    disp += '(' + str(len(shared_atoms[cluster_g1][cluster_g2])) + ') : '
                    for atom in shared_atoms[cluster_g1][cluster_g2]:
                        disp += atom[0] + '_' + str(atom[1]) + ', '
                    print(disp + '\n')

        #print(shared_atoms)

        return
