"""This script generates and analyzes overlap graphs of fasta files.
"""

# pylint: disable=invalid-name

from __future__ import with_statement, print_function, generators
import glob
from multiprocessing import Pool
import pprint
import os
import sys
import matplotlib.pyplot as plt
from pylab import rcParams
import numpy as np
import networkx as nx

from sv_pipeline import networkx_helpers as nx_helpers
from sv_pipeline import overlap
from sv_pipeline.data_io import *


## Download the compressed file from this site and extract it:
#https://xritza01.u.hpc.mssm.edu/trios/2016-05-12-data-for-networks/GR38/NA19240/dels/

def draw_community_bar_chart(graph):
    """draw bar chart of size of communities"""
    communities = nx_helpers.get_communities(graph)
    len_communities = sorted([len(c) for c in communities])
    plt.bar(np.arange(len(len_communities)), len_communities)


def read_paf(prefix, fasta_filename, min_matching_length, should_filter_paf=True):
    """ Reads the PAF file from minimap.  If should_filter_paf=True,
    filters the file according to Ali's specs.
    The returned list has three values.
    The zerotih value is read0, first value is read1
    and the third value is that quality of their mapping.
    See minimap for more info on quality rating.
    """

    def pass_filters(line):
        """returns True if the line passes Ali's filters, False otherwise.
        From filename of previous file:
        all_hg004.mapq200.alen3000.ovtol500
        """
        map_cutoff = 200
        alen_cutoff = 3000
        ov_tol = 500

        kept = False
        qname, ql, qs, qe, strand, tname, tl, ts, te, _, alen, mapq = line[0:12]
        qs, qe, ql = int(qs), int(qe), int(ql)
        ts, te, tl = int(ts), int(te), int(tl)
        mapq = int(mapq)
        alen = int(alen)
        ov = overlap.overlap(qname, tname, mapq, -1, "+", qs, qe, ql, strand, ts, te, tl)
        if mapq > map_cutoff and alen > alen_cutoff and ov.hasFullOverlap(ov_tol):
            #print l.strip()
            kept = True
        else:
            pass

        return kept

    # TODO put io parts in data_io.py

    paf_filename = temp_dir + "/tmp_" + prefix + ".paf"
    minimap_command = "./minimap -Sw5 -L{} -m0 {} {} > {}"\
                      .format(min_matching_length, fasta_filename, fasta_filename, paf_filename)
    os.system(minimap_command)

    # populate alignedreads by reading in paf (filter if neccessary)
    alignedreads = []
    with open(paf_filename) as fin:
        for line in fin:
            row = line.strip().split()
            if not should_filter_paf or pass_filters(row):
                alignedreads.append([row[0], row[5], int(row[10])])
    return alignedreads

def node_community_colors(graph, communities):
    """runs a community detection algorithm on graph and
    returns a coloring of the nodes based on the found communities
    """
    colors = nx_helpers.generate_colors(len(communities))

    def which_color(node):
        """finds which community node is in and returns
        its corresponding color
        """
        for i, com in enumerate(communities):
            if node in com:
                return colors[i]
        return nx_helpers.rgb_to_hex((0, 0, 0))

    node_colors = [which_color(node) for node in graph.nodes()]
    return node_colors

def node_set_colors(nodes, spanset, gapset, preset, postset):
    """returns a list of colors for coloring nodes based
    on which set each node is in"""

    node_colors = []
    for n in nodes:
        if n in preset:
            node_colors.append(nx_helpers.rgb_to_hex((255, 0, 0)))
        elif n in postset:
            node_colors.append(nx_helpers.rgb_to_hex((255, 255, 0)))
        elif n in gapset:
            node_colors.append(nx_helpers.rgb_to_hex((0, 10, 250)))
        elif n in spanset:
            node_colors.append(nx_helpers.rgb_to_hex((0, 250, 10)))
        else:
            # uncategorized
            node_colors.append(nx_helpers.rgb_to_hex((0, 0, 0)))
    return node_colors

def generate_graph(prefix, fasta_filename, min_matching_length):
    """generates and returns overlap graph from .paf files
    """
    alignedreads = read_paf(prefix, fasta_filename, min_matching_length, should_filter_paf=True)
    aligned = [(t, h) for t, h, _ in alignedreads]
    graph = nx.Graph()
    graph.add_edges_from(aligned)
    return graph

def drop_small_communities(graph, communities, n=4):
    """removes nodes from graph in they are in communities smaller than n
    """
    for community in communities:
        if len(community) < n:
            nx_helpers.remove_nodes(graph, community)
    communities = [c for c in communities if len(c) >= n]
    return graph, communities

def mapping_quality(graph, spanset, gapset):
    """Determines the quality of the mapping (assignment of edges)
    based on the "ground truth" of spanset and gapset.
    Sums up number of edges between spanset and gapset.
    Assumes undirected graph - see comments"""
    the_sum = sum(sum(1 for edge in graph.edges(node) if edge[1] in gapset) for node in spanset)
    # if directed graph, uncomment this:
    #the_sum += sum(sum(1 for edge in graph.edges(node) if edge[1] in spanset) for node in gapset)
    return the_sum

def community_quality(communities, spanset, gapset):
    """Determines the quality of the communities based
    on the "ground truth" of spanset and gapset.
    First, determines which community corresponds to gapset and spanset.
    Then, returns number of wrong nodes.
    """
    if len(communities) != 2:
        #print("not two communities")
        return -1

    com_sets = [set(c) for c in communities]
    spanset = set(spanset)
    gapset = set(gapset)

    spanset_0 = len(com_sets[0].difference(spanset))
    spanset_1 = len(com_sets[1].difference(spanset))
    gapset_0 = len(com_sets[0].difference(gapset))
    gapset_1 = len(com_sets[1].difference(gapset))

    # used for determining which community corresponds to gapset and spanset
    spanset_i = 1 - np.argmax([spanset_0, spanset_1])
    gapset_i = 1 - np.argmax([gapset_0, gapset_1])

    if spanset_i == gapset_i:
        #print("Error in finding community quality")
        return -1
    elif spanset_i == 0:
        return spanset_0 + gapset_1
    elif spanset_i == 1:
        return spanset_1 + gapset_0
    else:
        #print("Unexpected condition in finding community quality")
        return -1

def make_line_plot(bed_filename, the_sets):
    """Makes an IGV-style drawing with colors"""

    # TODO separate out this data grabbing into data_io, then make line plot in this function
    import csv

    # Pull coords from bed file
    coords = None # {name: (start, end)} for each read
    with open(bed_filename, 'r') as csvfile:
        bed_reader = csv.reader(csvfile, delimiter='\t')
        coords = {row[3]: (int(row[1]), int(row[2])) for row in bed_reader}

    # normalize values to between 0 and 1
    min_left_coord = min(x for x, _ in coords.values())
    the_range = float(max(y for _, y in coords.values()) - min_left_coord)
    coords = {read: (float(coord[0] - min_left_coord) / the_range,\
                     float(coord[1] - min_left_coord) / the_range)\
                     for read, coord in coords.items()}

    gapset, spanset, preset, postset = the_sets
    colors = node_set_colors(coords.keys(), gapset, spanset, preset, postset)

    y_increment = (1. / float(len(coords)))
    y_values = [float(i) * y_increment for i in range(0, len(coords))]
    for i, (coord, y_value) in enumerate(zip(coords.values(), y_values)):
        plt.plot(list(coord), [y_value, y_value], color=colors[i], linestyle='-', linewidth=1.5)
    plt.axis('off')
    plt.title("IGV style line plot")

def make_four_pdf(args):
    """
    Generates four graphs for the
    structural variant defined by merged_filename
    """
    merged_filename = args[0]
    the_dir = args[1]

    # if there are fewer than threshold reads then skip it
    threshold = 25 # threshold before plotting.
    if len(open(merged_filename).readlines()) < threshold:
        print('skipping %s' % (merged_filename))
        return

    rcParams['figure.figsize'] = 30, 30
    plt.clf()
    plt.figure(1)

    prefix = merged_filename[len(the_dir):-11]
    bed_filename = the_dir + prefix + '-refcoords.bed'
    fasta_filename = the_dir + prefix + ".fa"

    # Checks size of variant
    remove_punctuation = lambda x: ''.join(e for e in x if e.isdigit() or e == '.')
    coords = [int(remove_punctuation(a)) for a in prefix.split('_')[1:3]]
    dist = coords[1] - coords[0]
    if dist < 100:
        print('skipping %s' % (merged_filename))
        return




    min_matching_length = 600 # hard-code for now.
    graph = generate_graph(prefix, fasta_filename, min_matching_length)
    spanset, gapset, preset, postset = get_read_classifications(prefix,\
                                            bed_filename, merged_filename)

    # Draw Ground Truth
    plt.subplot(2, 2, 1)
    node_colors = node_set_colors(graph.nodes(), spanset, gapset, preset, postset)
    pos = nx.spring_layout(graph)

    assert(len(node_colors) == nx.number_of_nodes(graph))
    title = "Chr {0}; L={1}; Ground Truth Colors\n\
            Red=Preset, Yellow=Postset, Blue=GapSet, Green=SpanSet"\
            .format(prefix, min_matching_length)
    nx.draw(graph, node_color=node_colors, node_size=100, pos=pos)
    plt.title(title)

    # Draw Ground Truth with squashed nodes
    plt.subplot(2, 2, 2)
    # squash preset and postset nodes
    graph = nx_helpers.remove_nodes(graph, preset)
    graph = nx_helpers.remove_nodes(graph, postset)
    node_colors = node_set_colors(graph.nodes(), spanset, gapset, preset, postset)
    #pos = nx.spring_layout(graph)
    assert(len(node_colors) == nx.number_of_nodes(graph))
    title = "Chr {0}; L={1}; Ground Truth Colors \n\
            Removed Preset and Postsetnodes; Blue=GapSet, Green=SpanSet"\
            .format(prefix, min_matching_length)
    nx.draw(graph, node_color=node_colors, node_size=100, pos=pos)
    plt.title(title)

    # Drop Small Communities and Draw
    plt.subplot(2, 2, 3)
    communities = nx_helpers.get_communities(graph)
    graph, communities = drop_small_communities(graph, communities)
    node_colors = node_community_colors(graph, communities)
    assert(len(node_colors) == nx.number_of_nodes(graph))
    title = "Chr {0}; L={1}; After Removing Small Communities; NumCom={2}\n\
            ComQual={3}, MapQual={4}"\
            .format(prefix, min_matching_length, len(communities),\
                    community_quality(communities, spanset, gapset),\
                    mapping_quality(graph, spanset, gapset))
    nx.draw(graph, node_color=node_colors, node_size=100, pos=pos)
    plt.title(title)

    # IGV Line Plot
    plt.subplot(2, 2, 4)
    make_line_plot(bed_filename, (spanset, gapset, preset, postset))

    plt.savefig('figs/%s-communities.pdf' % (prefix))

def four_graphs(the_dir):
    """
    Generates four graphs for each structural variant in the directory
    formerly
    """

    files = get_files(the_dir)
    print('Looking in directory %s*merged.txt' % (the_dir))
    print('There are %d files' % (len(files)))
    the_dirs = [the_dir for _ in files]
    zipped = zip(files, the_dirs)
    p = Pool()
    p.map(make_four_pdf, zipped)

def sixteen_graphs(the_dir):
    """ generates graphs for each structual variant
    """
    rcParams['figure.figsize'] = 30, 30
    plt.clf()
    plt.figure(1)

    # should look like: read_data/all_files/chr4_124,017,492_124,029,032_merged.txt
    merged_files = glob.glob(the_dir + '*merged.txt')
    print("Running for {} regions".format(len(merged_files)))
    for merged_filename in merged_files:
        # get filenames
        prefix = merged_filename[len(the_dir):-11]
        fasta_filename = the_dir + prefix + ".fa"
        bed_filename = the_dir + prefix + "-refcoords.bed"
        print('Using ' + prefix)

        for min_matching_length in range(100, 1700, 100):
            print(min_matching_length)
            # used for ground truth
            spanset, gapset, preset, postset = get_read_classifications(prefix,\
                                                bed_filename, merged_filename)
            # Generate and prune graph
            graph = generate_graph(prefix, fasta_filename, min_matching_length)
            graph = nx_helpers.remove_nodes(graph, preset)
            graph = nx_helpers.remove_nodes(graph, postset)

            # Plot the graph
            plt.subplot(4, 4, min_matching_length/100)
            communities = nx_helpers.get_communities(graph)
            graph, communities = drop_small_communities(graph, communities)
            node_colors = node_community_colors(graph, communities)
            pos = nx.spring_layout(graph)
            title = "Chr {0};\n L={1}; NumCom={2}\nComQual = {3}, MapQual={4}"\
                    .format(prefix, min_matching_length, len(communities),\
                            community_quality(communities, spanset, gapset),\
                            mapping_quality(graph, spanset, gapset))
            nx.draw(graph, node_color=node_colors, node_size=100, pos=pos)
            plt.title(title)
        plt.savefig("figs/" + prefix + '-16-communities.pdf')
        plt.clf()


