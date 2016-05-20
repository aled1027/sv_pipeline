"""This script generates and analyzes overlap graphs of fasta files.
"""
from __future__ import with_statement, print_function, generators
import argparse
import glob
import os
import sys

import matplotlib.pyplot as plt
from pylab import rcParams
import numpy as np

import networkx as nx
import networkx_helpers as nx_helpers
import overlap

parser = argparse.ArgumentParser()
parser.add_argument("type", help="Choose type of computation")
parser.add_argument("dir", help="Choose directory")

## Download the compressed file from this site and extract it:
#https://xritza01.u.hpc.mssm.edu/trios/2016-05-12-data-for-networks/GR38/NA19240/dels/

def draw_community_bar_chart(graph):
    """draw bar chart of size of communities"""
    communities = nx_helpers.get_communities(graph)
    len_communities = sorted([len(c) for c in communities])
    plt.bar(np.arange(len(len_communities)), len_communities)

def get_ref_coords(bed_filename):
    """Given a bed_filename, looks up the refcoods.bed file
    and parses it."""

    lines = []
    with open(bed_filename) as fin:
        for line in fin:
            row = line.split()
            #row[3] = row[3][1:]
            for i in [1, 2, 4]:
                row[i] = int(row[i])
            lines.append(row)
    sorted(lines, key=lambda x: x[1])
    return lines

def get_overlaps(lines):
    """Given a set of lines from get_ref_coords (a bed file),
    gets the list of overlaps AND the coordinates of the
    leftmost position in the refernce.
    """

    def get_overlap(interval1, interval2):
        """Given two intervalse (numerical lists of length 2), return
        the amount of overlap between them."""
        if interval1[1] < interval2[0] or interval2[1] < interval1[0]:
            return 0
        li = interval1 + interval2
        sorted(li)
        return max(li[2]-li[1], li[1]-li[2])

    overlaps = []
    readleftcoords = {}
    for i, _ in enumerate(lines):
        readleftcoords[lines[i][3]] = int(lines[i][1])
        for j in range(i+1, len(lines)):
            ov = get_overlap(lines[i][1:3], lines[j][1:3])
            overlaps.append([lines[i][3], lines[j][3], ov])
    return overlaps, readleftcoords

def read_paf(fasta_filename, min_matching_length, should_filter_paf=True):
    """ Reads the PAF file from minimap.  If should_filter_paf=True,
    filters the file according to Ali's specs.
    The returned list has three values.
    The zerotih value is read0, first value is read1
    and the third value is that quality of their mapping.
    See minimap for more info on quality rating.
    """

    paf_filename = "tmp_paf_file.paf"
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

def get_read_lengths(prefix):
    """Gets the lengths of the reads by pulling from the original .fa file.
    """
    # TODO
    raise RuntimeError # deprecated --- needs updating
    fafile = 'read_data/original_reads/%s.fa' % (prefix)
    lengths = {}
    name = ''
    with open(fafile) as fin:
        for line in fin:
            if line[0] == '>':
                name = line.strip()[1:]
            else:
                lengths[name] = len(line.strip())
    return lengths

def get_read_classifications(prefix, bed_filename, merged_filename):
    """Reads read classifications from disk and returns them.
    Used for "ground truth" of graph.
    """
    refset = set()
    altset = set()
    preset = set()
    postset = set()

    ## get coordinates
    #coords = [int(a) for a in prefix.split('_')[-4:-2]]
    remove_punctuation = lambda x : ''.join(e for e in x if e.isdigit() or e == '.')
    coords = [int(remove_punctuation(a)) for a in prefix.split('_')[1:]]
    #buff = int(prefix.split('_')[-1])
    midpoint = sum(coords) / 2.0
    bed_lines = get_ref_coords(bed_filename) # get the ref coordinates from BED file
    _, left_coords = get_overlaps(bed_lines)

    with open(merged_filename, 'r') as fin:
        for line in fin:
            row = line.strip().split()
            read = row[0]
            if row[2] == 'b': # BOTH reference and alternate genomes
                if left_coords[read] < midpoint: # Pre node
                    preset.add(read)
                else: # post node
                    postset.add(read)
            elif row[2] == 'a': # ALTERNATE genome only
                altset.add(read)
            elif row[2] == 'r': # REFERENCE genome only
                refset.add(read)
            else:
                print('ERROR:', row)
                sys.exit()
    return refset, altset, preset, postset

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

def node_set_colors(graph, spanset, gapset, preset, postset):
    """returns a list of colors for coloring nodes based
    on which set each node is in"""

    node_colors = []
    for n in graph.nodes():
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

def find_shortest_paths(graph, out_filename, sources, targets, k_paths):
    """ Use pathlinker to find shortest paths """
    # TODO
    raise RuntimeError # this is probably broken
    assert(k_paths > 0)
    edgelist_filename = out_filename + "edgelist.temp"
    srctgt_filename = out_filename + "srctgt.temp"
    nx.write_edgelist(graph, edgelist_filename)

    with open(srctgt_filename, 'w') as f:
        for node in graph.nodes():
            if node in sources:
                f.write(str(node) + '\tsource\n')
            elif node in targets:
                f.write(str(node) + '\ttarget\n')

    s = "python PathLinker/PathLinker.py {} {} -o {} --write-paths --k-param={}"\
            .format(edgelist_filename, srctgt_filename, out_filename, k_paths)
    try:
        os.system(s)
        return read_paths(out_filename + "k_100-paths.txt")
    except Exception as e:
        print(e)
        return []

def read_paths(filename):
    """Read data that pathlinker wrote.
    TODO update doc
    Returns a bunch of subgraphs that consist of a single path
    """
    from itertools import tee, izip
    def pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = tee(iterable)
        next(b, None)
        return izip(a, b)

    paths = []
    with open(filename) as f:
        for line in f.readlines():
            # example line:
            # 92    0.00000e+00 read_43|read_45|read_46|read_71|read_73|read_64|read_96
            path = line.strip().split("\t")[2].split("|")
            path_graph = nx.Graph()
            for p0, p1 in pairwise(path):
                path_graph.add_edge(p0, p1)
            paths.append(path_graph)
    return paths[1:]

def generate_graph(fasta_filename, min_matching_length):
    """generates and returns overlap graph from .paf files
    """
    alignedreads = read_paf(fasta_filename, min_matching_length, should_filter_paf=True)
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

def mapping_quality(graph, _prefix):
    """Determines the quality of the mapping (assignment of edges)
    based on the "ground truth" of spanset and gapset.
    Sums up number of edges between spanset and gapset.
    Assumes undirected graph - see comments"""
    spanset, gapset, _, _ = get_read_classifications(_prefix)
    the_sum = sum(sum(1 for edge in graph.edges(node) if edge[1] in gapset) for node in spanset)
    # if directed graph, uncommen this:
    #the_sum += sum(sum(1 for edge in graph.edges(node) if edge[1] in spanset) for node in gapset)
    return the_sum

def community_quality(graph, _prefix):
    """Determines the quality of the communities based
    on the "ground truth" of spanset and gapset.
    First, determines which community corresponds to gapset and spanset.
    Then, returns number of wrong nodes.
    """
    communities = nx_helpers.get_communities(graph)
    spanset, gapset, _, _ = get_read_classifications(_prefix)

    if len(communities) != 2:
        print("not two communities")
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
        print("Error in finding community quality")
        return -1
    elif spanset_i == 0:
        return spanset_0 + gapset_1
    elif spanset_i == 1:
        return spanset_1 + gapset_0
    else:
        print("Unexpected condition in finding community quality")
        return -1

def four_graphs(_dir):
    """
    Generates four graphs for each structural variant in the directory
    formerly
    """
    print('Looking in directory %s/*merged.txt' % (_dir))
    files = glob.glob(_dir + '*merged.txt')
    print('There are %d files' % (len(files)))
    for infile in files:

        # if there are fewer than threshold reads,
        # then skip it.
        threshold = 25 # threshold before plotting.
        if len(open(infile).readlines()) < threshold:
            print('skipping %s' % (infile))
            continue

        _prefix = infile[:-11]
        print(_prefix)

        rcParams['figure.figsize'] = 30, 30
        plt.clf()
        plt.figure(1)

        min_matching_length = 100 # hard-code for now.

        graph = generate_graph(_prefix, min_matching_length)


        # squash preset and postset nodes
        # used for ground truth
        spanset, gapset, preset, postset = get_read_classifications(_prefix)

        # Draw Ground Truth
        plt.subplot(2, 2, 1)

        ### DRAW NETWORK ###
        node_colors = node_set_colors(graph, spanset, gapset, preset, postset)
        pos = nx.spring_layout(graph)
        assert(len(node_colors) == len(graph.nodes()))
        title = "Chr {0}; L={1}; Ground Truth Colors\n\
                Red=Preset, Yellow=Postset, Blue=GapSet, Green=SpanSet"\
                .format(_prefix, min_matching_length)
        nx.draw(graph, node_color=node_colors, node_size=100, pos=pos)
       # nx.draw_networkx_labels(graph, pos, font_size=6)
        plt.title(title)

        # squash preset and postset nodes
        ## used for ground truth
        graph = nx_helpers.remove_nodes(graph, preset)
        graph = nx_helpers.remove_nodes(graph, postset)

        # Draw Ground Truth with squashed nodes
        plt.subplot(2, 2, 2)

        ### DRAW NETWORK ###
        node_colors = node_set_colors(graph, spanset, gapset, preset, postset)
        #pos = nx.spring_layout(graph)
        assert(len(node_colors) == len(graph.nodes()))
        title = "Chr {0}; L={1}; Ground Truth Colors \n\
                Removed Preset and Postsetnodes; Blue=GapSet, Green=SpanSet"\
                .format(_prefix, min_matching_length)
        nx.draw(graph, node_color=node_colors, node_size=100, pos=pos)
       # nx.draw_networkx_labels(graph, pos, font_size=6)
        plt.title(title)

        # Draw Communities
        plt.subplot(2, 2, 3)

        communities = nx_helpers.get_communities(graph)
        node_colors = node_community_colors(graph, communities)
        assert(len(node_colors) == len(graph.nodes()))
        title = "Chr {0}; L={1}; NumCom={2}, \nComQual={3}; MapQual={4}"\
            .format(_prefix, min_matching_length, len(communities),\
                    community_quality(graph, _prefix), mapping_quality(graph, _prefix))
        nx.draw(graph, node_color=node_colors, node_size=100, pos=pos)
        #nx.draw_networkx_labels(graph, pos, font_size=6)
        plt.title(title)


        # Drop Small Communities and Draw
        plt.subplot(2, 2, 4)
        graph, communities = drop_small_communities(graph, communities)
        node_colors = node_community_colors(graph, communities)
        #pos = nx.spring_layout(graph)
        assert(len(node_colors) == len(graph.nodes()))
        title = "Chr {0}; L={1}; After Removing Small Communities; NumCom={2}, \n\
                ComQual={3}; MapQual={4}"\
                .format(_prefix, min_matching_length, len(communities),\
                    community_quality(graph, _prefix), mapping_quality(graph, _prefix))
        nx.draw(graph, node_color=node_colors, node_size=100, pos=pos)
        #nx.draw_networkx_labels(graph, pos, font_size=6)
        plt.title(title)

        print("saving fig %s-communities.pdf" % (_prefix))
        plt.savefig('%s-communities.pdf' % (_prefix))

def sixteen_graphs(_dir):
    """ generates graphs for each structual variant
    """
    rcParams['figure.figsize'] = 30, 30
    plt.clf()
    plt.figure(1)

    # should look like: read_data/all_files/chr4_124,017,492_124,029,032_merged.txt
    merged_files = glob.glob(_dir + '*merged.txt')
    print("Running for {} regions".format(len(merged_files)))
    for merged_filename in merged_files:
        # get filenames
        prefix = merged_filename[len(_dir):-11]
        fasta_filename = _dir + prefix + ".fa"
        bed_filename = _dir + prefix + "-refcoords.bed"
        print('Using ' + prefix)

        for min_matching_length in range(100, 1700, 100):
            print(min_matching_length)
            # used for ground truth
            _, _, preset, postset = get_read_classifications(prefix, bed_filename, merged_filename)

            # Generate and prune graph
            graph = generate_graph(fasta_filename, min_matching_length)
            graph = nx_helpers.remove_nodes(graph, preset)
            graph = nx_helpers.remove_nodes(graph, postset)

            # Plot the graph
            plt.subplot(4, 4, min_matching_length/100)
            communities = nx_helpers.get_communities(graph)
            graph, communities = drop_small_communities(graph, communities)
            node_colors = node_community_colors(graph, communities)
            pos = nx.spring_layout(graph)
            title = "Chr {0};\n L={1}; NumCom={2}"\
                    .format(prefix, min_matching_length, len(communities))
            nx.draw(graph, node_color=node_colors, node_size=100, pos=pos)
            plt.title(title)
        plt.savefig(_dir + prefix + '-communities.pdf')
        plt.clf()

if __name__ == '__main__':
    args = parser.parse_args()
    DIR = args.dir
    print('DIR = %s' % DIR)
    if DIR == '':
        sys.exit('ERROR: DIRECTORY with SV calls required.')
    if args.type == "four":
        four_graphs(DIR)
    elif args.type == "sixteen":
        sixteen_graphs(DIR)
    else:
        print("Unknown command")

