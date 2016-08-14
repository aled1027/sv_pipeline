"""This script generates and analyzes overlap graphs of fasta files.
"""

# pylint: disable=invalid-name

from __future__ import with_statement, print_function, generators
import collections
import glob
import itertools
from multiprocessing import Pool
import os
import time
import warnings
import pdb

import matplotlib.pyplot as plt
import pylab as plb
import numpy as np
import networkx as nx

from sv_pipeline import networkx_helpers as nx_helpers
from sv_pipeline.data_io import *
from sv_pipeline import smith_waterman
from sv_pipeline import utils # Timer

## Download the compressed file from this site and extract it:
#https://xritza01.u.hpc.mssm.edu/trios/2016-05-12-data-for-networks/GR38/NA19240/dels/

def smith_waterman_filter(graph, params):
    fasta_filename = params['fasta_filename']
    paf_filename = params['paf_filename']
    score_threshold = params['gap_score_threshold']
    fasta_dict = get_fasta_dict(fasta_filename)

    # Parse PAF into a dictionary
    paf_dict = {}
    with open(paf_filename) as fin:
        for line in fin:
            row = line.strip().split()
            query_name, _, query_start, query_end, _,\
                    target_name, _, target_start, target_end, _, _, _ = row[0:12]

            query_start = int(query_start)
            query_end = int(query_end)
            target_start = int(target_start)
            target_end = int(target_end)
            query_seq = fasta_dict[query_name][query_start:query_end]
            target_seq = fasta_dict[target_name][target_start:target_end]

            paf_dict[query_name+target_name] = {
                    'query_name': query_name,
                    'query_start': query_start,
                    'query_end': query_end,
                    'target_name': target_name,
                    'target_start': target_start,
                    'target_end': target_end,
            }

    # Generate scores dictionary
    scores = {}
    num_good_scores = 0
    num_bad_scores = 0
    edges_to_remove = []
    for query, target in nx.edges(graph):

        # 1. Get overlap info from the paf dictionary
        if str(query + target) in paf_dict:
            # get the info
            overlap_info = paf_dict[query+target]
        elif str(target + query) in paf_dict:
            # get info and swap them
            overlap_info = paf_dict[target+query]
            query, target = target, query
        else:
            overlap_info = None

        query_start = overlap_info['query_start']
        query_end = overlap_info['query_end']
        target_start = overlap_info['target_start']
        target_end = overlap_info['target_end']

        query_seq = fasta_dict[query][query_start:query_end]
        target_seq = fasta_dict[target][target_start:target_end]

        # 2. Align the sequences using the rolling method
        bad_score = False
        min_len = min(len(query_seq), len(target_seq))

        # This for loop is slow
        # Loop over windows of the subsequence and align
        window_size = params['window_size']
        # generate scores dictionary
        cur_scores = []

        for start, end in utils.pairwise(range(0, min_len, window_size)):
            qs = query_seq[start:end]
            ts = target_seq[start:end]
            score = smith_waterman.smith_waterman(qs, ts)
            cur_scores.append(score)

        scores[str(query + target)] = cur_scores

        # analyze scores
        min_score = min(cur_scores)
        if min_score < score_threshold:
            num_bad_scores += 1
            edges_to_remove.append((query, target))
        else:
            num_good_scores += 1

    # remove edges and isolated nodes
    graph.remove_edges_from(edges_to_remove)
    isolates = list(nx.isolates(graph))
    graph.remove_nodes_from(isolates)

    # the histogram of the data
    all_scores = list(utils.flatten(list(scores.values())))
    plt.hist(all_scores)
    plt.title("histogram of num_gaps / len(aligned_sequence)\n{} bad_scores {} good_scores\nthreshold = {}"
              .format(num_bad_scores, num_good_scores, score_threshold))
    return graph

    ### PLOTTING STUFF ###
    #max_list_size = max(len(x) for x in scores.values())
    #matrix = np.zeros((len(scores), 2 * max_list_size))
    #prefix = params['prefix']
    #for i, (query_target_pair, scores_list) in enumerate(scores.items()):
    #    matrix[i, :len(scores_list)] = np.array(sorted(scores_list))
    #plt.clf()
    #plt.matshow(matrix, cmap='cubehelix')
    #plt.title(prefix + "matrix")
    #plt.colorbar()
    #plt.savefig('figs/' + prefix +  '_matrix.pdf')

    #plt.clf()
    #stds = []
    #for i, (query_target_pair, scores_list) in enumerate(scores.items()):
    #    stds.append(np.array(scores_list).std())
    #plt.hist(stds)
    #plt.savefig('figs/' + prefix +  '_std.pdf')

    ## Why is it now plotting like I want it to?
    #plt.clf()
    #data = np.array([np.array(d) for d in list(scores.values())])
    #data = sorted(data, key=lambda x: x.mean())
    #plb.boxplot(data, whis='range', showbox=False, showmeans=True)
    #plt.savefig('figs/' + prefix +  '_box.pdf')
    ### END PLOTTING STUFF ###



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
        ## reads now may be missing the last set of numbers.  Account for this  in the node naming.
        elif n in gapset or any([g for g in gapset if n in g]):
            node_colors.append(nx_helpers.rgb_to_hex((0, 10, 250)))
        elif n in spanset or any([s for s in spanset if n in s]):
            node_colors.append(nx_helpers.rgb_to_hex((0, 250, 10)))
        else:
            # uncategorized
            node_colors.append(nx_helpers.rgb_to_hex((0, 0, 0)))
    return node_colors

def generate_graph(params):
    """generates and returns overlap graph from .paf files
    """

    alignedreads = read_paf(params, should_filter_paf=True)

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
        # Error in finding community quality
        return -1
    elif spanset_i == 0:
        return spanset_0 + gapset_1
    elif spanset_i == 1:
        return spanset_1 + gapset_0
    else:
        return -1

def make_line_plot(the_sets, params):
    """Makes an IGV-style drawing with colors"""
    bed_filename = params['bed_filename']

    coords = get_read_coordinates(bed_filename, normalize=True)
    gapset, spanset, preset, postset = the_sets
    colors = node_set_colors(coords.keys(), gapset, spanset, preset, postset)

    y_increment = (1. / float(len(coords)))
    y_values = [float(i) * y_increment for i in range(0, len(coords))]
    for i, (coord, y_value) in enumerate(zip(coords.values(), y_values)):
        plt.plot(list(coord), [y_value, y_value], color=colors[i], linestyle='-', linewidth=1.5)
    plt.axis('off')
    plt.title("IGV style line plot")

def make_four_params(args):
    m4_filename = args[0]
    the_dir = args[1]
    prefix = m4_filename.split('/')[-1].split('.')[0]
    bed_filename = the_dir + prefix + ".bed"
    fasta_filename = the_dir + prefix + ".fa"

    params = {
        'dir': the_dir,
        'min_matching_length': args[2],
        'minimap_call': './minimap',
        'prefix': prefix,
        'm4_filename': m4_filename,
        'bed_filename': bed_filename,
        'fasta_filename': fasta_filename,
        'paf_filename': temp_dir + "/tmp_" + prefix + ".paf",
        'gap_score_threshold': 0.12,
        'window_size': 500
    }
    return params

def make_four_pdf(args):
    """
    Generates four graphs for the
    structural variant defined by merged_filename
    """
    params = make_four_params(args)
    m4_filename = params['m4_filename']
    prefix = params['prefix']
    min_matching_length = params['min_matching_length']

    # if there are fewer than threshold reads then skip it
    threshold = 25 # threshold before plotting.
    if len(open(m4_filename).readlines()) < threshold:
        print('skipping %s because it has %d lines' % (
            m4_filename,
            len(open(m4_filename).readlines()))
        )
        return

    plb.rcParams['figure.figsize'] = 30, 30
    plt.clf()
    plt.figure(1)

    remove_punctuation = lambda x: ''.join(e for e in x if e.isdigit() or e == '.')
    coords = [int(remove_punctuation(a)) for a in prefix.split('_')[1:3]]
    dist = coords[1] - coords[0]

    graph = generate_graph(params)
    preset, postset, spanset, gapset = get_read_classifications(params)
    # Draw Ground Truth
    plt.subplot(2, 3, 1)
    node_colors = node_set_colors(graph.nodes(), spanset, gapset, preset, postset)
    pos = nx.spring_layout(graph)

    assert(len(node_colors) == nx.number_of_nodes(graph))
    title = "Chr {0}; L={1}; Ground Truth Colors\n\
            Red=Preset, Yellow=Postset, Blue=GapSet, Green=SpanSet\n\
            num_edges = {2}\
            "\
            .format(prefix, min_matching_length, nx.number_of_edges(graph))
    nx.draw_spring(graph, node_color=node_colors, node_size=100)
    #nx.draw(graph, node_color=node_colors, node_size=100, pos=pos)
    plt.title(title)

    # Draw histogram of smith waterman scores and remove bad edges
    plt.subplot(2, 3, 2)

    # squash preset and postset nodes
    graph = nx_helpers.remove_nodes(graph, preset)
    graph = nx_helpers.remove_nodes(graph, postset)

    # filter nodes by smith_waterman
    with utils.Timer("smith_waterman_filter"):
        graph = smith_waterman_filter(graph, params)

    # Draw groudn truth with squashed nodes
    plt.subplot(2, 3, 3)
    node_colors = node_set_colors(graph.nodes(), spanset, gapset, preset, postset)
    assert(len(node_colors) == nx.number_of_nodes(graph))
    title = "Chr {0}; L={1}; Ground Truth Colors \n\
            Removed Preset and Postsetnodes; Blue=GapSet, Green=SpanSet\n\
            number of edges = {2}"\
            .format(prefix, min_matching_length, nx.number_of_edges(graph))
    nx.draw_spring(graph, node_color=node_colors, node_size=100)
    #nx.draw(graph, node_color=node_colors, node_size=100, pos=pos)
    plt.title(title)

    # Drop Small Communities and Draw
    plt.subplot(2, 3, 4)
    communities = nx_helpers.get_communities(graph)
    graph, communities = drop_small_communities(graph, communities)
    node_colors = node_community_colors(graph, communities)
    assert(len(node_colors) == nx.number_of_nodes(graph))
    title = "Chr {0}; L={1}; After Removing Small Communities; NumCom={2}\n\
            ComQual={3}, MapQual={4}\n\
            number of edges = {5}"\
            .format(prefix, min_matching_length, len(communities),
                    community_quality(communities, spanset, gapset),
                    mapping_quality(graph, spanset, gapset),
                    nx.number_of_edges(graph))
    #nx.draw(graph, node_color=node_colors, node_size=100, pos=pos)
    nx.draw_spring(graph, node_color=node_colors, node_size=100)
    plt.title(title)

    # IGV Line Plot
    plt.subplot(2, 3, 5)
    make_line_plot((spanset, gapset, preset, postset), params)

    plt.savefig('figs/%s-communities.pdf' % (prefix))

    ret_string = '%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\tchr%s_slop5000.png\t%s-communities.pdf' % (
        prefix,
        prefix.split('_')[0],
        coords[0],coords[1],coords[1]-coords[0],
        len(communities),
        community_quality(communities, spanset, gapset),
        mapping_quality(graph, spanset, gapset),
        prefix,prefix
    )

    return ret_string

def sixteen_graphs(the_dir):
    """ generates graphs for each structual variant
    """
    # TODO change to deprecation warning
    warnings.warn("Does not call sv_pipeline functoins correctly", DeprecationWarning)

    plb.rcParams['figure.figsize'] = 30, 30
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
            preset, postset, spanset, gapset = get_read_classifications(prefix,\
                                                bed_filename, merged_filename=merged_filename)
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

def four_graphs(the_dir, min_matching_length):
    """
    Generates four graphs for each structural variant in the directory
    formerly
    """
    files = get_files(the_dir)
    print('Looking in directory %s*.m4' % (the_dir))
    print('There are %d files' % (len(files)))
    the_dirs = [the_dir for _ in files]
    min_matching_lengths = [min_matching_length for _ in files]
    zipped = zip(files, the_dirs, min_matching_lengths)

    ## print a header to screen: these values will be written at the end of the make_four_pdf()
    ## function for each input.
    header = 'prefix\tchr\tleftbp\trightbp\tdelsize\tnumcommunities\tcommunityquality\tmappingquality'

    # TODO uncomment this to make parallel
    #p = Pool()
    #results = p.map(make_four_pdf, zipped)

    # Comment this to make sychronous
    results = [make_four_pdf(z) for z in zipped]

    with open('results.txt', 'w') as results_file:
        results_file.write(header + '\n')
        results_file.write('\n'.join(results))
        results_file.write('\n')

    # TODO uncomment this to make parallel
    #p.close()

    # for testing purposes - only run one instance.
    #for z in zipped:
    #    make_four_pdf(z)

