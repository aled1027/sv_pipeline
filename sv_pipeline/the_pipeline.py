"""This script generates and analyzes overlap graphs of fasta files.
"""

# pylint: disable=invalid-name

from __future__ import with_statement, print_function, generators
import csv
import glob
from multiprocessing import Pool
import pprint
import os
import sys
import matplotlib.pyplot as plt
from pylab import rcParams
import numpy as np
import networkx as nx
import collections

from sv_pipeline import networkx_helpers as nx_helpers
from sv_pipeline import overlap
from sv_pipeline import temp_dir


## Download the compressed file from this site and extract it:
#https://xritza01.u.hpc.mssm.edu/trios/2016-05-12-data-for-networks/GR38/NA19240/dels/

def classify_reads(the_dir, save_to_disk=True):
    """For each "*_merged.m4" in directory the_dir, create a new file
    *_merged.txt. The file is a space separated file with three columns.
    Column 1 is the name of the base genome.
    Column 2 is the name of the read.
    Column 3 is an indicator of how the read aligns.
    We use the following code:
        b: the read aligns to both sequences equally well.
        r: the read aligns to the reference sequence substantially better
        a: the read aligns to the alternate sequence substantially better
        TODO compare their margins
        TODO we also, technically, don't need to write this data to disk.
        We can keep it in memory and go from their: combine this function with get_read classifications.
    """


    print("in classify_reads")
    m4_files = glob.glob(the_dir + '*merged.m4')


    M4Data = collections.namedtuple('M4Data', ['queryid', 'refid', 'score', 'percent_good'])

    for m4_file in m4_files:
        print("looking at {}".format(m4_file))

        """classifications maps a read to a classification; the primary
        role of this function is to correctly populate the classifications dictionary

        The temp_data_dict holds the information from the file, which we parse
        as we loop over it below.
        """

        temp_data_dict = collections.defaultdict(list)
        classifications = {} # maps a read to a classification

        # populate data_dict with the information about each read
        with open(m4_file, 'r') as csvfile:
            m4_reader = csv.reader(csvfile, delimiter=' ')
            for row in m4_reader:
                queryid = row[0]
                refid = row[1]
                score = row[2]
                percent_good = row[3]

                m4data = M4Data(*row[0:4])
                temp_data_dict[queryid].append(m4data)

        for queryid, m4data_li in temp_data_dict.items():
            if len(m4data_li) == 1:
                # since the read only mapped to one sequence, use that sequence
                the_class = 'a' if 'alt' in m4data_li[0].refid else 'r'
                classifications[queryid] = (queryid, m4data_li[0].refid, the_class)

            elif len(m4data_li) == 2:
                # determine which is alt and which is ref
                if 'alt' in m4data_li[0].refid:
                    alt_m4data = m4data_li[0]
                    ref_m4data = m4data_li[1]
                else:
                    alt_m4data = m4data_li[1]
                    ref_m4data = m4data_li[0]

                # compare alt and ref scores and classify accordingly
                if alt_m4data.score == ref_m4data.score:
                    classifications[queryid] = (queryid, m4data_li[0].refid, 'b')
                elif alt_m4data.score < ref_m4data.score:
                    classifications[queryid] = (queryid, m4data_li[0].refid, 'r')
                else:
                    classifications[queryid] = (queryid, m4data_li[0].refid, 'a')
            else:
                # TODO uncomment
                #raise ValueError("There should only be two rows in the .m4 file {} for {}".format(m4_file, queryid))
                pass

        if save_to_disk:
            # Write to the new _merged.txt
            new_filename = m4_file[:-3] + ".txt"
            with open(new_filename, 'w') as csvfile:
                writer = csv.writer(csvfile, delimiter=' ')
                writer.writerows(classifications.values())

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

def read_paf(prefix, fasta_filename, min_matching_length, should_filter_paf=True):
    """ Reads the PAF file from minimap.  If should_filter_paf=True,
    filters the file according to Ali's specs.
    The returned list has three values.
    The zerotih value is read0, first value is read1
    and the third value is that quality of their mapping.
    See minimap for more info on quality rating.
    """

    import tempfile

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
    remove_punctuation = lambda x: ''.join(e for e in x if e.isdigit() or e == '.')
    coords = [int(remove_punctuation(a)) for a in prefix.split('_')[1:3]]
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
                #print('ERROR:', row)
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

def get_files(the_dir):
    """Grabs all files from directory the_dir that
    are of form '*merged.txt'"""

    files = glob.glob(the_dir + '*merged.txt')

    #These regions causes problems
    #Zeroith path doesn't finish girvan newman community detection
    #Not sure about first and second
    #Third and fourth raise error in community detection: invalid value encountered in double scalars
    #According to stackoverflow, Nan might be returned from some function
    baddies = ["/data/mtsinai/2016_05_13_GR37_HG002_hapcalls/ambig_calls/14_22918113_22982906_buffer_10000_merged.txt",
               "/data/mtsinai/2016_05_13_GR37_HG002_hapcalls/homozygous_calls/14_105905509_105905510_buffer_10000_merged.txt",
               "/data/mtsinai/2016_05_13_GR37_HG002_hapcalls/homozygous_calls/14_100811401_100811402_buffer_10000_merged.txt",
               "/data/mtsinai/2016_05_13_GR37_HG002_hapcalls/heterozygous_hap1_calls/14_21805388_21805389_buffer_10000_merged.txt",
               "/data/mtsinai/2016_05_13_GR37_HG002_hapcalls/heterozygous_hap1_calls/14_104363030_104363031_buffer_10000_merged.txt"]

    for baddie in baddies:
        print(baddie)
        try:
            files.remove(baddie)
        except Exception as e:
            print(e, baddie)
    return files


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
        plt.savefig(prefix + '-16-communities.pdf')
        plt.clf()


