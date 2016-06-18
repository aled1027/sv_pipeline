""" Functions should be called read_*, get_*, write_*"""
import collections
import csv
import glob
import os
from itertools import tee
import warnings

from sv_pipeline import overlap
from sv_pipeline import temp_dir

def read_paf(prefix, fasta_filename, min_matching_length, should_filter_paf=True):
    """ Reads the PAF file from minimap.  If should_filter_paf=True,
    filters the file according to Ali's specs.
    The returned list has three values.
    The zerotih value is read0, first value is read1
    and the third value is that quality of their mapping.
    See minimap for more info on quality rating.
    """

    def pass_filters(row):
        """returns True if the row passes Ali's filters, False otherwise.
        From filename of previous file:
        all_hg004.mapq200.alen3000.ovtol500
        """
        map_cutoff = 200
        alen_cutoff = 3000
        ov_tol = 500

        kept = False
        qname, ql, qs, qe, strand, tname, tl, ts, te, _, alen, mapq = row[0:12]
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

def get_line_plot_coords(bed_filename):
    """Pulls coordinates of reads from bed file
    and returns a dictionary mapping a read to a tuple of starting coord, ending coord
    which are normalized from 0 to 1"""
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
    return coords

def get_ref_coords(bed_filename):
    """Given a bed_filename, looks up the refcoods.bed file
    and parses it.

    TODO this and get_line_plot_coords might be the same function?
    """
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

def read_paths(filename):
    """Reads data that pathlinker wrote, which is a list of paths.
    read_paths reads in the paths, transforms each path into a networkx graph,
    and returns a list of the graphs.
    """
    def pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = tee(iterable)
        next(b, None)
        return zip(a, b)

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

def get_read_classifications(prefix, bed_filename, m4_filename=None, merged_filename=None):
    """
    Returns the preset, postset, refset (aka spanset), altset (aka gapset)
    by reading data from bed_filename and the m4_filename

    TODO explain what refset, postset, gapset, etc. are

    Either m4_filename or merged_filename must be given. If both are given,
    we defer to using m4_filename
    """

    def parse_m4_file(_m4_filename):
        warnings.warn("Untested function")
        refset = set()
        altset = set()
        bothset = set()
        with open(merged_filename, 'r') as fin:
            for line in fin:
                row = line.strip().split()
                read = row[0]
                if row[2] == 'b': # BOTH reference and alternate genomes
                    bothset.add(read)
                elif row[2] == 'a': # ALTERNATE genome only
                    altset.add(read)
                elif row[2] == 'r': # REFERENCE genome only
                    refset.add(read)
                else:
                    raise ValueError("Unexpected value in column 2")
        return bothset, refset, altset

    def parse_merged_file(_merged_filename):
        refset = set()
        altset = set()
        bothset = set()

        MergedData = collections.namedtuple('MergedData', ['read', 'refid', 'set'])
        temp_data_dict = collections.defaultdict(list)
        with open(_merged_filename, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=' ')
            datas = [MergedData(*row) for row in reader]
            for data in datas:
                if data.set == 'a':
                    altset.add(data.read)
                elif data.set == 'r':
                    refset.add(data.read)
                elif data.set == 'b':
                    bothset.add(data.read)
                else:
                    raise ValueError("Unexpected value for set")
        return bothset, refset, altset

    if m4_filename:
        bothset, refset, altset = parse_m4_file(m4_filename)
    elif merged_filename:
        bothset, refset, altset = parse_merged_file(merged_filename)
    else:
        raise ValueError("Either m4_filename or merged_filename must not be None")

    # get coordinates by parsing bedfile
    remove_punctuation = lambda x: ''.join(e for e in x if e.isdigit() or e == '.')
    coords = [int(remove_punctuation(a)) for a in prefix.split('_')[1:3]]
    midpoint = sum(coords) / 2.0
    bed_lines = get_ref_coords(bed_filename) # get the ref coordinates from BED file
    _, left_coords = get_overlaps(bed_lines)

    # now loop over bothset to determine if a read belongs to preset or postset
    preset = set()
    postset = set()
    for read in bothset:
        if left_coords[read] < midpoint: # Pre node
            preset.add(read)
        else:
            postset.add(read)
    return preset, postset, refset, altset

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


