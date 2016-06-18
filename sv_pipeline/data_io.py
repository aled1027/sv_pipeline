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

def old_get_read_classifications(prefix, bed_filename, merged_filename):
    """Reads read classifications from disk and returns them.
    Used for "ground truth" of graph.
    """
    warnings.warn("Is deprecated")
    refset = set()
    altset = set()
    preset = set()
    postset = set()

    ## get coordinates
    remove_punctuation = lambda x: ''.join(e for e in x if e.isdigit() or e == '.')
    coords = [int(remove_punctuation(a)) for a in prefix.split('_')[1:3]]
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

        returns refset, altset, preset, postset
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
                raise ValueError("There should only be two rows in the .m4 file {} for {}".format(m4_file, queryid))

def get_read_classifications(prefix, bed_filename, m4_filename):
    """
    Returns the preset, postset, refset (aka spanset), altset (aka gapset)
    by reading data from bed_filename and the m4_filename

    TODO explain what refset, postset, gapset, etc. are
    """
    warnings.warn("Untested function")
    ## get coordinates
    remove_punctuation = lambda x: ''.join(e for e in x if e.isdigit() or e == '.')
    coords = [int(remove_punctuation(a)) for a in prefix.split('_')[1:3]]
    midpoint = sum(coords) / 2.0
    bed_lines = get_ref_coords(bed_filename) # get the ref coordinates from BED file
    _, left_coords = get_overlaps(bed_lines)

    M4Data = collections.namedtuple('M4Data', ['queryid', 'refid', 'score', 'percent_good'])

    refset = set()
    altset = set()
    preset = set()
    postset = set()
    temp_data_dict = collections.defaultdict(list)

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
                if 'alt' in m4data_li[0].refid:
                    altset.add(m4data_li[0].queryid)
                else:
                    refset.add(m4data_li[0].queryid)
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
                    # 'b'
                elif alt_m4data.score < ref_m4data.score:
                    # it's in the refset
                    refset.add(ref_m4data.queryid)
                else:
                    # it's in the altset
                    altset.add(alt_m4data.queryid)
            else:
                raise ValueError("There should only be two rows in the .m4 file {} for {}".format(m4_file, queryid))




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


