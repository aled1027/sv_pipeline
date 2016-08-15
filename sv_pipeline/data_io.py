""" Functions should be called read_*, get_*, write_*"""
import collections
import csv
import glob
import os
import sys
from itertools import tee
import warnings

from Bio import SeqIO

from sv_pipeline import overlap
from sv_pipeline import temp_dir
from sv_pipeline import smith_waterman
from pprint import pprint

def get_paf_dict(filename):
    """Generates a dictionary mapping
    str(read0+read1) to a dictionary of information
    about their alignent.
    """
    paf_dict = {}
    with open(filename) as fin:
        for line in fin:
            row = line.strip().split()
            query_name, _, query_start, query_end, _,\
                    target_name, _, target_start, target_end, _, _, _ = row[0:12]

            query_start = int(query_start)
            query_end = int(query_end)
            target_start = int(target_start)
            target_end = int(target_end)

            paf_dict[query_name+target_name] = {
                    'query_name': query_name,
                    'query_start': query_start,
                    'query_end': query_end,
                    'target_name': target_name,
                    'target_start': target_start,
                    'target_end': target_end,
            }
    return paf_dict

def get_fasta_dict(filename):
    """Generates a dictionary mapping read name to sequence

    SeqIO.parse returns an iterator to records,
    where each record is a Bio.SeqIO.SeqRecord

    Args:
        filename (string): The path of fasta file.
    Returns
        Dictionary: Maps read names to their string
    """

    with open(filename, 'r') as fasta_file:
        ret_dict = {record.id: str(record.seq) \
                for record in SeqIO.parse(fasta_file, "fasta")}
    return ret_dict

def read_paf(params, should_filter_paf=True):
    """ Reads the PAF file from minimap.  If should_filter_paf=True,
    filters the file according to Ali's specs.
    The returned list has three values.
    The zerotih value is read0, first value is read1
    and the third value is that quality of their mapping.
    See minimap for more info on quality rating.

    TODO this function could be cleaned up
    """

    prefix = params['prefix']
    fasta_filename = params['fasta_filename']
    min_matching_length = params['min_matching_length']
    minimap_call = params['minimap_call']
    paf_filename = params['paf_filename']

    def pass_minimap_filter(row):
        """returns True if the row passes Ali's filters, False otherwise.
        From filename of previous file:
        all_hg004.mapq200.alen3000.ovtol500
        """
        map_cutoff = 200
        alen_cutoff = 3000
        ov_tol = 500

        qname, ql, qs, qe, strand, tname, tl, ts, te, _, alen, mapq = row[0:12]
        qs, qe, ql = int(qs), int(qe), int(ql)
        ts, te, tl = int(ts), int(te), int(tl)
        mapq = int(mapq)
        alen = int(alen)
        ov = overlap.overlap(qname, tname, mapq, -1, "+", qs, qe, ql, strand, ts, te, tl)
        if mapq > map_cutoff and alen > alen_cutoff and ov.hasFullOverlap(ov_tol):
            return True
        return False

    minimap_command = "{} -Sw5 -L{} -m0 {} {} > {}"\
                      .format(minimap_call,
                              min_matching_length,
                              fasta_filename,
                              fasta_filename,
                              paf_filename
                      )

    os.system(minimap_command)

    # populate alignedreads by reading in paf (filter if neccessary)
    alignedreads = []
    with open(paf_filename) as fin:
        for line in fin:
            row = line.strip().split()
            if not should_filter_paf or pass_minimap_filter(row):
                alignedreads.append([row[0], row[5], int(row[10])])
    return alignedreads

def get_read_coordinates(bed_filename, normalize=False):
    """Pulls coordinates of reads from bed file
    and returns a dictionary mapping a read to a tuple of starting coord, ending coord
    which are normalized from 0 to 1"""
    coords = None # {name: (start, end)} for each read
    with open(bed_filename, 'r') as csvfile:
        bed_reader = csv.reader(csvfile, delimiter='\t')
        coords = {row[3]: (int(row[1]), int(row[2])) for row in bed_reader}

    if normalize:
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

    paths = []
    with open(filename) as f:
        for line in f.readlines():
            # example line:
            # 92    0.00000e+00 read_43|read_45|read_46|read_71|read_73|read_64|read_96
            path = line.strip().split("\t")[2].split("|")
            path_graph = nx.Graph()
            for p0, p1 in utils.pairwise(path):
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

def get_read_classifications(params):
    """
    Returns the preset, postset, refset (aka spanset), altset (aka gapset)
    by reading data from bed_filename and the m4_filename

    TODO explain what refset, postset, gapset, etc. are

    Either m4_filename or merged_filename must be given. If both are given,
    we defer to using m4_filename
    """
    prefix = params['prefix']
    m4_filename = params['m4_filename']
    bed_filename = params['bed_filename']
    merged_filename = None # TODO remove because deprecated

    def parse_m4_file(_m4_filename):
        refset = set()
        altset = set()
        bothset = set()

        M4Data = collections.namedtuple('M4Data', ['read', 'refid', 'score', 'percent'])
        with open(_m4_filename, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=' ')
            data_dict = collections.defaultdict(list)
            for row in reader:
                read = row[0]
                data_dict[read].append(M4Data(*row[0:4]))

        for read, m4data_li in data_dict.items():
            if len(m4data_li) == 1:
                # since the read only mapped to one sequence, use that sequence
                if 'alt' in m4data_li[0]:
                    altset.add(read)
                else:
                    refset.add(read)
            elif len(m4data_li) == 2:
                # determine which is alt and which is ref
                if 'alt' in m4data_li[0].refid:
                    alt_m4data = m4data_li[0]
                    ref_m4data = m4data_li[1]
                else:
                    alt_m4data = m4data_li[1]
                    ref_m4data = m4data_li[0]

                if alt_m4data.score == ref_m4data.score:
                    bothset.add(read)
                elif alt_m4data.score < ref_m4data.score:
                    refset.add(read)
                else:
                    altset.add(read)
            else:
                # TODO remove try except, and just raise error
                try:
                    raise ValueError("There should only be two rows in the\
                            .m4 file {} for {}".format(_m4_filename, read))
                except ValueError:
                    sys
        return bothset, refset, altset

    def parse_merged_file(_merged_filename):
        refset = set()
        altset = set()
        bothset = set()

        MergedData = collections.namedtuple('MergedData', ['read', 'refid', 'set'])
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

    # compute midpoint
    remove_punctuation = lambda x: ''.join(e for e in x if e.isdigit() or e == '.')
    coords = [int(remove_punctuation(a)) for a in prefix.split('_')[1:3]]
    midpoint = sum(coords) / 2.0

    ref_coords = get_read_coordinates(bed_filename)

    # now loop over bothset to determine if a read belongs to preset or postset
    preset = set()
    postset = set()
    for read in bothset:
        ## Sometimes the read name looks like
        ## m#####_####_#####_s1_p0/##/X_Y/X_Y
        ## instead of
        ## m#####_####_#####_s1_p0/##/X_Y
        ## Try removing the last pair of numbers.
        if read not in ref_coords:
            read = '/'.join(read.split('/')[:-1])

        if ref_coords[read][0] < midpoint: # Pre node
            preset.add(read)
        else:
            postset.add(read)
    return preset, postset, refset, altset

def get_files(the_dir):
    """Grabs all files from directory the_dir that
    are of form '*merged.txt'"""

    # Modified because the files are now renamed to simply be
    # chr:start-end.m4 rather than chr:start-end-merged.m4
    #files = glob.glob(the_dir + '*merged.m4')
    files = glob.glob(the_dir + '*.m4')

    #These regions causes problems
    #Zeroith path doesn't finish girvan newman community detection
    #Not sure about first and second
    #Third and fourth raise error in community detection: invalid value encountered in double scalars
    #According to stackoverflow, Nan might be returned from some function

    #baddies = ["/data/mtsinai/2016_05_13_GR37_HG002_hapcalls/ambig_calls/14_22918113_22982906_buffer_10000_merged.txt",
    #           "/data/mtsinai/2016_05_13_GR37_HG002_hapcalls/homozygous_calls/14_105905509_105905510_buffer_10000_merged.txt",
    #           "/data/mtsinai/2016_05_13_GR37_HG002_hapcalls/homozygous_calls/14_100811401_100811402_buffer_10000_merged.txt",
    #           "/data/mtsinai/2016_05_13_GR37_HG002_hapcalls/heterozygous_hap1_calls/14_21805388_21805389_buffer_10000_merged.txt",
    #           "/data/mtsinai/2016_05_13_GR37_HG002_hapcalls/heterozygous_hap1_calls/14_104363030_104363031_buffer_10000_merged.txt",
    #           "data/ambig_calls/14_22918113_22982906_buffer_10000_merged.m4"]

    #for baddie in baddies:
    #    print(baddie)
    #    try:
    #        files.remove(baddie)
    #    except Exception as e:
    #        sys.exc_info()[0]
    #        print(e, baddie)
    return files



