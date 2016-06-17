from __future__ import print_function
#import numpy as np
#import networkx as nx

class overlap:

    def __init__ (self, name1, name2, score, pctiden, strand1, start1, end1, len1, strand2, start2, end2, len2):
        self.name1, self.name2 = name1, name2
        self.score, self.pctiden = score, pctiden
        self.strand1, self.start1, self.end1, self.len1 = strand1, start1, end1, len1
        self.strand2, self.start2, self.end2, self.len2 = strand2, start2, end2, len2
        # added some code to help with alternative alignment output formats
        if self.strand1 == "+":
            self.strand1 = 0
        elif self.strand1 == "-":
            self.strand1 = 1

        if self.strand2 == "+":
            self.strand2 = 0
        elif self.strand2 == "-":
            self.strand2 = 1

    def __str__ (self):
        outstr = " ".join(map(str, [ self.name1, self.name2, self.score, self.pctiden, \
                                         self.strand1, self.start1, self.end1, self.len1, \
                                         self.strand2, self.start2, self.end2, self.len2]))
        return outstr

    def contained ( self, tolerance=200 ):
        if self.contained1() or self.contained2():
            return True
        else:
            return False

    def contained1 ( self, tolerance=200 ):
        if self.start1 < tolerance and (self.len1 - self.end1 < tolerance):
            return True
        else:
            return False

    def contained2 ( self, tolerance=200 ):
        if self.start2 < tolerance and (self.len2 - self.end2 < tolerance):
            return True
        else:
            return False

    def aligninfo (self, name):
        if name == self.name1:
            return self.strand1, self.start1, self.end1, self.len1
        elif name == self.name2:
            return self.strand2, self.start2, self.end2, self.len2
        else:
            raise Exception(name)

    def __str__ ( self ):
        print("%s %s %s %s" %(self.name1, self.name2, self.score, self.pctiden))

    def hasFullOverlap(self, tolerance=200):
        # query contained

        #print "checking contained:"
        #print self.start1, self.end1, self.len1, self.strand1
        #print self.start2, self.end2, self.len2, self.strand2
        if self.contained(tolerance): return True


        if self.strand1 == 0:
            if self.strand2 == 0:
                # --->
                #   --->
                if (self.len1 - self.end1 < tolerance) and (self.start2 < tolerance): return True
                #   --->
                # --->
                if (self.len2 - self.end2 < tolerance) and (self.start1 < tolerance): return True

            else:
                # --->
                #   <---
                if (self.len1 - self.end1 < tolerance) and (self.len2-self.end2 < tolerance): return True
                #   --->
                # <---
                if (self.start2 < tolerance) and (self.start1 < tolerance): return True

        else:
            if self.strand2 == 0:
                # <---
                #   --->
                if (self.start1 < tolerance) and (self.start2 < tolerance): return True
                #   <---
                # --->
                if (self.len2 - self.end2 < tolerance) and (self.len1 - self.end1 < tolerance): return True

        return False


    def __getitem__ (self, name):
        return aligninfo(name)

class read_overlaps:

    def __init__ (self, readname, readlen):
        self.read_name = readname
        self.read_len = readlen
        self.overlap_dict = {}

    def __len__ ( self ):
        return len(self.overlap_dict.keys())

    def add_overlap (self, o_id, overlap):
        if o_id in self.overlap_dict:
            if overlap.score > self.overlap_dict[o_id].score:
                self.overlap_dict[o_id] = overlap
        else:
            self.overlap_dict[o_id] = overlap

    def leftoverlaps ( self, thresh=200, containedCheck=False):
        left_overlaps = []
        for overlap in self.overlaps:
            if overlap.start2 < thresh and not (containedCheck and overlap.contained()):
                left_overlaps.append(overlap)
        return left_overlaps

    def rightoverlaps ( self, thresh=200, containedCheck=False):
        right_overlaps = []
        for overlap in self.overlaps:
            if overlap.len2-overlap.end2 < thresh and not (containedCheck and overlap.contained()):
                right_overlaps.append(overlap)
        return right_overlaps

    def __getitem__ (self, tname):
        return self.overlap_dict[tname]
    @property
    def overlaps ( self ):
        return self.overlap_dict.values()

#    def allBreakpointsOverlapped (self, edgeThresh=200, overlapThresh=5):
#        ''' Check to see if the read is a chimera by confirming all bp positions are covered'''
#        bpcov = np.zeros (self.read_len)
#        for overlap in self.overlaps:
#            start = overlap.start2 - overlapThresh
#            end = overlap.end2 - overlapThresh
#            bpcov[start:end] += 1
#        bpcov_internal = bpcov[edgeThresh:-edgeThresh]
#        if np.any(bpcov_internal == 0):
  #          return False
 #       else:
   #         return True
