"""Module that does smith-waterman"""

from __future__ import print_function
import numpy as np

import align # https://github.com/aled1027/align

from sv_pipeline import data_io

def smith_waterman(seq0, seq1):
    """
    gap_score: # of gaps / length of aligned sequence (lower number => better alignment)

    """
    score, gap_score, a0, a1 = align.align(seq0, seq1, local=False)
    return gap_score

if __name__ == '__main__':
    d = data_io.get_fasta_dict('example.fa')
    r0 = 'm150213_074729_42177R_c100777662550000001823160908051505_s1_p0/70715/9957_22166'
    r1 = 'm150126_093705_42156_c100779662550000001823165208251525_s1_p0/144605/28461_40297'
    sq0 = d[r0][:1000]
    sq1 = d[r1][:1000]

    score = smith_waterman(sq0, sq1)
    print(score)



