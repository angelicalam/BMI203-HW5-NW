# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    NW = NeedlemanWunsch(sub_matrix_file="./substitution_matrices/BLOSUM62.mat", 
                         gap_open=-10, gap_extend=-1)
    alignment_score, seqA_align, seqB_align = NW.align(seq1, seq2)
    # MYQR
    # M-QR
    # 5 + (-10-1) + 5 + 5
    assert alignment_score == 4, "Alignment score is not optimal"
    assert seqA_align == "MYQR", "aligned seq1 is not optimal"
    assert seqB_align == "M-QR", "aligned seq2 is not optimal"
    # Check that the alignment matrix, gap matrices, and backtrace matrix are correct
    assert np.array_equal( NW._align_matrix, np.array([ [  0, -11, -12, -13],
                                                        [-11,   5,  -6,  -7],
                                                        [-12,  -6,   4,  -7],
                                                        [-13,  -7,  -1,   5],
                                                        [-14,  -8,  -6,   4]  ]) )
    assert np.array_equal( NW._gapA_matrix, np.array([  [0, 0, 0, 0],
                                                        [1, 0, 1, 1],
                                                        [1, 0, 0, 1],
                                                        [1, 0, 0, 0],
                                                        [1, 0, 0, 0]  ]) )
    assert np.array_equal( NW._gapB_matrix, np.array([  [0, 1, 1, 1],
                                                        [0, 0, 0, 0],
                                                        [0, 1, 0, 0],
                                                        [0, 1, 0, 0],
                                                        [0, 1, 0, 0]  ]) )
    assert np.array_equal( NW._back, np.array([  [-1, 2, 2, 2],
                                                 [ 1, 0, 2, 2],
                                                 [ 1, 1, 0, 2],
                                                 [ 1, 1, 0, 0],
                                                 [ 1, 1, 0, 0]  ]) )

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    NW = NeedlemanWunsch(sub_matrix_file="./substitution_matrices/BLOSUM62.mat", 
                         gap_open=-10, gap_extend=-1)
    alignment_score, seqA_align, seqB_align = NW.align(seq4, seq3)
    # M---QLIRHP
    # MAVHQLIRRP
    # 5 + (-10-1) + (-1) + (-1) + 5 + 4 + 4 + 5 + 0 + 7
    assert alignment_score == 17, "Alignment score is not optimal"
    assert seqA_align == "M---QLIRHP", "aligned seq4 is not optimal"
    assert seqB_align == "MAVHQLIRRP", "aligned seq3 is not optimal"



