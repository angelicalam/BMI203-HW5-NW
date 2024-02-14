# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch
import numpy as np

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    NW = NeedlemanWunsch(sub_matrix_file="./substitution_matrices/BLOSUM62.mat", 
                         gap_open=-10, gap_extend=-1)
    species = ['Gallus_gallus', 'Mus_musculus', 'Balaeniceps_rex', 'tursiops_truncatus']
    align_v_human_scores = []
    # Align all species to human
    for seq in [gg_seq, mm_seq, br_seq, tt_seq]:
        alignment_score, seqA_align, seqB_align = NW.align(hs_seq, seq)
        align_v_human_scores.append(alignment_score)
    # Sort alignment scores and species
    sorted_scores = np.sort(align_v_human_scores)[::-1]
    sorted_species = [species[i] for i in np.argsort(align_v_human_scores)[::-1]]
    # Print species in order of most similar to human BRD
    print("Species BRD2 similarity to human BRD2, in order of most to least:")
    for spec in sorted_species:
        print(spec)

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    print("\nAlignment scores:")
    for speci in range(len(sorted_species)):
        print(sorted_species[speci], sorted_scores[speci])
    

if __name__ == "__main__":
    main()
