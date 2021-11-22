#!/usr/bin/env python3
"""
Author: Matthijs Pon
Date: 25-05-2020

Description: a script for the alignment of two protein sequences.
"""
# import statements here
import math

# functions between here and __main__
blosum = """
# http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
   I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
   L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
   K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
   M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
   F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
   P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
   S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
   T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
   W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
   Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
   V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
   B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
   Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
   * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
"""


# Not self-written, function was provided for the assignment.
def blosum62():
    """Return order and similarity scores from BLOSUM62 matrix

    order: dict of {res: idx_in_matrix}
    blosum_matrix: list of lists with similarity scores
    """
    order = {}
    blosum_matrix = []
    for line in blosum.split('\n'):
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) == 24:
            for idx, sym in enumerate(parts):
                order[sym] = idx
        else:
            # list around the map construction for python3 compatibility
            blosum_matrix.append(list(map(int, parts[1:])))
    return order, blosum_matrix


BLOSUM62_ORDER, BLOSUM62_MATRIX = blosum62()


# Not self-written, function was provided for the assignment.
def score(res1, res2):
    """Return similarity score from BLOSUM62 matrix for two residues
    
    res1: string, amino acid
    res2: string, amino acid
    """
    lookup1 = BLOSUM62_ORDER[res1]
    lookup2 = BLOSUM62_ORDER[res2]
    return BLOSUM62_MATRIX[lookup1][lookup2]


# write your own functions below here
def create_matrix(seq1, seq2):
    """Create a matrix which is the size: len(seq2)*len(seq1).

    :param seq1: sequence1 which is being aligned.
    :param seq2: sequence2 which is being aligned.
    :return: created matrix (list of lists).
    """
    matrix = [[0 for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1)]
    return matrix


def score_cell(matrix, traceback, seq1_i, seq2_j, seq1, seq2, gap_pen,
               end_gap_pen):
    """Function calculates and stores the score for a particular cell using
    the score function. Next to that, it also stores the direction from which
    the score was calculated in a matrix.

    :param matrix: list of lists containing all scores.
    :param traceback: list of lists containing steps taken to calc scores.
    :param seq1_i: the current index of sequence1.
    :param seq2_j: the current index of sequence2.
    :param seq1: sequence1.
    :param seq2: sequence2.
    :param gap_pen: penalty for creating a gap.
    :param end_gap_pen: penalty for creating an end-gap.
    :return: no return. The given matrices are altered.
    """
    if seq1_i == 0 and seq2_j == 0:
        matrix[seq1_i][seq2_j] = 0
        traceback[seq1_i][seq2_j] = (0, 0, 0)
    elif seq1_i == 0:
        matrix[seq1_i][seq2_j] = matrix[seq1_i][seq2_j - 1] - end_gap_pen
        traceback[seq1_i][seq2_j] = (0, 1, 0)
    elif seq2_j == 0:
        matrix[seq1_i][seq2_j] = matrix[seq1_i - 1][seq2_j] - end_gap_pen
        traceback[seq1_i][seq2_j] = (0, 0, 1)
    else:
        # score for diagonal parent
        diagonal = matrix[seq1_i - 1][seq2_j - 1] + \
                   score(seq1[seq1_i - 1], seq2[seq2_j - 1])
        # side parent
        if seq1_i < len(matrix) - 1:
            side = matrix[seq1_i][seq2_j - 1] - gap_pen
        else:
            side = matrix[seq1_i][seq2_j - 1] - end_gap_pen
        # top parent
        if seq2_j < len(matrix[0]) - 1:
            top = matrix[seq1_i - 1][seq2_j] - gap_pen
        else:
            top = matrix[seq1_i - 1][seq2_j] - end_gap_pen
        max_score = max(diagonal, side, top)

        # set score in matrix and put parent location in traceback.
        # diagonal = (1, 0, 0); side = (0, 1, 0); top = (0, 0, 1)
        if max_score == diagonal:
            matrix[seq1_i][seq2_j] = max_score
            traceback[seq1_i][seq2_j] = (1, 0, 0)
        elif max_score == side:
            matrix[seq1_i][seq2_j] = max_score
            traceback[seq1_i][seq2_j] = (0, 1, 0)
        else:
            matrix[seq1_i][seq2_j] = max_score
            traceback[seq1_i][seq2_j] = (0, 0, 1)


def score_matrix(matrix, traceback, seq1, seq2, gap_pen, end_gap_pen):
    """Fills both the score and traceback matrices with scores and movements
    using the score_cell function.

    :param matrix: score matrix.
    :param traceback: traceback matrix.
    :param seq1: sequence1.
    :param seq2: sequence2.
    :param gap_pen: penalty for creating a gap.
    :param end_gap_pen: penalty for creating an end-gap.
    :return: no return. The given matrices are altered.
    """
    for j in range(len(seq2) + 1):
        for i in range(len(seq1) + 1):
            score_cell(matrix, traceback, i, j, seq1, seq2, gap_pen,
                       end_gap_pen)


def traceback_alignment(move_i, move_j, traceback_matrix):
    """Creates a path from the traceback matrix.

    :param move_i: starting point of seq1.
    :param move_j: starting point of seq2.
    :param traceback_matrix: the traceback matrix.
    :return: traceback path as a list.
    """
    traceback_path = []
    movement = traceback_matrix[move_i][move_j]

    # check where
    if move_i < len(traceback_matrix) - 1:
        for i in range(len(traceback_matrix) - move_i - 1):
            traceback_path.append((move_i + i + 1, move_j))

    elif move_j < len(traceback_matrix[0]) - 1:
        for i in range(len(traceback_matrix[0]) - move_j - 1):
            traceback_path.append((move_i, move_j + i + 1))

    while movement != (0, 0, 0):
        traceback_path.insert(0, (move_i, move_j))
        movement = traceback_matrix[move_i][move_j]

        if movement[0] != 0:
            # diagonal
            move_i -= 1
            move_j -= 1

        elif movement[1] != 0:
            # side
            move_j -= 1

        else:
            # up
            move_i -= 1

    return traceback_path


def string_alignment(traceback_path, seq1, seq2):
    """Uses the traceback path to create the aligned strings.

    :param traceback_path: traceback_path created by traceback_alignment().
    :param seq1: sequence1.
    :param seq2: sequence2.
    :return: returns a tuple of (align_s1, align_s2, align_info).
    """
    string1 = ""
    string2 = ""
    alignment_information = ""

    # Removes the (0, 0) from the traceback path.
    if traceback_path[0] == (0, 0):
        traceback_path.pop(0)

    for i in range(len(traceback_path)):
        s1_int = traceback_path[i][0] - 1
        if s1_int == -1 or s1_int == traceback_path[i - 1][0] - 1:
            string1 += "-"
            char1 = None
        else:
            char1 = seq1[s1_int]
            string1 += char1

        s2_int = traceback_path[i][1] - 1
        if s2_int == -1 or s2_int == traceback_path[i - 1][1] - 1:
            string2 += "-"
            char2 = None
        else:
            char2 = seq2[s2_int]
            string2 += char2

        if char1 == char2 and char1 is not None:
            alignment_information += "|"
        else:
            alignment_information += " "

    return string1, string2, alignment_information


def calc_perc_identity(aligned_seq1, aligned_seq2):
    """Calculates the percentage identity.

    :param aligned_seq1: aligned sequence1.
    :param aligned_seq2: aligned sequence2.
    :return: returns the percentage identity.
    """
    total_length = len(aligned_seq1)
    equal_identities = 0
    for i in range(total_length):
        if aligned_seq1[i] == aligned_seq2[i] and aligned_seq1[i] != "-":
            equal_identities += 1

    return equal_identities / total_length * 100


def max_score_matrix(matrix):
    """Look up the max score in the outer borders of an alignment matrix.

    :param matrix: the alignment matrix
    :return: the coordinates of the max score, and the score itself.
    """
    int_j = None
    int_i = None
    n_rows = len(matrix)
    n_columns = len(matrix[0])
    # init a max score.
    max_score = matrix[0][n_columns - 1]

    for i in range(n_rows):

        if matrix[i][n_columns - 1] > max_score:
            max_score = matrix[i][n_columns - 1]
            int_i = i
            int_j = n_columns - 1

        if i == n_rows:
            for cell in matrix[i]:

                if cell > max_score:
                    int_i = i
                    int_j = matrix[i].index(cell)

    return int_i, int_j, max_score


def align_sequences(seq1, seq2, gap_pen=0, end_gap_pen=0):
    """Aligns two sequences using a linear gap_penalty.

    :param seq1: sequence1.
    :param seq2: sequence2.
    :param gap_pen: penalty for creating a gap.
    :param end_gap_pen: penalty for creating an end-gap.
    :return: a tuple of aligned strings, the percentage identity and the
    alignment score.
    """
    alignment_matrix = create_matrix(seq1, seq2)
    traceback_matrix = create_matrix(seq1, seq2)
    score_matrix(alignment_matrix, traceback_matrix, seq1, seq2, gap_pen,
                 end_gap_pen)

    int_i, int_j, max_score = max_score_matrix(alignment_matrix)

    traceback_path = traceback_alignment(int_i, int_j, traceback_matrix)
    aligned_seq = string_alignment(traceback_path, seq1, seq2)

    perc_id = calc_perc_identity(aligned_seq[0], aligned_seq[1])

    return aligned_seq, perc_id, max_score


def print_seqs(sequence_tuple):
    """Prints the sequences in an neat way.

    :param sequence_tuple: sequence tuple from align_sequences().
    """
    for i in range(math.ceil(len(sequence_tuple[0]) / 100)):
        print("{0}\t{1}".format(i, sequence_tuple[0][i * 100:(i + 1) * 100]))
        print("\t{0}".format(sequence_tuple[2][i * 100:(i + 1) * 100]))
        print("\t{0}".format(sequence_tuple[1][i * 100:(i + 1) * 100]))
    print()


def main():
    seq1 = "THISLINE"
    seq2 = "ISALIGNED"

    # seq3: GPA1_ARATH
    seq3 = ("MGLLCSRSRHHTEDTDENTQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIFKQIKLLFQ"
            "TGFDEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDSAKYMLSSESIAIGEKLSEIGGRL"
            "DYPRLTKDIAEGIETLWKDPAIQETCARGNELQVPDCTKYLMENLKRLSDINYIPTKEDVLYA"
            "RVRTTGVVEIQFSPVGENKKSGEVYRLFDVGGQRNERRKWIHLFEGVTAVIFCAAISEYDQTL"
            "FEDEQKNRMMETKELFDWVLKQPCFEKTSFMLFLNKFDIFEKKVLDVPLNVCEWFRDYQPVSS"
            "GKQEIEHAYEFVKKKFEELYYQNTAPDRVDRVFKIYRTTALDQKLVKKTFKLVDETLRRRNLL"
            "EA")
    # seq4: GPA1_BRANA
    seq4 = ("MGLLCSRSRHHTEDTDENAQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIFKQASS"
            "DKRKIIKLLFQTGFDEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDPAKYTLSSEN"
            "MAIGEKLSEIGARLDYPRLTKDLAEGIETLWNDPAIQETCSRGNELQVPDCTKYLMENLK"
            "RLSDVNYIPTKEDVLYARVRTTGVVEIQFSPVGENKKSGEVYRLFDVGGQRNERRKWIHL"
            "FEGVTAVIFCAAISEYDQTLFEDEQKNRMMETKELFDWVLKQPCFEKTSIMLFLNKFDIF"
            "EKKVLDVPLNVCEWFRDYQPVSSGKQEIEHAYEFVKKKFEELYYQNTAPDRVDRVFKIYR"
            "TTALDQKLVKKTFKLVDETLRRRNLLEAGLL")

    # Question 1
    aligned_seqs, perc_iden, align_score = align_sequences(seq1, seq2,
                                                           gap_pen=4)
    print("Question 1\nAlignment of {0} and {1}, gap penalty = {2}.\n"
          "Percentage identity: {3:.2f}%\t Alignment score: {4}."
          "\n\nAlignment:".format("seq1", "seq2", 4, perc_iden, align_score))
    print_seqs(aligned_seqs)

    aligned_seqs, perc_iden, align_score = align_sequences(seq1, seq2,
                                                           gap_pen=8)
    print("Question 1\nAlignment of {0} and {1}, gap penalty = {2}.\n"
          "Percentage identity: {3:.2f}%\t Alignment score: {4}."
          "\n\nAlignment:".format("seq1", "seq2", 8, perc_iden, align_score))
    print_seqs(aligned_seqs)

    # Question 3
    print("Question 3\n")
    for i in range(20):
        j = i + 1
        aligned_seqs, perc_iden, align_score = align_sequences(seq1, seq2,
                                                               gap_pen=j)
        print("Question 1\nAlignment of {0} and {1}, gap penalty = {2}.\n"
              "Percentage identity: {3:.2f}%\t Alignment score: {4}."
              "\n\nAlignment:".format("seq1", "seq2", j, perc_iden,
                                      align_score))
        print_seqs(aligned_seqs)

    # Question 4 & 5
    print("Question 4, 5\n")

    aligned_seqs, perc_iden, align_score = \
        align_sequences(seq3, seq4, gap_pen=5, end_gap_pen=1)
    print("Alignment of {0} and {1}, gap penalty = {2}, end gap penalty = {3}"
          ".\nPercentage identity: {4:.2f}%\t Alignment score: {5}.\n\n"
          "Alignment:".format("seq3", "seq4", 5, 1, perc_iden, align_score))
    print_seqs(aligned_seqs)

    aligned_seqs, perc_iden, align_score = \
        align_sequences(seq3, seq4, gap_pen=5, end_gap_pen=10)
    print("Alignment of {0} and {1}, gap penalty = {2}, end gap penalty = {3}"
          ".\nPercentage identity: {4:.2f}%\t Alignment score: {5}.\n\n"
          "Alignment:".format("seq3", "seq4", 5, 10, perc_iden, align_score))
    print_seqs(aligned_seqs)


if __name__ == "__main__":
    main()
