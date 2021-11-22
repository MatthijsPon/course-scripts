#!/usr/bin/env python3
"""
Author: Matthijs Pon
Date: 2020-06-08

Description: this is a script to train a HMM and create a sequence using the
model. It is used to answer assignment 4.
"""
# Import statements
from random import random


# Background amino acid probabilities
pa = {'A': 0.074, 'C': 0.025, 'D': 0.054, 'E': 0.054, 'F': 0.047, 'G': 0.074,
      'H': 0.026, 'I': 0.068, 'L': 0.099, 'K': 0.058, 'M': 0.025, 'N': 0.045,
      'P': 0.039, 'Q': 0.034, 'R': 0.052, 'S': 0.057, 'T': 0.051, 'V': 0.073,
      'W': 0.013, 'Y': 0.034}


# Function definitions
def parse_file(filename):
    """Parse a file to a list of sequences.

    :param filename: name of file to be parsed.
    :return: a list of sequences in the file.
    """
    data_dict = {}

    with open(filename) as file:
        for line in file:
            if line[0] == ">":
                key = line[1:].strip()
                data_dict[key] = []
            else:
                data_dict[key].append(line.strip())

    for keys in data_dict.keys():
        data_dict[keys] = "".join(data_dict[keys])

    return data_dict


def is_match_state(seq_dict, location):
    """Determines if the current location is a match state.

    :param seq_dict: dictionary of strings, sequences.
    :param location: (int) current position in sequences.
    :return: boolean.
    """
    count = 0

    for keys in seq_dict:
        if seq_dict[keys][location] != "-":
            count += 1

    if count > len(seq_dict) / 2:
        return True
    else:
        return False


def calc_match_states(seq_dict, len_align):
    """Creates a list of match states in the aligned sequences.

    :param seq_dict: dictionary of strings, sequences.
    :param len_align: int, length of the alignment.
    :return: match_states: list of booleans
             n_matches: int, amount of match states.
    """
    match_states = [False for i in range(len_align)]
    n_matches = 0

    for i in range(len_align):
        if is_match_state(seq_dict, i):
            match_states[i] = True
            n_matches += 1

    return match_states, n_matches


def calc_emission(seq_dict, len_align, match_states, n_matches):
    """Calculate the emission of match states.

    :param seq_dict: dictionary of sequences.
    :param len_align: length of the alignment.
    :param match_states: list of match states.
    :param n_matches: amount of match states.
    :return: A dictionary of emission based on match positions.
    """
    mat_aa_prob = [{} for i in range(n_matches)]
    count = 0

    for i in range(len_align):

        if match_states[i]:
            mat_dict = {}

            for keys in seq_dict:
                if seq_dict[keys][i] != "-":
                    if mat_dict.get(seq_dict[keys][i]) is not None:
                        mat_dict[seq_dict[keys][i]] += 1
                    else:
                        mat_dict[seq_dict[keys][i]] = 1

            mat_aa_prob[count] = mat_dict
            count += 1

    n_sequences = len(seq_dict)

    for i in range(len(mat_aa_prob)):
        for keys in mat_aa_prob[i]:
            mat_aa_prob[i][keys] = mat_aa_prob[i][keys] / n_sequences

    return mat_aa_prob


def transition_sequence(seq, len_align, match_states):
    """Transform a sequence into a transition sequence.

    :param seq: sequence to be transformed.
    :param len_align: length of the alignment.
    :param match_states: a list of match states in the alignment.
    :return: a transition sequence as list.
    """
    trans_seq = ["M"]

    for i in range(len_align):
        if match_states[i]:
            if seq[i] != "-":
                trans_seq.append("M")
            else:
                trans_seq.append("D")
        elif seq[i] != "-":
            trans_seq.append("I")

    trans_seq.append("M")

    return trans_seq


def calc_transition(seq_dict, len_align, match_states, n_matches):
    """Calculate the transition probabilities of the sequences.

    :param seq_dict: dictionary of sequences.
    :param len_align: length of the alignment.
    :param match_states: list of match states.
    :param n_matches: amount of match states.
    :return: a dictionary of lists with transition probabilities.
    """
    trans_dict = {("M", "M"): [0 for i in range(n_matches + 1)],
                  ("M", "I"): [0 for i in range(n_matches + 1)],
                  ("M", "D"): [0 for i in range(n_matches + 1)],
                  ("I", "M"): [0 for i in range(n_matches + 1)],
                  ("I", "I"): [0 for i in range(n_matches + 1)],
                  ("I", "D"): [0 for i in range(n_matches + 1)],
                  ("D", "M"): [0 for i in range(n_matches + 1)],
                  ("D", "I"): [0 for i in range(n_matches + 1)],
                  ("D", "D"): [0 for i in range(n_matches + 1)]}

    for seq in seq_dict.values():
        trans_seq = transition_sequence(seq, len_align, match_states)
        count = 0

        for i in range(len(trans_seq) - 1):
            trans_dict[(trans_seq[i], trans_seq[i + 1])][count] += 1

            if trans_seq[i + 1] == "M" or trans_seq[i + 1] == "D":
                count += 1

    for i in range(0, n_matches + 1):
        n_m = 0
        n_i = 0
        n_d = 0

        for keys in trans_dict:
            if keys[0] == "I":
                n_i += trans_dict[keys][i]
            elif keys[0] == "M":
                n_m += trans_dict[keys][i]
            else:
                n_d += trans_dict[keys][i]

        for keys in trans_dict:
            if keys[0] == "I":
                if n_i != 0:
                    temp_list = trans_dict[keys]
                    temp_list[i] = temp_list[i] / n_i
                    trans_dict[keys] = temp_list

            elif keys[0] == "M":
                if n_m != 0:
                    temp_list = trans_dict[keys]
                    temp_list[i] = temp_list[i] / n_m
                    trans_dict[keys] = temp_list

            else:
                if n_d != 0:
                    temp_list = trans_dict[keys]
                    temp_list[i] = temp_list[i] / n_d
                    trans_dict[keys] = temp_list

    return trans_dict


def train_hmm(filename):
    """Train a HMM using aligned sequences.

    :param filename: name of file containing the sequences.
    :return: the match state emission, insertion state emission and a dict
    containing the transition probabilities.
    """
    seq_dict = parse_file(filename)
    len_align = len(seq_dict[list(seq_dict.keys())[0]])
    match_states, n_matches = calc_match_states(seq_dict, len_align)

    match_emission = calc_emission(seq_dict, len_align, match_states,
                                   n_matches)
    insertion_emission = pa
    transition_dict = calc_transition(seq_dict, len_align, match_states,
                                      n_matches)

    return match_emission, insertion_emission, transition_dict


# Not self-written, function was provided for the assignment.
def sample_emission(events):
    """Return a key from dict based on the probabilities 

    events: dict of {key: probability}, sum of probabilities should be 1.0. 
    """
    key_options = list(events.keys())
    cum = [0.0 for i in key_options]

    cum[0] = events[key_options[0]]
    for i in range(1, len(events)):
        cum[i] = cum[i - 1] + events[key_options[i]]

    # Should not be necessary, but for safety
    cum[len(cum) - 1] = 1.0

    ref_point = random()
    pick = ""
    i = 0
    while pick == "" and i < len(cum):
        if ref_point < cum[i]:
            pick = key_options[i]
        i = i + 1
    return pick


def sample_transition(trans_dict, cur_state, mat_state):
    """Return a key from the transition dictionary.

    :param trans_dict: dict containing the transition probabilities.
    :param cur_state: current state, match, insertion or deletion.
    :param mat_state: number of current match state.
    :return: returns a transition.
    """
    option_list = []

    for keys in trans_dict:

        if keys[0] == cur_state:
            option_list.append(keys)

    pick = random()

    cum = [0.0 for i in range(len(option_list))]
    cum[0] = trans_dict[option_list[0]][mat_state]

    for i in range(1, len(option_list)):
        cum[i] = cum[i - 1] + trans_dict[option_list[i]][mat_state]

    for i in range(len(option_list)):
        if pick < cum[i]:
            return option_list[i]

    return False


def create_hmm_seq(mat_em, ins_em, trans_dict):
    """Create a sequence using a trained HMM.

    :param mat_em: list of dicts of the match state emission probabilities.
    :param ins_em: dict of the insertion state emission probability.
    :param trans_dict: dict containing the transition probabilities.
    :return: a sequence generated using a trained HMM.
    """
    mat_state = 0
    cur_state = "M"
    seq = []

    while mat_state < len(mat_em):
        trans = sample_transition(trans_dict, cur_state, mat_state)
        if trans[1] == "M":
            seq.append(sample_emission(mat_em[mat_state]))
            cur_state = "M"
            mat_state += 1
        elif trans[1] == "D":
            cur_state = "D"
            mat_state += 1
        else:
            seq.append(sample_emission(ins_em))
            cur_state = "I"

    return "".join(seq)


def main():
    """Main code."""

    # Question 1
    seq_dict = parse_file("test.fasta")
    len_align = len(seq_dict[list(seq_dict.keys())[0]])
    match_states, n_matches = calc_match_states(seq_dict, len_align)
    print("Q1\nNumber of match states {0}.\n".format(n_matches))

    # Question 2
    print("test.fasta: ")
    mat_em, ins_em, trans_dict = train_hmm("test.fasta")
    print("Q2\nMatch state emission probabilities: ")
    for i in range(len(mat_em)):
        print("State {0} = {1}".format(i + 1, mat_em[i]))

    print("\nTransition state probabilities: ")
    for keys in trans_dict:
        print("{0}: {1}".format(keys, trans_dict[keys]))

    # Question 3
    print("\nQ3\nEmission probability matrix:\n\n\t", end="")

    aa = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
          "L", "K", "M", "F", "P", "S", "T", "Q", "Y", "V"]
    for char in aa:
        print("  {0}\t".format(char), end="")

    print()

    for i in range(len(mat_em)):
        for char in aa:
            if mat_em[i].get(char) is None:
                mat_em[i][char] = 0.0

        print("{0}{1}\t".format(i + 1, " A"), end="")
        for char in aa:
            print("{0:.3f}\t".format(mat_em[i][char]), end="")
        print()

    # Question 4
    print("\nQ4")
    for i in range(10):
        print("Seq {0}:".format(i + 1))
        emission = create_hmm_seq(mat_em, ins_em, trans_dict)
        print(emission)

    # Question 5
    # Q5.2
    print("test_large.fasta: ")
    mat_em, ins_em, trans_dict = train_hmm("test_large.fasta")
    print("Q2\nMatch state emission probabilities: ")
    for i in range(len(mat_em)):
        print("State {0} = {1}".format(i + 1, mat_em[i]))

    print("\nTransition state probabilities: ")
    for keys in trans_dict:
        print("{0}: {1}".format(keys, trans_dict[keys]))

    # Q5.3
    print("\nQ3\nEmission probability matrix:\n\n\t", end="")

    aa = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
          "L", "K", "M", "F", "P", "S", "T", "Q", "Y", "V"]
    for char in aa:
        print("  {0}\t".format(char), end="")

    print()

    for i in range(len(mat_em)):
        for char in aa:
            if mat_em[i].get(char) is None:
                mat_em[i][char] = 0.0

        print("{0}{1}\t".format(i + 1, " A"), end="")
        for char in aa:
            print("{0:.3f}\t".format(mat_em[i][char]), end="")
        print()

    # Q5.4
    print("\nQ4")
    for i in range(10):
        print("Seq {0}:".format(i + 1))
        emission = create_hmm_seq(mat_em, ins_em, trans_dict)
        print(emission)


if __name__ == "__main__":
    main()
