#!/usr/bin/env python3
"""
Author: Matthijs Pon
Date: 2020-11-19
Description: Calculate the average alignment distance between protein families
             using blastp and write them to a comma-separated file.
Usage: python3 tf_family_distance_matrix.py <input.fasta> <output.csv>
    input.fasta: name of the input fasta file
    output.csv: name of file to output to
"""

from sys import argv
import os
import subprocess


def check_user_input():
    """Check the input given on the command line.

    output: None, function raises errors if input is incorrect.
    """
    if len(argv) == 3:
        return None
    else:
        raise ValueError("Please give two command line arguments.")


def blastp(input_file, database):
    """Run blastp using an input fasta file and a database file.

    input:
        input_file: string, name of file to blast against db
        database: string, name of file to be used as db, will be indexed if
                  the index cannot be found

    output: None, function raises errors
    """
    # Check if DB is indexed, otherwise index it.
    if not os.path.exists(database + ".phr"):
        subprocess.check_call("makeblastdb -in {} -dbtype prot"
                              "".format(database), shell=True)

    # Run blastp
    subprocess.check_call("blastp -evalue 1E-10 -outfmt 6 -max_hsps 1 "
                          "-db {database} -out {input}_blastp.tsv "
                          "-query {input}"
                          "".format(input=input_file, database=database),
                          shell=True)
    return None


def parse_blastp(filename):
    """Parse a tab separated blastp output file.

    input:
        filename: string, filename of .tsv file

    output: dict of tuples, dictionary structured as (query_id, subject_id):
            tuple containing all other tab-separated values
    """
    parse_dict = {}

    with open(filename) as file:
        for line in file:
            values = line.strip().split("\t")
            query = values[0]
            subject = values[1]
            # Only add the line if the query and subject end in .1
            if query.split("|")[0].endswith(".1") and \
                    subject.split("|")[0].endswith(".1"):
                parse_dict[(query, subject)] = tuple(values[2:])

    return parse_dict


def tf_family_distances(parsed_blastp):
    """Make a table of the average alignment length between TF families.

    input:
        parsed_blastp: dict of tuples, the output of the parse_blastp()
        function

    output: list of lists, a table of average alignment lengths between TF
            families
    """
    families = []
    for query, subject in parsed_blastp.keys():
        query, family = query.split("|")
        # Check which families are present in the data.
        if not family in families:
            families.append(family)

    # Init TF table in same size as families x families.
    tf_table = [[[0, 0] for item in families] for item in families]

    for key in parsed_blastp.keys():
        # Get the index of the families from the families list.
        fam_query_index = families.index(key[0].split("|")[1])
        fam_subj_index = families.index(key[1].split("|")[1])
        # Sum lengths and amount of alignments
        tf_table[fam_query_index][fam_subj_index][0] += \
            int(parsed_blastp[key][1])
        tf_table[fam_query_index][fam_subj_index][1] += 1

    # Average lengths and set empty cells to None
    for row in range(len(tf_table)):
        for column, values in enumerate(tf_table[row]):
            if values[1] == 0:
                tf_table[row][column] = None
            else:
                tf_table[row][column] = values[0] / values[1]
    return tf_table, families


def write_csv(data, headers, filename):
    """Write a list of lists (each sublist containing floats) to a csv file.

    input:
        data: list of lists which contain floats, data to write to file
        headers: list of strings containing, headers will be applied to both
                 rows and columns
        filename: string, name of file to write to

    output: None, function raises errors.
    """
    with open(filename, "w+") as file:
        # Write header
        file.write(",{}\n".format(",".join(headers)))
        for row_i, row in enumerate(data):
            # Write header for each rows
            file.write("{},".format(headers[row_i]))
            for column_i, cell in enumerate(row):
                if column_i < len(row) - 1:
                    if cell:
                        file.write("{:.1f},".format(cell))
                    else:
                        file.write(",")
                # Last cell should not have a comma
                else:
                    if cell:
                        file.write("{:.1f}".format(cell))
                    else:
                        file.write("")
            file.write("\n")


def main():
    """Main function."""
    check_user_input()
    # Only run blastp if needed.
    if not os.path.exists(argv[1] + "_blastp.tsv"):
        blastp(argv[1], argv[1])

    # Parse blastp into dict
    blastp_output = parse_blastp(argv[1] + "_blastp.tsv")

    # Make table of TF-family alignments
    tf_table, families = tf_family_distances(blastp_output)

    # Write table to csv
    write_csv(tf_table, families, argv[2])


if __name__ == "__main__":
    main()
