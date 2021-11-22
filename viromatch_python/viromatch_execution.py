#!/usr/bin/python3
"""
Author: Matthijs Pon
Date: 2021-08-06

Description: Script to run the ViroMatch pipeline for multiple samples in series.
The ViroMatch pipeline is described in the following paper: https://doi.org/10.1128/MRA.01468-20
"""


import ast
import shlex
import sys
import subprocess
import os.path
import re


class ViroMatchWrapper:
    def __init__(self):
        self.settings = {}
        self.samples = []
        self.current_sample = None
        self.reg_exp_present = False
        self.parse_settings()

    def parse_settings(self):
        """Parse the settings.txt file

        :return: All the settings are put into the self.settings dictionary, except for the samples variable, which is
        put in self.samples.
        """
        with open("settings.txt") as file:
            current_var = ""
            for line in file:
                split_line = shlex.split(line)
                # Skip comments or empty lines
                if line == "\n" or split_line[0] == "#":
                    pass
                # Work with VAR lines
                elif split_line[0] == "VAR":
                    current_var = split_line[1]
                    # Check if regular expression is present
                    if current_var == "regular_expression" and split_line[3] != "":
                        self.reg_exp_present = True
                        self.settings[current_var] = split_line[3]
                    elif current_var == "sample_names":
                        if self.reg_exp_present:
                            pass
                        else:
                            self.samples = ast.literal_eval(split_line[3])
                    # Put the var in the settings dict
                    else:
                        self.settings[current_var] = split_line[3]
                # If the line is none of the above, it must be a continuation of a previous line
                else:
                    self.settings[current_var] += line

    def check_settings(self):
        """A preliminary check to see if the settings have been correctly filled in

        :return: Errors when something is wrong.
        """
        # Check if samples or regular expression are given
        if not self.reg_exp_present:
            if not self.samples:
                # No samples or reg expression means not enough info
                sys.exit("No samples or regular expression given. Please give either one in settings.txt.")
            else:
                # Samples are given, so convert them to a list of tuples
                for i in range(len(self.samples)):
                    if len(self.samples[i]) != 2:
                        sys.exit("Something wrong with the sample format. Please give a list of tuples (max 2 items"
                                 " per tuple. Error in {}".format(self.samples[i]))
        # Check if the cores setting is an actual integer
        try:
            int(self.settings["cores"])
        except ValueError:
            sys.exit("Something is wrong with the cores parameter. It is not an integer.")

    def search_samples(self):
        """Search the samples automatically in the given directory.

        :return: A list of tuples, representing all the different sample pairs.
        """
        # Only search for samples if samples are not given and there is a regular expression. Otherwise give an error
        # or skip the functions, since samples are already given.
        if not self.samples and self.reg_exp_present:
            if os.path.isdir(self.settings["sample_directory"]):
                path = self.settings["sample_directory"]
                # Shell regular expression for finding the first files.
                reg_expr = self.settings["regular_expression"].format("1")
                # Python regular expression for checking if the 2nd file matches the 1st.
                python_reg_expr = re.compile(self.settings["regular_expression"].replace("*", ".+").format("2"))

                # Try the ls command. And decode if it works.
                try:
                    list_files = subprocess.check_output("ls {path}/{reg_expr}".format(path=path, reg_expr=reg_expr),
                                                         shell=True)
                except subprocess.CalledProcessError as ls_error:
                    sys.exit("Something went wrong with listing the files in the sample directory using the given "
                             "regular expression. Error: {}".format(ls_error.output))
                list_files = list_files.decode("UTF-8")

                # Gather only the filenames from the ls output and search a partner file for them.
                for file in str(list_files).split():
                    file = file.split("/")[-1]
                    for i in range(len(file)):
                        if file[i] == "1":
                            temp_file = file[:i] + "2" + file[i + 1:]
                            if bool(re.match(python_reg_expr, temp_file)):
                                self.samples.append((file, temp_file))

        elif not self.samples and not self.reg_exp_present:
            sys.exit("No samples or regular expression given. Please give either one in settings.txt.")
        else:
            pass

    def next_sample(self):
        """Function which gives the next sample to the wrapper for analysis.

        :return: either a tuple of samples which are to be analyzed next, or "None" which indicates that all samples
        have been analysed
        """
        if len(self.samples) > 0:
            self.current_sample = self.samples.pop(-1)
        else:
            self.current_sample = None

    def run_viromatch(self):
        """Run the actual viromatch command using the parsed settings."""

        command = "docker container run -it " + \
                  "-v {sample_dir}:/data " + \
                  "-v {out_dir}:/outdir " + \
                  "-v {nt_dir}:/nt " + \
                  "-v {nr_dir}:/nr " + \
                  "-v {viral_nt_dir}:/viralfna " + \
                  "-v {viral_nr_dir}:/viralfaa " + \
                  "-v {host_dir}:/host " + \
                  "-v {adaptor_dir}:/adaptor " + \
                  "-v {taxonomy_dir}:/taxonomy " + \
                  "twylie/viromatch:latest viromatch " + \
                  "--smkcores {threads} " + \
                  "--sampleid {current_sample_1} " + \
                  "--input /data/{current_sample_1} /data/{current_sample_2} " + \
                  "--outdir /outdir/{current_sample_1}_python_automated " + \
                  "--nt /nt/{nt_file} " + \
                  "--nr /nr/{nr_file} " + \
                  "--viralfna /viralfna/{viral_nt_file} " + \
                  "--viralfaa /viralfaa/{viral_nr_file} " + \
                  "--host /host/{host_file} " + \
                  "--adaptor /adaptor/{adaptor_file} " + \
                  "--taxid /taxonomy/{taxonomy_file};"

        command = command.format(sample_dir=self.settings["sample_directory"],
                                 out_dir=self.settings["output_directory"],
                                 nt_dir=self.settings["ncbi_nt"], nr_dir=self.settings["ncbi_nr"],
                                 viral_nt_dir=self.settings["viral_nt"], viral_nr_dir=self.settings["viral_nr"],
                                 host_dir=self.settings["host"], adaptor_dir=self.settings["adaptor"],
                                 taxonomy_dir=self.settings["taxonomy"], threads=self.settings["cores"],
                                 current_sample_1=self.current_sample[0], current_sample_2=self.current_sample[1],
                                 nt_file=self.settings["ncbi_nt_file"], nr_file=self.settings["ncbi_nr_file"],
                                 viral_nt_file=self.settings["viral_nt_file"],
                                 viral_nr_file=self.settings["viral_nr_file"],
                                 host_file=self.settings["host_file"], adaptor_file=self.settings["adaptor_file"],
                                 taxonomy_file=self.settings["taxonomy_file"])
        subprocess.run(command, shell=True)

    def viromatch_wrapper(self):
        """The actual wrapper function which runs all the functions needed to run the ViroMatch samples."""
        self.parse_settings()
        self.check_settings()
        self.search_samples()
        print("samples:", self.samples)
        for samples in self.samples:
            (print(samples))
        self.next_sample()
        while self.current_sample:
            print("\nRunning samples {} and {}\n".format(self.current_sample[0], self.current_sample[1]))
            self.run_viromatch()
            self.next_sample()


def test_creation():
    dummy = ViroMatchWrapper()
    return dummy


def test_parse_settings(dummy):
    dummy.parse_settings()
    for settings in dummy.settings:
        print(settings + ":", dummy.settings[settings])
    print("\n", dummy.reg_exp_present, "\n", dummy.samples)


def test_check_settings(dummy):
    dummy.check_settings()


def test_viromatch_wrapper(dummy):
    dummy.viromatch_wrapper()


if __name__ == "__main__":
    ViroMatchWrapper().viromatch_wrapper()
