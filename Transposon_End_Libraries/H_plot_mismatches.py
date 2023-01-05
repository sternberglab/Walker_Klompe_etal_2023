#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg Laboratory
# Last updated: 2023-01-05
#----------------------------------
"""Plot % coupled reads (Sup Fig 1B) and % most abundant incorrect seq (Sup Fig 1C)"""
#----------------------------------
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

def main():
    plot_percent_coupled_reads()
    plot_percent_most_abundant_incorrect_seq()

def plot_percent_coupled_reads():
    libraries = ["right", "left"]
    sequencings = ["MiSeq"]
    for library in libraries:
        for sequencing in sequencings:
            input_dir = "8_matching_seqs/" + sequencing + "/" + library + "/sequence_mismatches_high_quality.csv"
            output_dir = "7_graphs/mismatches/percents/percent_correct/" + sequencing + "/high_quality/" + library + ".pdf"
            counts = pd.read_csv(input_dir, index_col = 0)
            plt.rcParams['figure.dpi'] = 300
            plt.figure()
            ax = sns.displot(counts, x="% matches", height=3, binwidth = 2) # hue=library, #kind="kde"
            #ax.figure.set_size_inches(10,10)
            ax.set(xlabel="% Coupled reads", ylabel = "Library members (#)")
            ax.set(xlim=(1, 100))
            plt.savefig(output_dir, transparent=True)
            print(output_dir)

def plot_percent_most_abundant_incorrect_seq():
    library_info_dict = {
        "MiSeq-right" : {
        "input_dir" : "8_matching_seqs/MiSeq/right/"},
        "MiSeq-left" : {
        "input_dir" : "8_matching_seqs/MiSeq/left/"}
    }
    for library in library_info_dict:
        input_dir = library_info_dict[library]['input_dir'] + "sequence_mismatches_high_quality.csv"
        output_dir = "7_graphs/mismatches/percent_most_common/high_quality/" + library + ".png"
        counts = pd.read_csv(input_dir, index_col = 0)

        plt.rcParams['figure.dpi'] = 300
        plt.figure()
        ax = sns.displot(counts, x="% Most abundant incorrect sequence", height=3, binwidth=1) # hue=library, #kind="kde"
        #ax.figure.set_size_inches(10,10)
        ax.set(xlabel="max % most common", ylabel = "Library members (#)")
        ax.set(xlim=(1, 100))
        plt.savefig(output_dir, transparent=True)
        print(output_dir)

main()