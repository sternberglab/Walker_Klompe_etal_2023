#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg Laboratory
# Last updated: 2023-01-05
#----------------------------------
"""Plot read counts for input samples (SF1A), comparison b/w replicates (SF1E), % coupled reads (SF1B) and % most abundant incorrect seq (SF1C)"""
#----------------------------------
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os

def main():
    plot_percent_coupled()
    plot_percent_most_abundant_incorrect_seq()
    plot_inputs()
    compare_reps()

def plot_inputs():
    os.makedirs("7_graphs/input_counts/", exist_ok=True)
    libraries = ["Right", "Left"]
    for library in libraries:
        input_dir = "3_sum_counts/rep_1/" + library + ".csv"
        output_dir = "7_graphs/input_counts/" + library + "_input.pdf"
        counts = pd.read_csv(input_dir, index_col = 0)

        plt.rcParams['figure.dpi'] = 300
        plt.figure()
        ax = sns.displot(counts, x="input", height=3)
        ax.set(xlabel="Reads (#)", ylabel = "Library Members (#)")
        plt.savefig(output_dir, transparent=True)
        print(output_dir)

def compare_reps():
    os.makedirs("7_graphs/compare_reps/", exist_ok=True)
    libraries = ["Right", "Left"]
    for library in libraries:
        orientations = ["tRL", "tLR"]
        for orientation in orientations:
            input_dir_1 = "5_enrichment/log2(fold-change)/rep_1/normalized/" + library + ".csv"
            input_dir_2 = "5_enrichment/log2(fold-change)/rep_2/normalized/" + library + ".csv"
            output_dir = "7_graphs/compare_reps/" + library + "_" + orientation + ".pdf"
            counts_rep_1 = pd.read_csv(input_dir_1, index_col = 0)
            counts_rep_2 = pd.read_csv(input_dir_2, index_col = 0)
            counts_rep_1.rename(columns = {'tRL': 'tRL_rep1', 'tLR': 'tLR_rep1'}, inplace=True)
            counts_rep_2.rename(columns = {'tRL': 'tRL_rep2', 'tLR': 'tLR_rep2'}, inplace=True)
            counts = pd.concat([counts_rep_1, counts_rep_2], axis=1)

            # style the plot
            plt.rcParams['figure.dpi'] = 300
            plt.figure()
            plt.rcParams.update({'font.size': 20})
            ax = sns.regplot(data=counts, x=counts[orientation + "_rep1"], y=counts[orientation + "_rep2"], line_kws = {'color':'blue'}, scatter_kws={"color":"black",'s':10}) # binwidth=1 hue=library, #kind="kde"
            ax.figure.set_size_inches(10,10)
            ax.spines['left'].set_position('zero')
            ax.spines['bottom'].set_position('zero')
            xlab = " " #"rep 1 log2 fold-change"
            ylab = " " #rep 2 log2 fold-change"
            ax.set(xlabel=xlab, ylabel = ylab)
            ax.set(xlim=(-18, 8))
            ax.set(ylim=(-18, 8))
            plt.savefig(output_dir, transparent=True)
            print(output_dir)

def plot_percent_coupled():
    os.makedirs("7_graphs/mismatches/percent_correct/", exist_ok=True)
    libraries = ["right", "left"]
    for library in libraries:
        input_dir = "6_uncoupling/" + library + "/sequence_mismatches.csv"
        output_dir = "7_graphs/mismatches/percent_correct/" + library + ".pdf"
        counts = pd.read_csv(input_dir, index_col = 0)
        plt.rcParams['figure.dpi'] = 300
        plt.figure()
        ax = sns.displot(counts, x="% matches", height=3, binwidth = 2) # hue=library, #kind="kde"
        # ax.figure.set_size_inches(10,10)
        ax.set(xlabel="% Coupled reads", ylabel = "Library members (#)")
        ax.set(xlim=(1, 100))
        
        plt.savefig(output_dir, transparent=True)
        print(output_dir)

def plot_percent_most_abundant_incorrect_seq():
    os.makedirs("7_graphs/mismatches/percent_most_common/", exist_ok=True)
    library_info_dict = {
        "MiSeq-right" : {
            "input_dir" : "6_uncoupling/right/"},
        "MiSeq-left" : {
            "input_dir" : "6_uncoupling/left/"}
    }
    for library in library_info_dict:
        input_dir = library_info_dict[library]['input_dir'] + "sequence_mismatches.csv"
        output_dir = "7_graphs/mismatches/percent_most_common/" + library + ".png"
        counts = pd.read_csv(input_dir, index_col = 0)

        plt.rcParams['figure.dpi'] = 300
        plt.figure()
        ax = sns.displot(counts, x="% most abundant incorrect seq", height=3, binwidth=1) # hue=library, #kind="kde"
        ax.set(xlabel="% Most abundance incorrect sequence", ylabel = "Library members (#)")
        ax.set(xlim=(1, 100))
        ax.cla()
        plt.savefig(output_dir, transparent=True)
        print(output_dir)

main()