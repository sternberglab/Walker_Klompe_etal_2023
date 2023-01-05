#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg Laboratory
# Last updated: 2023-01-05
#----------------------------------
"""Plot comparison between biological replicates (Supplemental Figure 1E)"""
#----------------------------------
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

def main():
    compare_reps_output_fold_change()

def compare_reps_output_fold_change():
    libraries = ["Right", "Left"]
    for library in libraries:
        orientations = ["tRL", "tLR"]
        for orientation in orientations:
            input_dir_1 = "5_fold-change/rep_1/normalized/" + library + ".csv"
            input_dir_2 = "5_fold-change/rep_2/normalized/" + library + ".csv"
            output_dir = "7_graphs/compare_reps/output_fold-change/" + library + "_" + orientation + ".pdf"
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

main()