#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg Laboratory
# Last updated: 2023-01-05
#----------------------------------
"""Plot the read counts for the input samples (Supplementary Figure 1A)"""
#----------------------------------
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

def main():
    plot_inputs()

def plot_inputs():
    libraries = ["Right", "Left"]
    for library in libraries:
        input_dir = "3_sum_counts/rep_1/" + library + ".csv"
        output_dir = "7_graphs/input_counts/" + library + "_input.pdf"
        counts = pd.read_csv(input_dir, index_col = 0)

        plt.rcParams['figure.dpi'] = 300
        plt.figure()
        ax = sns.displot(counts, x="input", height=3) # hue=library, #kind="kde"
        #ax.figure.set_size_inches(10,10)
        ax.set(xlabel="Reads (#)", ylabel = "Library Members (#)")
        #plt.show()
        plt.savefig(output_dir, transparent=True)
        print(output_dir)

main()