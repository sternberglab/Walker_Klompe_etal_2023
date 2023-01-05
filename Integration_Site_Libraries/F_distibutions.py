#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg and Landweber Laboratories
# Last updated: 2022-01-05
#----------------------------------
"""Make .csv files with integration distance distributions"""
#----------------------------------

import pandas as pd
import csv
from decimal import *

def main():
    distributions()

def distributions():
    # first get all sequences in each input library 
    target_A_df = pd.read_csv("3_sum_counts/target-A/43bp.csv", usecols = ["N8", "input"])
    target_B_df = pd.read_csv("3_sum_counts/target-B/43bp.csv", usecols = ["N8", "input"])

    input_count_filter = 0 # only use sequences that are represented in the input. 0 = discards sequences that are not represented
    target_A_filter = target_A_df[target_A_df["input"]>input_count_filter] # > value
    target_A_filter_list = target_A_filter["N8"].tolist()
    target_B_filter = target_B_df[target_B_df["input"]>input_count_filter] # > value
    target_B_filter_list = target_B_filter["N8"].tolist()
    Ns_in_both_libraries = set(target_B_filter_list).intersection(target_A_filter_list)

    # iterate through each library
    libraries = ["target-A","target-B"]
    for library in libraries:
        # make dictionaries at each integration position with the read counts for all degenerate sequences at that position
        N8_dict = dict()
        for i in range(43,57):
            input_dir = "3_sum_counts/" + library + "/" + str(i) + "bp.csv"
            dict_name_RL = str(i) + "_RL"
            df_RL = pd.read_csv(input_dir, usecols = ["N8", "RL"], index_col=[0], squeeze=True).to_dict()
            N8_dict[dict_name_RL] = df_RL
            dict_name_LR = str(i) + "_LR"
            df_LR = pd.read_csv(input_dir, usecols = ["N8", "LR"], index_col=[0], squeeze=True).to_dict()
            N8_dict[dict_name_LR] = df_LR
        # print(N8_dict["49_RL"])

        # iterate through each sequence in the intersecting libraries. Save to csv in terms of proportions. Do this for RL and LR for both libraries
        for end in ["RL", "LR"]:
            output_dir = "7_distributions/proportions/unweighted/" + library + "/" + end + ".csv"
            with open(output_dir, 'w') as out:
                writer = csv.writer(out)
                writer.writerow(["N8", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56"])
                for library_member in Ns_in_both_libraries:
                    row = [library_member]
                    sum_all_integrations = 0
                    # get sum across all integration positions
                    for i in range(43,57):
                        dict_name = str(i) + "_" + end
                        sum_all_integrations +=  N8_dict[dict_name][library_member]
                    # calculate proportion of total reads at each integration position, and append to csv
                    if sum_all_integrations != 0:
                        for i in range(43,57):
                            dict_name = str(i) + "_" + end
                            count_value = N8_dict[dict_name][library_member]
                            row.append(Decimal(count_value/sum_all_integrations))
                        writer.writerow(row)
                    else: continue

    # open the newly created csvs and take the average across all columns! done.
    libraries = ["target-A", "target-B"]
    for library in libraries:
        input_dir_RL = "7_distributions/proportions/unweighted/" + library + "/RL.csv"
        input_dir_LR = "7_distributions/proportions/unweighted/" + library + "/LR.csv"
        output_dir_2 = "7_distributions/proportions/unweighted/" + library + "/average.csv"
        df_RL = pd.read_csv(input_dir_RL)
        df_LR = pd.read_csv(input_dir_LR)
        means_RL = ["RL"]
        means_LR = ["LR"]
        for i in range(43,57):
            means_RL.append(df_RL[str(i)].mean())
            means_LR.append(df_LR[str(i)].mean())
        with open (output_dir_2, 'w') as out:
            writer = csv.writer(out)
            writer.writerow(["orientation", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56"])
            writer.writerow(means_RL)
            writer.writerow(means_LR)

main()