#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg Laboratory
# Last updated: 2023-01-05
#----------------------------------
"""Generate .csv files of enriched sequences to use for weblogos (Fig 3C & Sup Fig 3A)"""
#----------------------------------
import pandas as pd
import csv
import math
from decimal import *

def main():
    all_enriched_seqs()
    weblogo_seqs_target_A()
    weblogo_seqs_target_B()
    sort_enrichments()

def all_enriched_seqs():
    libraries = ["target-B", "target-A"]
    orientations = ["RL", "LR"]
    # threshold
    thresholds = [1, 1.5, 2, 3, 4, 5]   
    for library in libraries:
        for threshold in thresholds:
            for i in range(43,57):
                for orientation in orientations:
                    input_dir = "6_fold-change/" + library + "/" + str(i) + "bp.csv"
                    output_dir = "7_miscellaneous_analysis/A_top_sequences/" + library + "/threshold/" + str(threshold) + "/" + orientation + "_" + str(i) + ".csv"
                    df = pd.read_csv(input_dir, usecols = ["N8", orientation])
                    df = df[df[orientation].notna()]
                    df_thresholded = df[df[orientation]>threshold]
                    df_sorted = df_thresholded.sort_values(orientation, ascending=False)
                    df_sorted.to_csv(output_dir, index=False)
    # summarize thresholds
    for library in libraries:
        for threshold in thresholds:
            for i in range(43,57):
                input_dir_RL = "7_miscellaneous_analysis/A_top_sequences/" + library + "/threshold/" + str(threshold) + "/RL_" + str(i) + ".csv"
                input_dir_LR = "7_miscellaneous_analysis/A_top_sequences/" + library + "/threshold/" + str(threshold) + "/LR_" + str(i) + ".csv"
                output_dir = "7_miscellaneous_analysis/A_top_sequences/" + library + "/threshold/" + str(threshold) + "/summary/" + str(i) + ".csv"
                with open(input_dir_RL) as input_RL, open(input_dir_LR) as input_LR, open(output_dir, 'w') as out:
                    reader_RL = csv.reader(input_RL)
                    reader_LR = csv.reader(input_LR)
                    writer = csv.writer(out)
                    writer.writerow(["orientation", "N8", "enrichment"])
                    next(reader_RL)
                    next(reader_LR)
                    for row in reader_RL:
                        writer.writerow(["RL", row[0], row[1]])
                    for row in reader_LR:
                        writer.writerow(["LR", row[0], row[1]])

def weblogo_seqs_target_A():
    # load info
    info = pd.read_excel("0_info/weblogo_seq_info.xlsx", index_col=0).to_dict('index')
    master_list = []
    # iterate through positions -5, tsd1-5, to +5
    for key in info:
        MAX = 5000 # 5785
        # find max seqs
        max_seqs = []
        num_seqs = math.floor(MAX/(info[key]["draw_from"]*2))
        for orientation in ["RL", "LR"]:
            for i in range(0,info[key]["draw_from"]):
                dir_to_open = info[key]["starting_from"] + i
                input_dir = "8_miscellaneous_analysis/A_top_sequences/target-A/threshold/2/" + orientation + "_" + str(dir_to_open) + ".csv"
                df = pd.read_csv(input_dir)
                max_seqs.append(len(df))
        num_seqs_to_extract = [0] * len(max_seqs)
        counter = 0
        while counter < MAX:
            for j in range(0,len(max_seqs)):
                if counter < MAX:
                    if max_seqs[j] > 0:
                        num_seqs_to_extract[j] += 1
                        max_seqs[j] -= 1
                        counter += 1
                    elif max_seqs[j] <=0:
                        continue
                elif counter >= MAX:
                    break

        nucleotide_list = []
        j = 0
        for orientation in ["RL", "LR"]:
            for i in range(0,info[key]["draw_from"]):
                # open up the RL
                dir_to_open = info[key]["starting_from"] + i
                # print(dir_to_open)
                input_dir = "8_miscellaneous_analysis/A_top_sequences/target-A/threshold/2/" + orientation + "_" + str(dir_to_open) + ".csv"
                df = pd.read_csv(input_dir)
                # take appropriate num of seqs
                if isinstance(key, int):
                    if key>1:
                        nucleotide_list.extend(df["N8"].head(num_seqs_to_extract[j]).str[i+(int(key)-1)].tolist()) # CHANGE I HERE
                    else:
                        nucleotide_list.extend(df["N8"].head(num_seqs_to_extract[j]).str[i].tolist())
                else:
                    nucleotide_list.extend(df["N8"].head(num_seqs_to_extract[j]).str[i].tolist()) # CHANGE I HERE
                j += 1
        master_list.append(nucleotide_list)
    #print(len(master_list))
    output_dir = "8_miscellaneous_analysis/C_weblogo_seqs/target-A/all.csv"
    with open(output_dir, 'w') as out:
        writer = csv.writer(out)
        writer.writerow(["-5", "-4","-3","-2","-1","TSD1","TSD2","TSD3","TSD4","TSD5","+1","+2","+3","+4","+5"])
        for i in range(0,len(master_list[1])):
            row = []
            for j in range(0,len(master_list)):
                row.append(master_list[j][i])
            writer.writerow(row)

def weblogo_seqs_target_B():
    # load info
    info = pd.read_excel("0_info/weblogo_seq_info.xlsx", index_col=0).to_dict('index')
    master_list = []
    # iterate through positions -5, tsd1-5, to +5
    for key in info:
        MAX = 500 #9681
        # find max seqs
        max_seqs = []
        num_seqs = math.floor(MAX/(info[key]["draw_from"]*2))
        for orientation in ["RL", "LR"]:
            for i in range(0,info[key]["draw_from"]):
                dir_to_open = info[key]["starting_from"] + i
                input_dir = "8_miscellaneous_analysis/A_top_sequences/target-B/threshold/2/" + orientation + "_" + str(dir_to_open) + ".csv"
                df = pd.read_csv(input_dir)
                max_seqs.append(len(df))
        # print(key)
        # print(max_seqs)

        num_seqs_to_extract = [0] * len(max_seqs)
        counter = 0
        while counter < MAX:
            for j in range(0,len(max_seqs)):
                if counter < MAX:
                    if max_seqs[j] > 0:
                        num_seqs_to_extract[j] += 1
                        max_seqs[j] -= 1
                        counter += 1
                    elif max_seqs[j] <=0:
                        continue
                elif counter >= MAX:
                    break

        nucleotide_list = []
        j = 0
        for orientation in ["RL", "LR"]:
            for i in range(0,info[key]["draw_from"]):
                # open up the RL
                dir_to_open = info[key]["starting_from"] + i
                # print(dir_to_open)
                input_dir = "8_miscellaneous_analysis/A_top_sequences/target-B/threshold/2/" + orientation + "_" + str(dir_to_open) + ".csv"
                df = pd.read_csv(input_dir)
                # take appropriate num of seqs
                if isinstance(key, int):
                    if key>1:
                        nucleotide_list.extend(df["N8"].head(num_seqs_to_extract[j]).str[i+(int(key)-1)].tolist()) # CHANGE I HERE
                    else:
                        nucleotide_list.extend(df["N8"].head(num_seqs_to_extract[j]).str[i].tolist())
                else:
                    nucleotide_list.extend(df["N8"].head(num_seqs_to_extract[j]).str[i].tolist()) # CHANGE I HERE
                j += 1
        master_list.append(nucleotide_list)
    #print(len(master_list))
    output_dir = "8_miscellaneous_analysis/C_weblogo_seqs/target-B/all.csv"
    with open(output_dir, 'w') as out:
        writer = csv.writer(out)
        writer.writerow(["-5", "-4","-3","-2","-1","TSD1","TSD2","TSD3","TSD4","TSD5","+1","+2","+3","+4","+5"])
        for i in range(0,len(master_list[1])):
            row = []
            for j in range(0,len(master_list)):
                row.append(master_list[j][i])
            writer.writerow(row)

def sort_enrichments():
    libraries = ["target-B", "target-A"]
    orientations = ["RL", "LR"]
    for library in libraries:
        for i in range(43,57):
            for orientation in orientations:
                    input_dir = "6_fold-change/" + library + "/" + str(i) + "bp.csv"
                    output_dir = "7_miscellaneous_analysis/B_sorted_enrichments/" + library + "/" + orientation + "/" + str(i) + ".csv"
                    df = pd.read_csv(input_dir)
                    df.sort_values([orientation], axis=0,  inplace=True)
                    df.dropna(subset=[orientation], inplace=True)
                    if orientation == "RL":
                        df_new = df.drop("LR", axis=1)
                    elif orientation == "LR":
                        df_new = df.drop("RL", axis=1)
                    df_new.to_csv(output_dir, index=False)

main()