#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg and Landweber Laboratories
# Last updated: 2022-01-05
#----------------------------------
"""Sum output reads and determine the fold-change compared to input"""
#----------------------------------
from Bio.Seq import Seq
import pandas as pd
import csv
import math
from decimal import *
import os

def main():
    sum_output_count_summary_P1_P2()
    sum_counts_P1_P2()
    discard_reads_with_low_representation_in_input()
    abundance()
    log2_fold_change()

def sum_output_count_summary_P1_P2():
    sample_info = pd.read_excel("0_info/sample_info.xlsx", index_col=0).to_dict('index')
    libraries = ["target-A", "target-B"]
    for library in libraries:
        os.makedirs("3_sum_counts/"+library, exist_ok=True)
        output_dir = "3_sum_counts/" + library + "/summary.csv"
        P1 = pd.read_csv("2_sample_counts/220401/output_" + library + "_P1/counts/counts_summary.csv")
        P2 = pd.read_csv("2_sample_counts/220401/output_" + library + "_P2/counts/counts_summary.csv")
        P1_P2 = pd.concat([P1["distance"], P1["orientation"], P2["counts"] + P1["counts"]], axis=1)
        P1_P2.to_csv(output_dir, index=False)

def sum_counts_P1_P2(): 
    sample_info = pd.read_excel("0_info/sample_info.xlsx", index_col=0).to_dict('index')
    libraries = ["target-A", "target-B"]
    for library in libraries:
        for i in range(43,57):
            output_dir = "3_sum_counts/" + library + "/" + str(i) + "bp.csv"
            input_P1 = pd.read_csv("2_sample_counts/220203/input_digested_P1/counts.csv", index_col=0, usecols = ["N8", "N8_count"], squeeze=True).to_dict()
            input_P2 = pd.read_csv("2_sample_counts/220203/input_digested_P2/counts.csv", index_col=0, usecols = ["N8", "N8_count"], squeeze=True).to_dict()
            RL_P1 = pd.read_csv("2_sample_counts/220401/output_" + library + "_P1/counts/RL_" + str(i) + "bp.csv", index_col=0, usecols = ["N", "count"], squeeze=True).to_dict()
            RL_P2 = pd.read_csv("2_sample_counts/220401/output_" + library + "_P2/counts/RL_" + str(i) + "bp.csv", index_col=0, usecols = ["N", "count"], squeeze=True).to_dict()
            LR_P1 = pd.read_csv("2_sample_counts/220401/output_" + library + "_P1/counts/LR_" + str(i) + "bp.csv", index_col=0, usecols = ["N", "count"], squeeze=True).to_dict()
            LR_P2 = pd.read_csv("2_sample_counts/220401/output_" + library + "_P2/counts/LR_" + str(i) + "bp.csv", index_col=0, usecols = ["N", "count"], squeeze=True).to_dict()
            with open(output_dir, 'w') as out:
                writer = csv.writer(out)
                writer.writerow(["N8","input", "RL", "LR"])
                for barcode in input_P1:
                    if library == "target-A":
                        writer.writerow([barcode, input_P1[barcode]+input_P2[barcode], RL_P1[barcode]+RL_P2[barcode], LR_P1[barcode]+LR_P2[barcode]])
                    elif library == "target-B":
                        barcode_seq = Seq(barcode)
                        writer.writerow([barcode, input_P1[barcode_seq.reverse_complement()]+input_P2[barcode_seq.reverse_complement()], RL_P1[barcode]+RL_P2[barcode], LR_P1[barcode]+LR_P2[barcode]])

def discard_reads_with_low_representation_in_input():
    sample_info = pd.read_excel("0_info/sample_info.xlsx", index_col=0).to_dict('index')
    libraries = ["target-A", "target-B"]
    for library in libraries:
        os.makedirs("4_discard_low_represented/"+library, exist_ok=True)
        for i in range(43,57):
            output_dir = "4_discard_low_represented/" + library + "/" + str(i) + "bp.csv"
            input_P1 = pd.read_csv("2_sample_counts/220203/input_digested_P1/counts.csv", index_col=0, usecols = ["N8", "N8_count"], squeeze=True).to_dict()
            input_P2 = pd.read_csv("2_sample_counts/220203/input_digested_P2/counts.csv", index_col=0, usecols = ["N8", "N8_count"], squeeze=True).to_dict()
            RL_P1 = pd.read_csv("2_sample_counts/220401/output_" + library + "_P1/counts/RL_" + str(i) + "bp.csv", index_col=0, usecols = ["N", "count"], squeeze=True).to_dict()
            RL_P2 = pd.read_csv("2_sample_counts/220401/output_" + library + "_P2/counts/RL_" + str(i) + "bp.csv", index_col=0, usecols = ["N", "count"], squeeze=True).to_dict()
            LR_P1 = pd.read_csv("2_sample_counts/220401/output_" + library + "_P1/counts/LR_" + str(i) + "bp.csv", index_col=0, usecols = ["N", "count"], squeeze=True).to_dict()
            LR_P2 = pd.read_csv("2_sample_counts/220401/output_" + library + "_P2/counts/LR_" + str(i) + "bp.csv", index_col=0, usecols = ["N", "count"], squeeze=True).to_dict()
            with open(output_dir, 'w') as out:
                writer = csv.writer(out)
                writer.writerow(["N8","input", "RL", "LR"])
                for barcode in input_P1:
                    if library == "target-A":
                        input_count =input_P1[barcode]+input_P2[barcode] 
                        if input_count>15:
                            writer.writerow([barcode, input_count, RL_P1[barcode]+RL_P2[barcode], LR_P1[barcode]+LR_P2[barcode]])
                    elif library == "target-B":
                        barcode_seq = Seq(barcode)
                        input_count = input_P1[barcode_seq.reverse_complement()]+input_P2[barcode_seq.reverse_complement()]
                        if input_count>15:
                            writer.writerow([barcode, input_count, RL_P1[barcode]+RL_P2[barcode], LR_P1[barcode]+LR_P2[barcode]])

def abundance():
    libraries = ["target-A", "target-B"]
    for library in libraries:
        os.makedirs("5_abundance/"+library, exist_ok=True)
        for i in range(43,57):
            input_dir = "4_discard_low_represented/" + library + "/" + str(i) + "bp.csv"
            output_dir = "5_abundance/" + library + "/" + str(i) + "bp.csv"
            df = pd.read_csv(input_dir)
            input_sum = df["input"].sum()
            RL_sum = df["RL"].sum()
            LR_sum = df["LR"].sum()
            with open(input_dir) as barcode_counts, open(output_dir, 'w') as out:
                reader = csv.reader(barcode_counts)
                writer = csv.writer(out)
                writer.writerow(["N8", "input", "RL", "LR"])
                next(reader)
                for row in reader:
                    writer.writerow([row[0], (Decimal(row[1])/input_sum*100),(Decimal(row[2])/RL_sum*100), (Decimal(row[3])/LR_sum*100)])

def log2_fold_change():
    libraries = ["target-A", "target-B"]
    for library in libraries:
        os.makedirs("6_fold-change/"+library, exist_ok=True)
        for i in range(43,57):
            input_dir = "5_abundance/" + library + "/" + str(i) + "bp.csv"
            output_dir = "6_fold-change/" + library + "/" + str(i) + "bp.csv"
            with open(input_dir) as input_abundance, open(output_dir, 'w') as out:
                reader = csv.reader(input_abundance)
                writer = csv.writer(out)
                next(reader)
                writer.writerow(["N8", "RL", "LR"])
                for rows in reader:
                    input = Decimal(rows[1])
                    tRL = Decimal(rows[2])
                    tLR = Decimal(rows[3])
                    row = [rows[0]]
                    try:
                        row.append(math.log2(tRL/input))
                    except:
                        row.append("NA") # zeros are converted to NAs
                    try:
                        row.append(math.log2(tLR/input))
                    except:
                        row.append("NA")
                    writer.writerow(row)

main()