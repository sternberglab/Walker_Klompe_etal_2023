#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg Laboratory
# Last updated: 2023-01-05
#----------------------------------
"""Sum barcode counts across samples, then calculate barcode enrichment scores"""
#----------------------------------

from pickle import TRUE
import pandas as pd
import csv
import math
from decimal import *

def main():
    sum_counts()
    abundance()
    enrichment()
    log2_fold_change()
    normalized_log2_fold_change()

def sum_counts(): 
    sample_info = pd.read_excel("0_info/sample_info.xlsx", index_col=0).to_dict('index')
    replicates = ["rep_1", "rep_2"]
    for replicate in replicates:
        libraries = ["Right", "Left"]
        for library in libraries:
            if replicate == "rep_1":
                output_date = "210604"
            if replicate == "rep_2":
                output_date = "220401"
            output_dir = "3_sum_counts/" + replicate + "/" + library + ".csv"
            input_1 = pd.read_csv(sample_info["210520_" + library + "_pDNA_pre-trans"]["output_dir"]+'counts.csv', index_col=0, usecols = ["barcode", "barcode_count"], squeeze=True).to_dict()
            input_2 = pd.read_csv(sample_info["210604_" + library + "_pDNA_pre-trans"]["output_dir"]+'counts.csv', index_col=0, usecols = ["barcode", "barcode_count"], squeeze=True).to_dict()
            tRL = pd.read_csv(sample_info[output_date + "_" + library + "_gDNA_tRL"]["output_dir"]+'counts.csv', index_col=0, usecols = ["barcode", "barcode_count"], squeeze=True).to_dict()
            tLR = pd.read_csv(sample_info[output_date + "_" + library + "_gDNA_tLR"]["output_dir"]+'counts.csv', index_col=0, usecols = ["barcode", "barcode_count"], squeeze=True).to_dict()
            with open(output_dir, 'w') as out:
                writer = csv.writer(out)
                writer.writerow(["barcode","input", "tRL", "tLR"])
                for barcode in input_1:
                    writer.writerow([barcode, input_1[barcode]+input_2[barcode], tRL[barcode], tLR[barcode]])

def abundance():
    replicates = ["rep_1", "rep_2"]
    for replicate in replicates:
        libraries = ["Right", "Left"]
        for library in libraries:
            input_dir = "3_sum_counts/" + replicate + "/" + library + ".csv"
            output_dir = "4_abundance/" + replicate + "/" + library + ".csv"
            df = pd.read_csv(input_dir)
            input_sum = df["input"].sum()
            tRL_sum = df["tRL"].sum()
            tLR_sum = df["tLR"].sum()
            with open(input_dir) as barcode_counts, open(output_dir, 'w') as out:
                reader = csv.reader(barcode_counts)
                writer = csv.writer(out)
                writer.writerow(["barcode", "input", "tRL", "tLR"])
                next(reader)
                for row in reader:
                    writer.writerow([row[0], (Decimal(row[1])/input_sum*100),(Decimal(row[2])/tRL_sum*100), (Decimal(row[3])/tLR_sum*100)])

def enrichment():
    replicates = ["rep_1", "rep_2"]
    for replicate in replicates:
        libraries = ["Right", "Left"]
        for library in libraries:
            input_dir = "4_abundance/" + replicate + "/" + library + ".csv"
            output_dir = "5_fold-change/fold-change/" + replicate + "/not-normalized/" + library + ".csv"
            with open(input_dir) as input_abundance, open(output_dir, 'w') as out:
                reader = csv.reader(input_abundance)
                writer = csv.writer(out)
                next(reader)
                writer.writerow(["barcode", "tRL", "tLR"])
                for rows in reader:
                    input = Decimal(rows[1])
                    tRL = Decimal(rows[2])
                    tLR = Decimal(rows[3])
                    row = [rows[0]]
                    row.append(tRL/input)
                    row.append(tLR/input)
                    writer.writerow(row)
            print("output directory = ", output_dir)

def log2_fold_change():
    replicates = ["rep_1", "rep_2"]
    for replicate in replicates:
        libraries = ["Right", "Left"]
        for library in libraries:
            input_dir = "4_abundance/" + replicate + "/" + library + ".csv"
            output_dir = "5_fold-change/log2(fold-change)/" + replicate + "/not-normalized/" + library + ".csv"
            with open(input_dir) as input_abundance, open(output_dir, 'w') as out:
                reader = csv.reader(input_abundance)
                writer = csv.writer(out)
                next(reader)
                writer.writerow(["barcode", "tRL", "tLR"])
                for rows in reader:
                    input = Decimal(rows[1])
                    tRL = Decimal(rows[2])
                    tLR = Decimal(rows[3])
                    row = [rows[0]]
                    try:
                        row.append(math.log2(tRL/input))
                    except:
                        row.append("-15") # zeros are converted to -10s
                    try:
                        row.append(math.log2(tLR/input))
                    except:
                        row.append("-15")
                    writer.writerow(row)
            print("output directory = ", output_dir)

def normalized_log2_fold_change():
    WT_barcodes = pd.read_excel("0_info/WT_barcodes.xlsx", index_col=0).to_dict('index')
    replicates = ["rep_1", "rep_2"]
    for replicate in replicates:
        libraries = ["Right", "Left"]
        for library in libraries:
            input_dir = "4_abundance/" + replicate + "/" + library + ".csv"
            output_dir = "5_fold-change/" + replicate + "/normalized/" + library + ".csv"
            abundances = pd.read_csv(input_dir, index_col = 0, squeeze=TRUE)

            # determine average WT fold-change
            WT_fold_changes_tRL= []
            WT_fold_changes_tLR = []
            for keys in WT_barcodes[library]:
                input_abundance = abundances["input"][WT_barcodes[library][keys]]
                tRL_abundance = abundances["tRL"][WT_barcodes[library][keys]]
                tLR_abundance = abundances["tLR"][WT_barcodes[library][keys]]
                tRL_fold_change = tRL_abundance/input_abundance
                tLR_fold_change = tLR_abundance/input_abundance
                WT_fold_changes_tRL.append(tRL_fold_change)
                WT_fold_changes_tLR.append(tLR_fold_change)
            average_WT_fold_change_tRL = Decimal(sum(WT_fold_changes_tRL)/len(WT_fold_changes_tRL))
            average_WT_fold_change_tLR = Decimal(sum(WT_fold_changes_tLR)/len(WT_fold_changes_tLR))
            
            # normalize tRL and tLR to the average WT fold-change
            with open(input_dir) as input_abundance, open(output_dir, 'w') as out:
                reader = csv.reader(input_abundance)
                writer = csv.writer(out)
                next(reader)
                writer.writerow(["barcode", "tLR", "tRL"])
                for rows in reader:
                    input_abundance = Decimal(rows[1])
                    tRL_abundance = Decimal(rows[2])
                    tLR_abundance = Decimal(rows[3])
                    row = [rows[0]]
                    try:
                        row.append(math.log2((tLR_abundance/input_abundance)/average_WT_fold_change_tLR))
                    except:
                        row.append("-15") # zeros are converted to -10s
                    try:
                        row.append(math.log2((tRL_abundance/input_abundance)/average_WT_fold_change_tRL))
                    except:
                        row.append("-15")
                    writer.writerow(row)

main()