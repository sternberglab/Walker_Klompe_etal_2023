#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg Laboratory
# Last updated: 2023-01-05
#----------------------------------
"""Assign enrichment scores (fold-change) to code(s) associated with each barcode"""
#----------------------------------
from pickle import TRUE
import pandas as pd
import csv
from decimal import *

def main():
    assign_fold_change_to_codes()

def assign_fold_change_to_codes():
    for rep in ["rep_1", "rep_2", "rep1+rep2"]:
        for library in ["Left", "Right"]:
            # 3 main variables: (1) codes:barcodes dict, (2) fold-changes of all barcodes, (3) mutations directory of codes to assign enrichments.
            codes_barcodes_dir = "0_info/codes_barcodes_" + library + ".csv"
            enrichment_dir = "5_fold-change/" + rep + "/normalized/" + library + ".csv"
            mutation_categories = ["one_bp", "two_bp", "four_bp", "spacing", "truncations", "binding_sites", "extra_binding_site", "all"]
            codes_barcodes_dict = pd.read_csv(codes_barcodes_dir, index_col = 0, squeeze=True).to_dict()
            enrichments = pd.read_csv(enrichment_dir, index_col = 0, squeeze=TRUE)
            LR_enrichments_dict = enrichments["tLR"].to_dict()
            RL_enrichments_dict = enrichments["tRL"].to_dict()
            for mutations in mutation_categories:
                input_codes_dir = "0_info/mutation_codes/" + library + "/" + mutations + ".csv"
                output_enrichments_dir = "6_code_assignment/" + rep + "/" + library + "/" + mutations + ".csv"
                with open(input_codes_dir) as input_codes, open(output_enrichments_dir, 'w') as out:
                    reader = csv.reader(input_codes)
                    next(reader)
                    writer = csv.writer(out)
                    writer.writerow(["code", "barcode", "tLR", "tRL"])
                    for codes in reader:
                        barcode = codes_barcodes_dict[codes[0]]
                        writer.writerow([codes[0], barcode, LR_enrichments_dict[barcode], RL_enrichments_dict[barcode]])
                print(output_enrichments_dir)

main()