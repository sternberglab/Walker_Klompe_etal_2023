#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg Laboratory
# Last updated: 2023-01-05
#----------------------------------
"""Sum barcode counts across samples and biological reps, calculating abundances and enrichment scores"""
#----------------------------------
from pickle import TRUE
from Bio import SeqIO
import pandas as pd
from os import path
import numpy as np
import csv
from decimal import *
import os

def main():
  info = pd.read_excel("0_info/sample_info.xlsx", index_col=0).to_dict('index')

  datasets = [
    "220401_WT",
    "220401_ihfA",
    "220401_ihfB",
    "220401_fis",
    "220401_hupA",
    "220401_hupB",
    "220401_hns",
    "220401_ycbG",
    "220401_neg_control"
  ]

  for dataset in datasets:
    os.makedirs(info[dataset]['output_dir'], exist_ok=True) 

    def analyze():
      create_info_csv()
      count_reads()
      count_library_members()
      abundance()

    def create_info_csv():
      output_dir_filtering = info[dataset]['output_dir'] + "info.csv"
      df = pd.DataFrame([info[dataset]['fastq']], columns = ["sample_name"])
      df.to_csv(output_dir_filtering, index=False)
      print("Created info csv at ", output_dir_filtering)

    def count_reads():
      input_dir = "1_seq_data/" + info[dataset]['fastq'] + ".fastq"
      count = 0
      for rec in SeqIO.parse(input_dir, "fastq"):
        count += 1
      add_to_info_csv("total_reads", count)
      print("Total reads = ", count)

    def add_to_info_csv(colname, value):
      output_dir_info = info[dataset]['output_dir'] + "info.csv"
      df = pd.read_csv(output_dir_info)
      df.insert(len(df.columns), colname, value)
      df.to_csv(output_dir_info, index=False)

    def count_library_members():
      input_dir = "1_seq_data/" + info[dataset]['fastq'] + ".fastq"
      library_dir = "0_info/library_members.csv"
      output_dir_temp = info[dataset]['output_dir'] + "counts_temp.csv"
      output_dir_reads_with_barcodes = info[dataset]['output_dir'] + "with_barcodes.fastq"

      seq_description_dict = pd.read_csv(library_dir, header = None, index_col = 1, squeeze=TRUE).to_dict()
      seq_counter = dict.fromkeys(seq_description_dict, 0)

      # Count matches per barcode
      total_reads_with_matches = 0
      fastq_parser = SeqIO.parse(input_dir, "fastq")
      for rec in fastq_parser:
        record_seq = rec.seq[:75]
        if record_seq in seq_counter:
          total_reads_with_matches += 1
          seq_counter[record_seq] += 1
      add_to_info_csv("total_reads_with_matches", total_reads_with_matches)
      print("barcodes with matches = ", total_reads_with_matches)

      description_dict = pd.read_csv(library_dir, header = None, index_col = 0, squeeze=TRUE).to_dict()
      description_counter = dict.fromkeys(description_dict, 0)
      for description in description_dict:
        description_counter[description] = seq_counter[description_dict[description]]

      # Save barcode counter as csv
      s = pd.Series(description_counter, name="seq_count")
      s.index.name = "seq"
      s.to_csv(output_dir_temp) 
      print("output dir = ", output_dir_temp)

    def abundance():
      input_dir_counts = info[dataset]['output_dir'] + "counts_temp.csv"
      input_dir_info = info[dataset]['output_dir'] + "info.csv"
      output_dir_counts = info[dataset]['output_dir'] + "counts.csv"
      filtering_info = pd.read_csv(input_dir_info, index_col = "sample_name")
      reads_with_barcodes = filtering_info["total_reads_with_matches"][0]

      with open(input_dir_counts) as barcode_counts, open(output_dir_counts, 'w') as out:
          reader = csv.reader(barcode_counts)
          writer = csv.writer(out)
          next(reader)
          writer.writerow(["description", "count", "abundance"])
          for row in reader:
              writer.writerow([row[0], row[1], (Decimal(row[1])/reads_with_barcodes*100)])

    analyze()

main()

