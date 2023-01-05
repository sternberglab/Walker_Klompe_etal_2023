#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg Laboratory
# Last updated: 2023-01-05
#----------------------------------
"""Extract and count transposon end library barcodes for input and output samples"""
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
    "210520_Right_pDNA_pre-trans",
    "210520_Left_pDNA_pre-trans",

    "210604_Right_gDNA_tRL",
    "210604_Right_gDNA_tLR", 
    "210604_Left_pDNA_pre-trans", 
    "210604_Left_pDNA_post-trans", 
    "210604_Right_pDNA_pre-trans",
    "210604_Right_pDNA_post-trans",
    "210604_Left_gDNA_tRL", 
    "210604_Left_gDNA_tLR",
    
    "210625_Left_pDNA_post-trans_spike-in_R1", 
    "210625_Left_pDNA_post-trans_spike-in_R2", 
    "210625_Right_pDNA_post-trans_spike-in_R1", 
    "210625_Right_pDNA_post-trans_spike-in_R2",

    "220401_Right_gDNA_tRL",
    "220401_Right_gDNA_tLR",
    "220401_Right_pDNA_post-trans",
    "220401_Left_gDNA_tLR",
    "220401_Left_gDNA_tRL",
    "220401_Left_pDNA_post-trans",

    "220526_Left_gDNA_tRL_R1"
    "220526_Right_gDNA_tRL_R1"
  ]

  for dataset in datasets: 

    def analyze():
      create_info_csv()
      count_reads()
      trim_PBS()
      count_barcodes()
      stats()
      abundance()
      remove_temp_files()

    def create_info_csv():
      output_dir_filtering = info[dataset]['output_dir'] + "info.csv"
      df = pd.DataFrame([info[dataset]['fastq']], columns = ["sample_name"])
      df.to_csv(output_dir_filtering, index=False)
      print("Created info csv at ", output_dir_filtering)

    def count_reads():
      input_dir = "1_seq_data/" + str(info[dataset]['date']) + "/" + info[dataset]['fastq'] + ".fastq"
      count = 0
      for rec in SeqIO.parse(input_dir, "fastq"):
        count += 1
      add_to_info_csv("total_reads", count)
      print("Total reads = ", count)

    def trim_PBS():
      input_dir = "1_seq_data/" + str(info[dataset]['date']) + "/" + info[dataset]['fastq'] + ".fastq"
      output_dir_trimmed = info[dataset]['output_dir'] + "trimmed.fastq"
      PBS_LENGTH = 19
      PBS_INFO = pd.read_csv("0_info/pbs_info.csv", header = None, index_col = 0, squeeze=True).to_dict()

      trimmed = (
          rec[PBS_LENGTH:] for rec in SeqIO.parse(input_dir, "fastq") if rec.seq.startswith(PBS_INFO[info[dataset]['pbs']])
      )
      count_trimmed = SeqIO.write(trimmed, output_dir_trimmed, "fastq")
      add_to_info_csv("reads_with_PBS", count_trimmed)
      print("Reads with PBS = ", count_trimmed)

    def count_barcodes():
      input_dir = info[dataset]['output_dir'] + "trimmed.fastq"
      library_dir = "0_info/library_design/" + info[dataset]['lib'] + ".csv"
      output_dir_temp = info[dataset]['output_dir'] + "counts_temp.csv"
      output_dir_reads_with_barcodes = info[dataset]['output_dir'] + "with_barcodes.fastq"

      barcode_length = 10
      barcode_start = info[dataset]['barcode_start']

      barcode_seq_dict = pd.read_csv(library_dir, header = None, index_col = 1, squeeze=TRUE).to_dict()
      barcode_counter = dict.fromkeys(barcode_seq_dict, 0)

      # Save reads with library barcodes
      if info[dataset]['fastq'] == "A4352_062521_R2" or info[dataset]['fastq'] == "A4347_062521_R2": # These are for the reverse reads!
        reads_with_barcodes = (
          rec for rec in SeqIO.parse(input_dir, "fastq") if rec.seq[barcode_start:(barcode_start+barcode_length)] in barcode_counter
        )
      else: 
        reads_with_barcodes = (
          rec for rec in SeqIO.parse(input_dir, "fastq") if rec.seq[:barcode_length].reverse_complement() in barcode_counter
        )
      count_reads_with_barcodes = SeqIO.write(reads_with_barcodes, output_dir_reads_with_barcodes, "fastq")
      add_to_info_csv("barcodes_with_matches", count_reads_with_barcodes)
      print("barcodes with matches = ", count_reads_with_barcodes)
      
      # Count matches per barcode
      reads_with_barcodes_parser = SeqIO.parse(output_dir_reads_with_barcodes, "fastq") 
      for rec in reads_with_barcodes_parser:
        if info[dataset]['fastq'] == "A4352_062521_R2" or info[dataset]['fastq'] == "A4347_062521_R2": seq_barcode = rec.seq[barcode_start:(barcode_start+barcode_length)]
        else: seq_barcode = rec.seq[:barcode_length].reverse_complement()
        if seq_barcode in barcode_counter:
          barcode_counter[seq_barcode] += 1

      # Save barcode counter as csv
      s = pd.Series(barcode_counter, name="barcodes_count")
      s.index.name = "barcode"
      s.to_csv(output_dir_temp)
      print("output dir = ", output_dir_temp)

    def stats():
      input_dir = info[dataset]['output_dir'] + "counts_temp.csv"
      counts = pd.read_csv(input_dir)

      # MEAN, STD, MIN, MAX
      mean = np.mean(counts.loc[:,"barcodes_count"],0)
      std = np.std(counts.loc[:,"barcodes_count"],0)
      min = np.percentile(counts.loc[:,"barcodes_count"],0)
      max = np.percentile(counts.loc[:,"barcodes_count"],100)

      # 90:10
      ninety = np.percentile(counts.loc[:,"barcodes_count"],90)
      ten = np.percentile(counts.loc[:,"barcodes_count"],10)
      sum_above_ninety = sum(e for e in counts.loc[:,"barcodes_count"] if e >= ninety)
      sum_below_ten = sum(e for e in counts.loc[:,"barcodes_count"] if e <= ten)
      ninety_ten = sum_above_ninety/sum_below_ten

      add_to_info_csv("mean", mean)
      add_to_info_csv("min", min)
      add_to_info_csv("max", max)
      add_to_info_csv("std", std)
      add_to_info_csv("90th %", ninety)
      add_to_info_csv("10th %", ten)
      add_to_info_csv("ninety-ten ratio", ninety_ten)

    def abundance():
      input_dir_counts = info[dataset]['output_dir'] + "counts_temp.csv"
      input_dir_info = info[dataset]['output_dir'] + "info.csv"
      output_dir_counts = info[dataset]['output_dir'] + "counts.csv"
      filtering_info = pd.read_csv(input_dir_info, index_col = "sample_name")
      reads_with_barcodes = filtering_info["barcodes_with_matches"][0]

      with open(input_dir_counts) as barcode_counts, open(output_dir_counts, 'w') as out:
          reader = csv.reader(barcode_counts)
          writer = csv.writer(out)
          next(reader)
          writer.writerow(["barcode", "barcode_count", "barcode_abundance"])
          for row in reader:
              writer.writerow([row[0], row[1], (Decimal(row[1])/reads_with_barcodes*100)])

    def add_to_info_csv(colname, value):
      output_dir_info = info[dataset]['output_dir'] + "info.csv"
      df = pd.read_csv(output_dir_info)
      df.insert(len(df.columns), colname, value)
      df.to_csv(output_dir_info, index=False)

    def remove_temp_files():
      counts_temp_dir = info[dataset]['output_dir'] + "counts_temp.csv"
      os.remove(counts_temp_dir)

    analyze()

main()