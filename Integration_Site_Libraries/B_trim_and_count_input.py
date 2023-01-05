#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg and Landweber Laboratories
# Last updated: 2022-08-22
#----------------------------------
"""Extract and count degenerate sequences from single-end input reads"""
#----------------------------------
from Bio import SeqIO
import pandas as pd
import numpy as np
import csv
from decimal import *

def main():
  info = pd.read_excel("0_info/sample_info.xlsx", index_col=0).to_dict('index')

  datasets = [
    "220203_input_digested_P1",
    "220203_input_digested_P2"
  ]

  for dataset in datasets: 

    def analyze():
      create_info_csv()
      count_reads()
      trim_pre_N8()
      count_N8s()
      stats()
      abundance()
    
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

    def trim_pre_N8():
      input_dir = "1_seq_data/" + str(info[dataset]['date']) + "/" + info[dataset]['fastq'] + ".fastq"
      output_dir_trimmed = info[dataset]['output_dir'] + "trimmed.fastq"
      PRE_N8_SEQ = info[dataset]['pre_N8_seq']

      trimmed = (
          rec[len(PRE_N8_SEQ):] for rec in SeqIO.parse(input_dir, "fastq") if rec.seq.startswith(PRE_N8_SEQ)
      )
      count_trimmed = SeqIO.write(trimmed, output_dir_trimmed, "fastq")
      add_to_info_csv("reads_with_pre_N8_seq", count_trimmed)
      print("Reads with pre_N8_seq = ", count_trimmed)

    def count_N8s():
      input_dir = info[dataset]['output_dir'] + "trimmed.fastq"
      output_dir = info[dataset]['output_dir'] + "counts_temp.csv"

      # create dictionary to use as counter for degenerate sequences
      N8_counter_dict = dict()
      N8_file = open("0_info/N8s.csv")
      for N8s in N8_file:
        N8s = N8s.strip('\n')
        (key) = N8s
        N8_counter_dict[key] = 0

      # Count matches per degenerate seq
      N8_parser = SeqIO.parse(input_dir, "fastq") 
      for rec in N8_parser:
        N8_seq = rec.seq[:8]
        if N8_seq in N8_counter_dict:
          N8_counter_dict[N8_seq] += 1
      sum_N8_counts = sum(N8_counter_dict.values())
      add_to_info_csv("reads_with_N8", sum_N8_counts)
      print("Reads with N8 = ", sum_N8_counts)

      # Save N8 counter as csv
      s = pd.Series(N8_counter_dict, name="N8_count")
      s.index.name = "N8"
      s.to_csv(output_dir)

    def stats():
      input_dir = info[dataset]['output_dir'] + "counts_temp.csv"
      counts = pd.read_csv(input_dir)

      # MEAN, STD, MIN, MAX
      mean = np.mean(counts.loc[:,"N8_count"],0)
      std = np.std(counts.loc[:,"N8_count"],0)
      min = np.percentile(counts.loc[:,"N8_count"],0)
      max = np.percentile(counts.loc[:,"N8_count"],100)

      # 90:10
      ninety = np.percentile(counts.loc[:,"N8_count"],90)
      ten = np.percentile(counts.loc[:,"N8_count"],10)
      sum_above_ninety = sum(e for e in counts.loc[:,"N8_count"] if e >= ninety)
      sum_below_ten = sum(e for e in counts.loc[:,"N8_count"] if e <= ten)
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
      reads_with_N8s = filtering_info["reads_with_N8"][0]

      with open(input_dir_counts) as N8_counts, open(output_dir_counts, 'w') as out:
          reader = csv.reader(N8_counts)
          writer = csv.writer(out)
          next(reader)
          writer.writerow(["N8", "N8_count", "N8_abundance"])
          for row in reader:
              writer.writerow([row[0], row[1], (Decimal(row[1])/reads_with_N8s*100)])

    def add_to_info_csv(colname, value):
      output_dir_info = info[dataset]['output_dir'] + "info.csv"
      df = pd.read_csv(output_dir_info)
      df.insert(len(df.columns), colname, value)
      df.to_csv(output_dir_info, index=False)

    analyze()

main()