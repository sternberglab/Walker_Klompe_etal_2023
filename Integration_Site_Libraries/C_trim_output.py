#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg Laboratory
# Last updated: 2023-01-05
#----------------------------------
"""Trim paired-end output reads"""
#----------------------------------
from Bio import SeqIO
import pandas as pd
from decimal import *

def main():
  info = pd.read_excel("0_info/sample_info.xlsx", index_col=0).to_dict('index')

  datasets = [
    "220401_Output_target-A_integrated_P1",
    "220401_Output_target-A_integrated_P2",
    "220401_Output_target-B_integrated_P1",
    "220401_Output_target-B_integrated_P2"
  ]

  for dataset in datasets: 

    def analyze():
      create_info_csv()
      count_reads()
      trim_pre_N8()
    
    def create_info_csv():
      pe = ["R1", "R2"]
      for end in pe:
        description = dataset + "_" + end
        input_dir = info[description]['fastq']
        output_dir_filtering = info[description]['output_dir'] + "info.csv"
        df = pd.DataFrame([input_dir], columns = ["sample_name"])
        df.to_csv(output_dir_filtering, index=False)
        print("Created info csv at ", output_dir_filtering)

    def count_reads():
      pe = ["R1", "R2"]
      for end in pe:
        description = dataset + "_" + end
        input_dir = "1_seq_data/" + str(info[description]['date']) + "/" + info[description]['fastq'] + ".fastq"
        count = 0
        for rec in SeqIO.parse(input_dir, "fastq"):
          count += 1
        add_to_info_csv("total_reads", count, end)
        print("Total reads = ", count)

    def trim_pre_N8():
      pe = ["R1", "R2"]
      for end in pe:
        description = dataset + "_" + end
        input_dir = "1_seq_data/" + str(info[description]['date']) + "/" + info[description]['fastq'] + ".fastq"
        output_dir_trimmed = info[description]['output_dir'] + "trimmed.fastq"
        PRE_N8_SEQ = info[description]['pre_N8_seq']

        trimmed = (
            rec[len(PRE_N8_SEQ)-5:] for rec in SeqIO.parse(input_dir, "fastq") if rec.seq.startswith(PRE_N8_SEQ)
        )
        count_trimmed = SeqIO.write(trimmed, output_dir_trimmed, "fastq")
        print("Reads with pre_N8_seq = ", count_trimmed)
        add_to_info_csv("reads_with_pre_N8_seq", count_trimmed, end)

    def add_to_info_csv(colname, value, end):
      description = dataset + "_" + end
      output_dir_info = info[description]['output_dir'] + "info.csv"
      df = pd.read_csv(output_dir_info)
      df.insert(len(df.columns), colname, value)
      df.to_csv(output_dir_info, index=False)

    analyze()

main()