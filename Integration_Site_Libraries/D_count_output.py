#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg Laboratory
# Last updated: 2023-01-05
#----------------------------------
"""Extract and count degenerate sequences from paired-end output reads"""
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
    "220401_Output_target-A_integrated_P1"
    "220401_Output_target-A_integrated_P2",
    "220401_Output_target-B_integrated_P1",
    "220401_Output_target-B_integrated_P2"
  ]

  for dataset in datasets: 

    def analyze():
      count_Ns()

    def count_Ns():
      print("step 1: make a dictionary of R2 ids : sequences")
      description_R2 = dataset + "_R2"
      input_dir_R2 = info[description_R2]['output_dir'] + "trimmed.fastq"
      R2_parser = SeqIO.parse(input_dir_R2, "fastq")
      R2_id_seq_dict = {}
      for rec in R2_parser:
          rec_id = rec.id
          rec_seq = rec.seq
          R2_id_seq_dict[rec_id] = rec_seq

      print("step 2: make nested master dictionary")
      N8_nested_dict = dict()
      for end in ["LE", "RE"]:
        for i in range(43,57): # 13 dictionaries, from 43-55 bp (inclusive)
          N8_file = open("0_info/N8s.csv")
          dict_name = str(i) + "_" + end
          N8_nested_dict[dict_name] = dict()
          for N8_seq in N8_file:
            N8_seq = N8_seq.strip('\n')
            (key) = N8_seq[:8]
            N8_nested_dict[dict_name][key]=0

      print("step 3: define LE & RE transposon sequence, and the post N8 seq")
      RE = "TGTTGATACAACCATAAAAT"
      LE = "TGTTGATGCAACCATAAAGT"

      print("step 4: iterate through each R1 read. Check the TSD and integration distance and compare to R2. If a match, save the N in appropriate dict")
      description_R1 = dataset + "_R1"
      input_dir_R1 = info[description_R1]['output_dir'] + "trimmed.fastq"
      R1_parser = SeqIO.parse(input_dir_R1, "fastq")
      if info[description_R1]['lib'] == "target-A":
        for rec in R1_parser:
          for i in range(0,14):
            R1_seq = rec.seq[i+5:i+25]
            
            if i == 0:
              if R1_seq == RE:
                R1_tsd = rec.seq[0:5]
                if rec.id in R2_id_seq_dict:
                  R2_tsd = R2_id_seq_dict[rec.id][13:18].reverse_complement()
                  if R1_tsd == R2_tsd:
                    try:
                      N_seq = R2_id_seq_dict[rec.id][5:13].reverse_complement()
                      N8_nested_dict["43_RE"][N_seq] += 1
                    except:
                      break
              elif R1_seq == LE:
                R1_tsd = rec.seq[0:5]
                if rec.id in R2_id_seq_dict:
                  R2_tsd = R2_id_seq_dict[rec.id][13:18].reverse_complement()
                  if R1_tsd == R2_tsd:
                    try:
                      N_seq = R2_id_seq_dict[rec.id][5:13].reverse_complement()
                      N8_nested_dict["43_LE"][N_seq] += 1
                    except:
                      break
            
            if 0 < i <= 8:
              if R1_seq == RE:
                R1_tsd = rec.seq[i:i+5]
                if rec.id in R2_id_seq_dict:
                  R2_tsd = R2_id_seq_dict[rec.id][(13-i):(18-i)].reverse_complement()
                  if R1_tsd == R2_tsd:
                    dict_name = str(i+43) + "_RE"
                    try:
                      N_seq = str(rec.seq[5:5+i]) + R2_id_seq_dict[rec.id][5:(13-i)].reverse_complement()
                      N8_nested_dict[dict_name][N_seq] += 1
                    except:
                      break
              elif R1_seq == LE:
                R1_tsd = rec.seq[i:i+5]
                if rec.id in R2_id_seq_dict:
                  R2_tsd = R2_id_seq_dict[rec.id][(13-i):(18-i)].reverse_complement()
                  if R1_tsd == R2_tsd:
                    dict_name = str(i+43) + "_LE"
                    try:
                      N_seq = str(rec.seq[5:5+i]) + R2_id_seq_dict[rec.id][5:(13-i)].reverse_complement()
                      N8_nested_dict[dict_name][N_seq] += 1
                    except:
                      break
            
            if i > 8:
              if R1_seq == RE:
                R1_tsd = rec.seq[i:i+5]
                if rec.id in R2_id_seq_dict:
                  R2_tsd = R2_id_seq_dict[rec.id][(13-i):(18-i)].reverse_complement()
                  if R1_tsd == R2_tsd:
                    dict_name = str(i+43) + "_RE"
                    try:
                      N_seq = str(rec.seq[5:13])
                      N8_nested_dict[dict_name][N_seq] += 1
                    except:
                      break
              elif R1_seq == LE:
                R1_tsd = rec.seq[i:i+5]
                if rec.id in R2_id_seq_dict:
                  R2_tsd = R2_id_seq_dict[rec.id][(13-i):(18-i)].reverse_complement()
                  if R1_tsd == R2_tsd:
                    dict_name = str(i+43) + "_LE"
                    try:
                      N_seq = str(rec.seq[5:13])
                      N8_nested_dict[dict_name][N_seq] += 1
                      break
                    except:
                      break

      elif info[description_R1]['lib'] == "target-B":
        for rec in R1_parser:
          for i in range(0,14):
            sequence = rec.seq[i+5:i+25]

            if i == 0:
              if sequence == RE:
                R1_tsd = rec.seq[0:5]
                if rec.id in R2_id_seq_dict:
                  R2_tsd = R2_id_seq_dict[rec.id][13:18].reverse_complement()
                  if R1_tsd == R2_tsd:
                    try:
                      N_seq = R2_id_seq_dict[rec.id][5:13]
                      N8_nested_dict["56_LE"][N_seq] += 1
                    except:
                      break
              elif sequence == LE:
                R1_tsd = rec.seq[0:5]
                if rec.id in R2_id_seq_dict:
                  R2_tsd = R2_id_seq_dict[rec.id][13:18].reverse_complement()
                  if R1_tsd == R2_tsd:
                    try:
                      N_seq = R2_id_seq_dict[rec.id][5:13]
                      N8_nested_dict["56_RE"][N_seq] += 1
                    except:
                      break

            if 0 < i <=8:
              if sequence == RE:
                # check that R1_tsd=R2_tsd
                R1_tsd = rec.seq[i:i+5]
                if rec.id in R2_id_seq_dict:
                  R2_tsd = R2_id_seq_dict[rec.id][(13-i):(18-i)].reverse_complement()
                  if R1_tsd == R2_tsd:
                    dict_entry_name = str(56-i) + "_LE"
                    try:
                      N_seq = R2_id_seq_dict[rec.id][5:(13-i)] + str(rec.seq[5:(5+i)].reverse_complement())
                      N8_nested_dict[dict_entry_name][N_seq] += 1
                      break
                    except:
                      break
              elif sequence == LE:
                # check that R1_tsd=R2_tsd
                R1_tsd = rec.seq[i:i+5]
                if rec.id in R2_id_seq_dict:
                  R2_tsd = R2_id_seq_dict[rec.id][(13-i):(18-i)].reverse_complement()
                  if R1_tsd == R2_tsd:
                    dict_entry_name = str(56-i) + "_RE"
                    try:
                      N_seq = R2_id_seq_dict[rec.id][5:(13-i)] + str(rec.seq[5:(5+i)].reverse_complement())
                      N8_nested_dict[dict_entry_name][N_seq] += 1
                      break
                    except:
                      break

            if i > 8:
              if sequence == RE:
                R1_tsd = rec.seq[i:i+5]
                if rec.id in R2_id_seq_dict:
                  R2_tsd = R2_id_seq_dict[rec.id][(13-i):(18-i)].reverse_complement()
                  if R1_tsd == R2_tsd:
                    dict_entry_name = str(56-i) + "_LE"
                    try:
                      N8_nested_dict[dict_entry_name][rec.seq[5:13].reverse_complement()] += 1
                      break
                    except:
                      break
              elif sequence == LE:
                R1_tsd = rec.seq[i:i+5]
                if rec.id in R2_id_seq_dict:
                  R2_tsd = R2_id_seq_dict[rec.id][(13-i):(18-i)].reverse_complement()
                  if R1_tsd == R2_tsd:
                    dict_entry_name = str(56-i) + "_RE"
                    try:
                      N8_nested_dict[dict_entry_name][rec.seq[5:13].reverse_complement()] += 1
                      break
                    except:
                      break
        
      print("step 5: save counters as csv")
      output_dir = info[description_R1]['output_dir'][:-3] + "counts/"
      end_convert = {"LE":"LR", "RE":"RL"}
      for end in ["LE", "RE"]:
        for i in range(0,14):
          dict_name_for_csv = str(i+43) + "_" + end
          s = pd.Series(N8_nested_dict[dict_name_for_csv], name="count")
          s.index.name = "N"
          out_dir = output_dir + end_convert[end] + "_" + str(i+43) + "bp.csv"
          s.to_csv(out_dir)

      print("step 6: save summary count totals as csv")
      output_dir = info[description_R1]['output_dir'][:-3] + "counts/counts_summary.csv"
      with open(output_dir, 'w') as out:
        writer = csv.writer(out)
        writer.writerow(["distance", "orientation", "counts"])
        for end in ["LE", "RE"]:
          for i in range(0,14):
            dict_name = str(i+43) + "_" + end
            s = pd.Series(N8_nested_dict[dict_name], name="count")
            sum_s = sum(s)
            row = [str(i+43), end_convert[end], sum_s]
            writer.writerow(row)

    analyze()

main()