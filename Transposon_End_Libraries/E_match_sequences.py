#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg Laboratory
# Last updated: 2023-01-05
#----------------------------------
"""Determine extent to which barcodes are correctly associated with transposon end seqences"""
#----------------------------------
from pickle import TRUE
from Bio import SeqIO
import pandas as pd
import csv
from collections import defaultdict, Counter

def main():
  info = pd.read_excel("0_info/sample_info.xlsx", index_col=0).to_dict('index')

  library_info_dict = {
    "MiSeq-right" : {
      "dataset" : "210625_Right_pDNA_post-trans_spike-in_R1",
      "end_length" : 127,
      "output_dir" : "8_matching_seqs/MiSeq/right/"},
    "MiSeq-left" : {
      "dataset" : "210625_Left_pDNA_post-trans_spike-in_R1",
      "end_length" : 152,
      "output_dir" : "8_matching_seqs/MiSeq/left/"}
  }

  for library in library_info_dict:
    def analyze():
      extract_barcode_and_transposon_end()
      filter_by_phred_score()
      count_matching_sequences()

    def extract_barcode_and_transposon_end():
      print("extracting barcode and transposon end for: " + library)
      dataset = library_info_dict[library]["dataset"]
      input_dir = info[dataset]['output_dir'] + "with_barcodes.fastq"
      barcode_and_end_seq = (
        rec[:10+library_info_dict[library]["end_length"]] for rec in SeqIO.parse(input_dir, "fastq")
      )
      output_dir = library_info_dict[library]['output_dir'] + "barcode_and_end_seqs.fastq"
      count = SeqIO.write(barcode_and_end_seq, output_dir, "fastq")
      print("Extracted " + str(count) + " reads")

    def filter_by_phred_score():
      dataset = library_info_dict[library]["dataset"]
      print("filtering reads for: " + dataset)
      input_dir = library_info_dict[library]['output_dir'] + "barcode_and_end_seqs.fastq"
      high_quality_reads = (
        rec
        for rec in SeqIO.parse(input_dir, "fastq")
        if min(rec.letter_annotations["phred_quality"]) >= 20
      )
      output_dir = library_info_dict[library]['output_dir'] + "high_quality.fastq"
      count = SeqIO.write(high_quality_reads, output_dir, "fastq")
      print("Saved " + str(count) + " reads")

    def count_matching_sequences():
      info = pd.read_excel("0_info/sample_info.xlsx", index_col=0).to_dict('index')

      print("counting matches for: " + library)
      dataset = library_info_dict[library]["dataset"]
      barcode_seq_dict = pd.read_csv("0_info/library_design/" + info[dataset]['lib'] + ".csv", header = None, index_col = 1, squeeze=TRUE).to_dict()
      barcode_matches_seq_counter = dict.fromkeys(barcode_seq_dict,0)
      barcode_mismatches_seq_counter = dict.fromkeys(barcode_seq_dict,0)
      barcode_mismatches_greater_than_2_counter = dict.fromkeys(barcode_seq_dict,0)
      mismatched_sequences = {keys : [] for keys in barcode_seq_dict}
      barcode_all_number_of_mismatches = defaultdict(list)
      input_dir = library_info_dict[library]['output_dir'] + "high_quality.fastq"
      fastq_with_barcodes = SeqIO.parse(input_dir, "fastq") 

      for read in fastq_with_barcodes:
        barcode = read.seq[:10].reverse_complement()
        sequence = read.seq[10:].reverse_complement()
        corresponding_library_seq = barcode_seq_dict[barcode]
        start = len(corresponding_library_seq)-len(sequence)

        # if sequence matches expected
        if sequence == corresponding_library_seq[start:]:
          barcode_matches_seq_counter[barcode] += 1

        # if sequence mismatches
        elif sequence != corresponding_library_seq[start:]:
          barcode_mismatches_seq_counter[barcode] += 1
          number_mismatches = sum(c1!=c2 for c1,c2 in zip(sequence,corresponding_library_seq[start:])) # check # mismatches
          barcode_all_number_of_mismatches[barcode].append(number_mismatches)
          mismatched_sequences[barcode].append(str(sequence))
          if number_mismatches > 2:
            barcode_mismatches_greater_than_2_counter[barcode] += 1

      # save csv
      output_dir = library_info_dict[library]['output_dir'] + "sequence_mismatches.csv"
      with open(output_dir, 'w') as out:
        writer = csv.writer(out)
        writer.writerow(["barcode", "total reads", "# matches", "# mismatches", "average nt-difference", "% matches", "# of unique mismatched seqs", "most common", "% max incorrect / correct", "% max incorrect / total"])
        for key in barcode_matches_seq_counter:
          matches = barcode_matches_seq_counter[key]
          mismatches = barcode_mismatches_seq_counter[key]
          total_reads = barcode_matches_seq_counter[key] + barcode_mismatches_seq_counter[key]
          try: average_number_mismatches = sum(barcode_all_number_of_mismatches[key])/len(barcode_all_number_of_mismatches[key])
          except: average_number_mismatches = 0
          try: percent_matches =  matches/(matches+mismatches)*100
          except: percent_matches = "NA"
          counter_object = Counter(mismatched_sequences[key])
          try: most_common = counter_object.most_common(1)[0][1]
          except: most_common = "NA"
          try: incorrect_divided_by_correct = (most_common / matches)*100
          except: incorrect_divided_by_correct = "NA"
          try: incorrect_divided_by_total = (most_common / total_reads)*100
          except: incorrect_divided_by_total = "NA"
          writer.writerow([key, total_reads, matches, mismatches, average_number_mismatches, percent_matches, len(set(mismatched_sequences[key])), most_common, incorrect_divided_by_correct, incorrect_divided_by_total])

    analyze()

main()