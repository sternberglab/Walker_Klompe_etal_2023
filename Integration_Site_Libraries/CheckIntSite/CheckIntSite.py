
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg and Landweber Laboratories
# Last updated: 2022-01-05
#----------------------------------
"""Generate distance distribution prediction for a given 8-mer. See README.md for details"""
#----------------------------------
import os
import csv
import pandas as pd

def main():
    check_int_site()

def check_int_site():
    user_seq_temp = input("\nEnter 8-bp sequence to check integration site preference: ")
    user_seq = user_seq_temp.upper()

    if len(user_seq) != 8:
        print("\nError. Please ensure that you enter an 8-bp sequence.\n")
        exit(1)

    two_nt_dict = {
        "R" : ["C", "T"],
        "Y" : ["G", "A"],
        "K" : ["C", "A"],
        "M" : ["T", "G"],
        "S" : ["A", "T"],
        "W" : ["G", "C"]
    }
    
    three_nt_dict = {
        "B" : "A",
        "D" : "C",
        "H" : "G",
        "V" : "T"
    }

    for target in ["target-A", "target-B"]:
        output_dir = "results/" + target + "/"
        output_dir_counts = output_dir + user_seq_temp + ".csv"
        os.makedirs(output_dir, exist_ok=True)
        with open(output_dir_counts, 'w') as out:
            writer = csv.writer(out)
            writer.writerow(["distance (bp)", "mean_proportion", "end"])
            for end in ["RL", "LR"]:
                input_dir = "data/" + target + "_" + end + ".csv"
                df = pd.read_csv(input_dir)
                for j, v in enumerate(user_seq):
                    if v == "N":
                        continue
                    elif v in two_nt_dict:
                        df = df[(df["N8"].str[j] != two_nt_dict[v][0]) & (df["N8"].str[j] != two_nt_dict[v][1])]
                    elif v in three_nt_dict:
                        df = df[(df["N8"].str[j] != three_nt_dict[v])]
                    else:
                        df = df[(df["N8"].str[j] == v)]
                for i in range(43,57):
                    row = [i, df[str(i)].mean(), end]
                    writer.writerow(row)

        print("\nSuccess. Count table saved as: %s \n" %(output_dir_counts))

main()