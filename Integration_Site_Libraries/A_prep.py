#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg Laboratory
# Last updated: 2023-01-05
#----------------------------------
"""Generate all possible degenerate sequences"""
#----------------------------------
import csv
import itertools

def main():
    generate_N8s()

def generate_N8s():
    Ns = 'ACTG'
    output_dir_deg = "0_info/N8s.csv"

    with open(output_dir_deg, 'w') as out:
        write = csv.writer(out)
        for output in itertools.product(Ns, repeat=8):
            deglist = [''.join(output)]
            write.writerow(deglist)

main()