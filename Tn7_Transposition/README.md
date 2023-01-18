# Tn7 Transposition NGS Analysis

These scripts are used to analysis NGS data related to the e EcoTn7 transposition experiments from Walker, Klompe et al., 2023

## Installation

Use a package manager, i.e. [pip] (https://pip.pypa.io/en/stable/), to install the required packages listed in requirements.txt, if you do not have them already installed.

Ensure that NGS data are in the directory, "1_seq_data". The scripts require the following fastq files (all of which may be downloaded from the SRA. BioProject Accession #PRJNA919078). Details about fastq files may be found in the SRA and in the sample_info.xlsx file, located in the directory, "0_info".

A4329_R1_040122.fastq
A4321_R1_040122.fastq
A4372_R1_040122.fastq
A4356_R1_040122.fastq
A4345_R1_040122.fastq
A4349_R1_040122.fastq
A4383_R1_040122.fastq
A4301_R1_040122.fastq
A4318_R1_040122.fastq

## Usage

Execute main.py to analyze the sequencing data.