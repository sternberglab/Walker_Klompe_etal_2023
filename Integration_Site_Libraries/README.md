# Integration Site Library Analysis

These scripts are used to analysis NGS data related to the integration site libraries from Walker, Klompe et al., 2023

Note that this directory contains a subdirectory, "CheckIntSite", with the script that was used to check the integration site distribution of a given degenerate sequence. This script, CheckIntSite.py, has its own README.md file in the CheckIntSite directory.

## Installation

Use a package manager, i.e. [pip] (https://pip.pypa.io/en/stable/), to install the required packages listed in requirements.txt, if you do not have them already installed.

Ensure that NGS data are in the directory, "1_seq_data". The scripts require the following fastq files (all of which may be downloaded from the SRA. BioProject Accession #PRJNA919078). Details about fastq files may be found in the SRA and in the sample_info.xlsx file, located in the directory, "0_info".

A4378_R1_020322.fastq
A4375_R1_020322.fastq
A4364_R1_040122.fastq
A4319_R1_040122.fastq
A4348_R1_040122.fastq
A4313_R1_040122.fastq
A4364_R2_040122.fastq
A4319_R2_040122.fastq
A4348_R2_040122.fastq
A4313_R2_040122.fastq

## Usage

Simply execute main.sh. This script will run all python scripts in the appropriate order for analysis. 