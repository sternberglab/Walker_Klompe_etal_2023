# Transposon End Library Analysis

These scripts are used to analysis NGS data related to the transposon end libraries from Walker, Klompe et al., 2023

## Installation

Use a package manager, i.e. [pip] (https://pip.pypa.io/en/stable/), to install the required packages listed in requirements.txt, if you do not have them already installed.

Ensure that NGS data are in the directory, "1_seq_data". The scripts require the following fastq files (all of which may be downloaded from the SRA. BioProject Accession #PRJNA919078). Details about fastq files may be found in the SRA and in the sample_info.xlsx file, located in the directory, "0_info".

A2000_060421.fastq
A2001_060421.fastq
A2003_060421.fastq
A2004_060421.fastq
A2006_060421.fastq
A2007_060421.fastq
A4337_040122.fastq
A4338_040122.fastq
A4340_040122.fastq
A4347_062521.fastq
A4352_062521.fastq
A4382_052021.fastq
A4390_052021.fastq
A4395_040122.fastq

## Usage

Simply execute main.sh. This script will run all python scripts in the appropriate order for analysis. 