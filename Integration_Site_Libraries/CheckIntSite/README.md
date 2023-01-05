# CheckIntSite

CheckIntSite checks the integration site preference for a given 8-bp sequence, using data from the integration site library.

## Installation

Use a package manager, i.e. [pip] (https://pip.pypa.io/en/stable/), to install the required package (pandas), if you do not have this already installed.

## Usage

In terminal, cd to the parent directory "CheckIntSite" and execute CheckIntSite.py

You will be prompted to enter an 8-bp integration site sequence. The first bp should correspond to the nucleotide 43-bp downstream of your target sequence.

Sequence guidelines: the sequence may contain any nucleotide, including those representing multiple possibile nucleotides (i.e. R, Y, K, etc.). For a full list of supported nucleotide symbols, see http://www.hgmd.cf.ac.uk/docs/nuc_lett.html

After entering an 8-bp sequence, the program will check the integration site library data (in "data/") and will output a csv with results to "results/".

If you receive any errors, double check that you have installed the required package (pandas).