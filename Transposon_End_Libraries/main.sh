#!/bin/sh

################################################################################################
#
# Script Name   :   main.sh
# Description   :   Transposon end library analysis
# Author        :   Matt W.G. Walker
# Institution   :   Columbia University, Sternberg Laboratory
# Last updated  :   January 17, 2023
# Requirements  :   All python packages in requirements.txt
#
#################################################################################################

python3 $PWD/samples_analysis.py

python3 $PWD/enrichment_analysis.py

python3 $PWD/uncoupling_analysis.py

python3 $PWD/plot.py