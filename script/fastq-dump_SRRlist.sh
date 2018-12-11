#! /bin/bash

# Alfonso Saera Vila

# V 1.0
# 26/11/2018

# find https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP065494
# using view all run option
# select (all) desired SRR entries.
# Click accessionlist to generate list file

# Usage:
# Run the script in the destination folder
#      $ ./fastq-dump_SRRlist.sh

# Modify SRA_toolkit_path as needed

files=`cat SRR_Acc_List.txt`
SRA_toolkit_path="/mnt/518D6BCF3ECC578E/sratoolkit.current-ubuntu64/sratoolkit.2.9.2-ubuntu64/bin"

for i in $files
  do
    # $SRA_toolkit_path/fastq-dump --split-3 --gzip --accession $i
    $SRA_toolkit_path/fasterq-dump $i
  done
