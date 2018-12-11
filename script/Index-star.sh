#! /bin/bash

# Alfonso Saera Vila
# 02/12/2018

# STAR index generation for Arabidopsis

##################
# CREATING INDEX #
##################

# FOLDERS
fastaPATH=/mnt/518D6BCF3ECC578E/indexes/Arabidopsis_thaliana_Ensembl_TAIR10/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta
indexPATH=/mnt/518D6BCF3ECC578E/indexes/Arabidopsis_thaliana_Ensembl_TAIR10/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/index_STAR
GTFpath=/mnt/518D6BCF3ECC578E/indexes/Arabidopsis_thaliana_Ensembl_TAIR10/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Archives/archive-2015-07-17-14-28-46/Genes
STAR_path=/mnt/518D6BCF3ECC578E/STAR/source

# RUN STAR
$STAR_path/STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $indexPATH \
  --genomeFastaFiles $fastaPATH/genome.fa --sjdbGTFfile $GTFpath/genes.gtf
