#! /bin/bash

# Alfonso Saera Vila
# 02/12/2018

# V 1.0


#####################
# TIMING the SCRIPT #
#####################

start=`date +%s`

#############
# Functions #
#############

function tiempo {
  min=$(( $1 / 60 ))
  seconds=$(( $1 - $(( $min * 60 )) ))
  minutes=$min
  hours=0
  if [ "$min" -gt 60 ]
  then
    hours=$(( $min / 60 ))
    minutes=$(( $min - $(( $hours * 60 )) ))
  fi
  echo $hours" h., "$minutes" min., "$seconds" secs."
}

###################
# PATHS and FILES #
###################
STAR_path=/mnt/518D6BCF3ECC578E/STAR/source
featureCounts_path=/mnt/518D6BCF3ECC578E/subread-1.6.3-Linux-x86_64/bin
dataPATH=/mnt/518D6BCF3ECC578E/RNAseq/results/tr_fastq

outDIR=/mnt/518D6BCF3ECC578E/RNAseq/results/STAR_mapping

cores=4

indexPATH=/mnt/518D6BCF3ECC578E/indexes/Arabidopsis_thaliana_Ensembl_TAIR10/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/index_STAR
GTFpath=/mnt/518D6BCF3ECC578E/indexes/Arabidopsis_thaliana_Ensembl_TAIR10/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Archives/archive-2015-07-17-14-28-46/Genes

###############################################################################
# Process samples
###############################################################################
for fastq in ${dataPATH}/*.fastq
do
  sample=${fastq##*/}
  name=${sample%.trim.fastq}
  echo `date +"%T"`" - Processing sample $name"

  # create result folder
  mkdir -p ${outDIR}/${name}

  ##########################
  #      STAR MAPPING      #
  ##########################
  start_4=`date +%s`
  echo `date +"%T"`" - STAR MAPPING"
  cd ${outDIR}/${name}
  ${STAR_path}/STAR --runThreadN $cores \
                    --genomeDir $indexPATH \
                    --readFilesIn $fastq \
                    --outSAMtype BAM SortedByCoordinate
                    # outputs Aligned.sortedByCoord.out.bam
  end_4=`date +%s`

  ###########################
  #  Assign reads to genes  #
  ###########################
  start_5=`date +%s`
  echo `date +"%T"`" - FeatureCounts"
  ${featureCounts_path}/featureCounts -T $cores -a $GTFpath/genes.gtf \
                -o gene_assigned_P Aligned.sortedByCoord.out.bam
  end_5=`date +%s`
  runtime_4=$(( end_4 - start_4 ))
  runtime_5=$(( end_5 - start_5 ))
  echo `date +"%T"`"    STAR MAPPING:"
  echo `date +"%T"`"        "`tiempo $runtime_4`
  echo `date +"%T"`"    Assign reads to genes:"
  echo `date +"%T"`"        "`tiempo $runtime_5`
  echo " "
done

#####################
# TIMING the SCRIPT #
#####################
end=`date +%s`
runtime=$(( end - start ))
# show timing in terminal
echo "$0 runtime:"
echo "    "`tiempo $runtime`
