#!/usr/bin/perl

# Alfonso Saera Vila
# 30/11/2018

# script to read fastq files and perform quality preprocessing using bbduk and
# fastqc

use File::Path qw(make_path);

# These must be adapted appropriately
my $BBMap_path="/mnt/518D6BCF3ECC578E/BBMap_38.32/bbmap";
my $FASTQC_path="/mnt/518D6BCF3ECC578E/fastqc_v0.11.8/FastQC";
my $Working_dir_path="/mnt/518D6BCF3ECC578E/RNAseq";
my $Fastq_path="${Working_dir_path}/fastq";
my $result_path="${Working_dir_path}/results";
my $TR_path="${result_path}/tr_fastq";

# create directories
make_path ("$result_path/fastqc") or die "cannot create folder $result_path/fastqc\n";
make_path ("$TR_path/fastqc") or die "cannot create folder $TR_path/fastqc\n";

# Process sample by sample
opendir DIR, $Fastq_path or die "cannot open folder $Fastq_path\n";
while ($file = readdir DIR) {
  if ($file !~ /^\./) { # avoid hidden files and the dirs . and ..
    if ($file =~ /(.+).fastq/) {
      $sample_name=$1;
      print "\nProcessing $sample_name\n";
      # QC before Trimming
      print "\nFASTQC before trimming\n\n";
      system "$FASTQC_path/fastqc $Fastq_path/${sample_name}.fastq -o $result_path/fastqc";
      # Trimming
      print "\nTrimming step\n\n";
      system "$BBMap_path/bbduk.sh in=$Fastq_path/${sample_name}.fastq out=$TR_path/${sample_name}.trim.fastq minlen=25 qtrim=r trimq=10";
      # QC after Trimming
      print "\nFASTQC after trimming\n\n";
      system "$FASTQC_path/fastqc $TR_path/${sample_name}.trim.fastq -o $TR_path/fastqc";
    } else {
      print "$file is not a proper sample file. It will not be used.\n";
    }
  }
}
closedir DIR;
