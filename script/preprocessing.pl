#!/usr/bin/perl

# Alfonso Saera Vila
# 26/11/2018

# script to read fastq files and perform quality preprocessing using bbduk and
# fastqc

use File::Path qw(make_path);

# These must be adapted appropriately
my $BBMap_path="/mnt/518D6BCF3ECC578E/BBMap_38.32/bbmap";
my $FASTQC_path="/mnt/518D6BCF3ECC578E/fastqc_v0.11.8/FastQC";
my $Working_dir_path="./";
my $Fastq_path="${Working_dir_path}/fastq";
my $result_path="${Working_dir_path}/results";
my $TR_path="${result_path}/tr_fastq";

# create directories
make_path ("$result_path/fastqc") or die "cannot create folder $result_path/fastqc\n";
make_path ("$TR_path/fastqc") or die "cannot create folder $TR_path/fastqc\n";

# make list of sample names
opendir DIR, $Fastq_path or die "cannot open folder $Fastq_path\n";
while ($file = readdir DIR) {
  if ($file !~ /^\./) { # avoid hidden files and the dirs . and ..
    if ($file =~ /(.+)_\d.fastq/) {
      $sample_names{$1} = 1;
    } else {
      print "$file is not a proper sample name. It will not be used.\n";
    }
  }
}
closedir DIR;

@sample_names = keys %sample_names;

# Process sample by sample
foreach $sample_name (@sample_names) {
  print "\nProcessing $sample_name\n";
  # QC before Trimming
  print "\nFASTQC before trimming\n\n";
  system "$FASTQC_path/fastqc $Fastq_path/${sample_name}_1.fastq $Fastq_path/${sample_name}_2.fastq -o $result_path/fastqc";
  # Trimming
  print "\nTrimming step\n\n";
  system "$BBMap_path/bbduk.sh in1=$Fastq_path/${sample_name}_1.fastq in2=$Fastq_path/${sample_name}_2.fastq out1=$TR_path/${sample_name}_1.fastq out2=$TR_path/${sample_name}_2.fastq minlen=25 qtrim=r trimq=10";
  # QC after Trimming
  print "\nFASTQC after trimming\n\n";
  system "$FASTQC_path/fastqc $TR_path/${sample_name}_1.fastq $TR_path/${sample_name}_2.fastq -o $TR_path/fastqc";
}
