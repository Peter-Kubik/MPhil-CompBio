#!/bin/bash

# FG3 assignment 
# pipeline for RNA-seq analysis
# paper: DOI: 10.1038/ncomms8776

# pipeline adopted from https://github.com/eturro

# Step 1a: Index the reference transcript sequences with Bowtie 1 
bowtie-build --offrate 3 Homo_sapiens.GRCh37.70.ref_transcripts.fa Homo_sapiens.GRCh37.70.ref_transcripts 

# loop over all fastq files 
# trim -> align -> annotate alignment -> mmseq expression quantification 

for f in *.fastq
do 
  echo "Analysing file $f ..."
  # Step 1c: Trim out adapter sequences if necessary
  trim_galore -q 15 --stringency 3 -e 0.05 --length 36 --trim1 $f 
  # Step 2a: Align reads with Bowtie 1 (not Bowtie 2) 
  bowtie -a --best --strata -S -m 100 --chunkmbs 256 -p 8 Homo_sapiens.GRCh37.70.ref_transcripts \
  $f | samtools view -F -bS - | \
  samtools sort -n - $f.namesorted
  # Step 3: Map reads to transcript sets 
  bam2hits Homo_sapiens.GRCh37.70.ref_transcripts.fa $f.namesorted.bam > $f.hits
  # Step 4: Obtain expression estimates
  mmseq $f.hits $f.quant
done 

# Step 5: DE analysis using mmdiff 
