#! /bin/bash

#Define locations
#CHR="17"
#REGION="chr17:45500000-46500000"

SV_DIR=chrs_union;
GENO_DIR=output_union;

VCF=$SV_DIR/${CHR}_no_BND.vcf.gz;

CRAM_DIR=crams;

#Define files
REF=GRCh38_full_analysis_set_plus_decoy_hla.fa;

SAM=$CRAM_DIR/cram.list;

cd $GENO_DIR;

graphtyper genotype_sv $REF $VCF --max_files_open=32 --sams=$SAM --threads=6 --region=$REGION --output $GENO_DIR --advanced --force_use_input_ref_for_cram_reading;
