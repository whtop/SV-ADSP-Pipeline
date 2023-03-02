#! /bin/bash

#Define locations
#CHR="17"
#REGION="chr17:45500000-46500000"

SV_DIR=/home/hui_wang/20210326_sv_cnv/data/graphtyper/chrs_union;
GENO_DIR=/home/hui_wang/20210326_sv_cnv/data/graphtyper/output_union;

VCF=$SV_DIR/${CHR}_no_BND.vcf.gz;

CRAM_DIR=/home/hui_wang/20210326_sv_cnv/data/svi_gra;

#Define files
REF=/mnt/data3/old-master/leew/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa;

SAM=$CRAM_DIR/cram.list;

cd $GENO_DIR;

/home/hui_wang/tools/graphtyper genotype_sv $REF $VCF --max_files_open=32 --sams=$SAM --threads=6 --region=$REGION --output $GENO_DIR --advanced --force_use_input_ref_for_cram_reading;
