#!/bin/bash

#SAMPLE=$1

M_VCF=s3://wanglab-play/leew8/17k_sv/manta_vcf/$SAMPLE.manta.diploidSV.vcf.gz
S_VCF=s3://wanglab-play/leew8/17k_sv/smoove_vcf/$SAMPLE\-smoove.genotyped.vcf.gz
WD=/mnt/data3/old-master/leew/17k_sv/scripts/tmp
OUT=/mnt/data3/old-master/leew/17k_sv/svimmer

aws s3 cp $M_VCF $WD
aws s3 cp $S_VCF $WD

bcftools index $WD/$SAMPLE.manta.diploidSV.vcf.gz
bcftools index $WD/$SAMPLE\-smoove.genotyped.vcf.gz

echo $WD/$SAMPLE.manta.diploidSV.vcf.gz > $WD/$SAMPLE.txt
echo $WD/$SAMPLE\-smoove.genotyped.vcf.gz >> $WD/$SAMPLE.txt

/mnt/data3/old-master/leew/tools/svimmer/svimmer $WD/$SAMPLE.txt $(seq -f 'chr%g' 1 22) \
  | bgzip -c > $OUT/$SAMPLE.vcf.gz

rm $WD/$SAMPLE.txt
rm $WD/$SAMPLE.manta.diploidSV.vcf.gz*
rm $WD/$SAMPLE\-smoove.genotyped.vcf.gz*
