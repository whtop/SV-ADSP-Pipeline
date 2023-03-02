#! /bin/bash

for ((i = 1; i < 23; i++))
do
bcftools filter -e 'INFO/SVTYPE="BND"' chrs_union/${i}.vcf -Oz -o chrs_union/${i}_no_BND.vcf.gz
bcftools index chrs_union/${i}_no_BND.vcf.gz
done
