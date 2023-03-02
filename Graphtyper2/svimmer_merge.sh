for ((i = 1; i < 23; i++))
do echo "svimmer vcf_list_union chr$i > chrs_union/$i.vcf" | qsub -V -cwd -l h_vmem=60G
done
