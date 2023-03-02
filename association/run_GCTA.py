import os
import pandas as pd


w = open('adsp.phe', 'w')
for i in open("adsp.fam"):
    items = i.split()
    w.write("{}\t{}\t{}\n".format(items[0], items[1], int(items[-1])-1))
w.close()

adsp = pd.read_csv("adsp.cov")

adsp[['#FID', 'sample.id'] + list(adsp.columns[8:])].to_csv("adsp.ccov", index=False, sep="\t", header=True)
cols = ['#FID', 'sample.id', 'age'] + ["PC{}".format(i) for i in range(1,6)]
adsp[cols].to_csv("adsp.qcov", sep="\t", index=False, header=True)

#os.system("gcta --bfile adsp --grm-sparse ../data/snp/adsp_sp --fastGWA-mlm --pheno adsp.phe --qcovar adsp.qcov --covar adsp.ccov --geno 1 --maf 0 --thread-num 12 --out adsp")
w = open("run.sh", 'w')
w.write("""#!/bin/bash
#
#SBATCH --job-name=gcta
#SBATCH -o gcta.out
#SBATCH -e gcta.err
#SBATCH --ntasks=10
#SBATCH -N 1
#SBATCH --qos=high
#SBATCH --mem=50g
#
#
gcta --bfile adsp --grm-sparse ../data/snp/adsp_sp --fastGWA-mlm --pheno adsp.phe --qcovar adsp.qcov --covar adsp.ccov --geno 1 --maf 0 --thread-num 10 --out adsp
""")
w.close()
os.system("chmod a+x run.sh")
os.system("sbatch ./run.sh")
