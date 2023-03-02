import pandas as pd

def normalize(age, norm=True):
    if norm:
        return (age - age.mean())/(age.std())
    else:
        return (age - age.mean())/(age.max() - age.min())

outliers = [i.strip() for i in open("outliers.txt")]
pheno = pd.read_csv("17k_phenotype_pc_der_20220317.csv")
samples = sorted([i for i in set(pheno['sample.id'].values) if i not in outliers])
pheno = pheno[pheno['sample.id'].isin(samples)]  # remove outlier samples
w = open('adsp.fam', 'w')  # same as in extract.py
for _, row in pheno.iterrows():
    w.write("0\t{}\t0\t0\t{}\t{}\n".format(row['sample.id'], row['sex']+1, row['AD']+1))
w.close()

p1 = pheno[['sample.id', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'age', 'sex']].copy()
p1['age'] = p1.age.fillna(p1.age.median())  # 524 without age status
p1.insert(0, "#FID", 0)
# get dummies, using any combinations with count < 5 as reference
# 1. consider Platform + PCR.free X Platform
dumm1 = pd.get_dummies(pheno['Sequencing_Center']) # USE USUHS as reference
dumm2 = pd.get_dummies(pheno['Platform']) # use NovaSeq as reference
dumm3 = pd.get_dummies(pheno['PCR.free'])
combis1 = [['Yes', 'HiSeq 2000/2500'], ['No', 'HiSeq X/Ten'], ['Yes', 'HiSeq X/Ten'],  ['No', 'HiSeq 2000/2500']]
for i in dumm1.columns:
  if i != "USUHS":
    if i == "Illumina":
      p1[i] = dumm1[i].values*dumm3['No']
    else:
      p1[i] = dumm1[i].values
for x,y in combis1:
    name = "{}|{}".format(x,y)
    p1[name] = (dumm3[x]*dumm2[y]).values
    print(name, sum(p1[name]))
colnames = [i.replace(" ", "_") for i in p1.columns]
p1.columns = colnames
for i in ["PC{}".format(i) for i in range(1,6)]:
    p1[i] = normalize(p1[i].values)
p1['age'] = normalize(p1['age'].values)
p1.to_csv("adsp.cov", index=False)
