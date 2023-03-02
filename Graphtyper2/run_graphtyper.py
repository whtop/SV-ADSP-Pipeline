import os
import pandas as pd

chrs = """
1	248,956,422	CM000663.2	NC_000001.11
2	242,193,529	CM000664.2	NC_000002.12
3	198,295,559	CM000665.2	NC_000003.12
4	190,214,555	CM000666.2	NC_000004.12
5	181,538,259	CM000667.2	NC_000005.10
6	170,805,979	CM000668.2	NC_000006.12
7	159,345,973	CM000669.2	NC_000007.14
8	145,138,636	CM000670.2	NC_000008.11
9	138,394,717	CM000671.2	NC_000009.12
10	133,797,422	CM000672.2	NC_000010.11
11	135,086,622	CM000673.2	NC_000011.10
12	133,275,309	CM000674.2	NC_000012.12
13	114,364,328	CM000675.2	NC_000013.11
14	107,043,718	CM000676.2	NC_000014.9
15	101,991,189	CM000677.2	NC_000015.10
16	90,338,345	CM000678.2	NC_000016.10
17	83,257,441	CM000679.2	NC_000017.11
18	80,373,285	CM000680.2	NC_000018.10
19	58,617,616	CM000681.2	NC_000019.10
20	64,444,167	CM000682.2	NC_000020.11
21	46,709,983	CM000683.2	NC_000021.9
22	50,818,468	CM000684.2	NC_000022.11
X	156,040,895	CM000685.2	NC_000023.11
Y	57,227,415	CM000686.2	NC_000024.10
"""

cyto = pd.read_csv("../masks/unused/cytoBand.txt", delim_whitespace=True, header=None)
cyto = cyto[cyto[4]=="acen"]
chrs = [k.split()[:2] for k in chrs.strip().split(sep="\n")]
seg = 1000000

# check centromere location
for i, j in chrs:
    df = cyto.loc[cyto[0]=="chr{}".format(i), [1,2]]
    cs = df.min().min()  # start of centromere
    ce = df.max().max()  # end of centromere
    print(i, cs/seg, ce/seg, (ce - cs)/seg)


output_dir="/home/hui_wang/20210326_sv_cnv/data/graphtyper/output_union"
for i, j in chrs:
    if i=="X" or i =="Y": continue

    j = int(j.replace(",", ""))
    n = int(j/seg)
    df = cyto.loc[cyto[0]=="chr{}".format(i), [1,2]]
    cs = df.min().min()  # start of centromere
    ce = df.max().max()  # end of centromere

    for k in range(n+1):

        st = k*seg + 1
        ed = (k+1)*seg
        if (cs >= ed) or (ce <= st):
            pass # outside of centromere
        elif (cs > st) and (ce > ed):
            continue
            #ed = cs  # intersect with the start of centromere
        elif (cs < st) and (ce < ed):
            continue
            #st = ce  # intersect with the end of centromere
        #elif (cs > st) and (ce < ed): # region includes centromere, this only happens on chromosome Y, so ignore this
        else:
            continue
        region = "chr{}:{}-{}".format(i, st, ed)
        if os.path.exists(os.path.join(output_dir, "chr{}".format(i), "{:09}-{:09}.vcf.gz".format(st, ed))): 
            continue
        # mayge do not split by region, but by short/long arm of chromosome
        # then assign resources depending on the length of short/long arm
        os.system('qsub -cwd -v CHR="{}",REGION="{}" -pe DJ 6 -N {} -l h_vmem=3G run_graphtyper.sh'.format(
            i, region, region.replace(":", "_")
            ))
    