# copy svimmer merged calls for each sample from wan-ping's folder to my folder
import os

samples = [i.strip() for i in open("../../rare_CNV_analysis/ms_merged/summary_analysis/samples.txt")]
fdir = "/mnt/data3/old-master/leew/17k_sv/svimmer"
odir = "union"
if not os.path.exists(odir): os.mkdir(odir)

w = open("vcf_list_union", 'w')
for i in os.listdir(fdir):
    if i[:-7] in samples:
        os.system("cp {} {}/{}".format(
            os.path.join(fdir, i), odir, i
        ))
        os.system("bcftools index {}/{}".format(odir, i))
        w.write(os.path.join(odir, i) + "\n")

w.close()
