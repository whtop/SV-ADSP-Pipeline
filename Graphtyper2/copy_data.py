# copy svimmer merged calls for each sample to my folder
import os

samples = [i.strip() for i in open("samples.txt")]
fdir = "svimmer_result_directory"
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
