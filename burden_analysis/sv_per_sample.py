import pandas as pd

# note that only CNVs have high/moderate SVs
vep = pd.read_csv("SVs.csv", sep="|", index_col=0)   # vep annotation of each SV
hi = set(vep.index[vep.impact])  # define high-impact SVs
def get_count(f, svtype, hi):
    svs = {}
    svs_ = {}
    with open(f) as f:
        next(f)
        for i in f:
            _, s, c, st, ed, tp, _, _ = i.split()
            l = int(ed) - int(st)
            if svtype=="CNV":
                if int(tp) > 2:
                    svtype="DUP"
                else:
                    svtype="DEL"
            name = "chr{}:{}-{}:{}".format(c, st, ed, svtype)
            if name in hi:
                svs_.setdefault(s, []).extend([l]*abs(int(tp)-2))
            svs.setdefault(s, []).extend([l]*abs(int(tp)-2))
    return svs, svs_

# Input data should be in Plink CNV format
ccnv, ccnv_ = get_count("adsp_comm_cnv.cnv", "CNV", hi)
rcnv, rcnv_ = get_count("adsp_rare_cnv.cnv", "CNV", hi)
cins, cins_ = get_count("adsp_comm_ins.cnv", "INS", hi)
rins, rins_ = get_count("adsp_rare_ins.cnv", "INS", hi)
cinv, cinv_ = get_count("adsp_comm_inv.cnv", "INV", hi)
rinv, rinv_ = get_count("adsp_rare_inv.cnv", "INV", hi)


w = open("17k_count_hq.txt", 'w')
w.write("sample\trcnv_count\trcnv_len\tccnv_count\tccnv_len\trins_count\trins_len\tcins_count\tcins_len\trinv_count\trinv_len\tcinv_count\tcinv_len")
w.write("\trcnv_count_\trcnv_len_\tccnv_count_\tccnv_len_\n")
for i in [i.split()[1] for i in open("../covariate/adsp.fam")]:
  w.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
      i,
      len(rcnv[i]), sum(rcnv[i]), 
      len(ccnv[i]), sum(ccnv[i]),
      len(rins[i]), sum(rins[i]),
      len(cins[i]), sum(cins[i]),
      len(rinv.get(i, [])), sum(rinv.get(i, [0])), 
      len(cinv.get(i, [])), sum(cinv.get(i, [0])),

      len(rcnv_.get(i, [])), sum(rcnv_.get(i, [0])), 
      len(ccnv_.get(i, [])), sum(ccnv_.get(i, [0])),
      ))
w.close()
