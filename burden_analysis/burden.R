# first calculate sv per sample using python script
library(ggplot2)
require(MASS)

burden <- function(df, vari, num_cores) {
  # burden test: linear regression
  # permutation reference: https://stats.stackexchange.com/questions/292597/how-to-do-permutation-test-on-model-coefficients-when-including-an-interaction-t
  # df <- df[df[,vari]!=0,]
  df[,vari] <- (df[,vari] - mean(df[,vari]))/sd(df[,vari])
  model <- sprintf("AD~%s+PC1+PC2+PC3+PC4+PC5+age+sex+Baylor+Broad+GENENTECH+Illumina+NYGC+WashU+Yes.HiSeq_2000.2500+No.HiSeq_X.Ten+Yes.HiSeq_X.Ten+No.HiSeq_2000.2500", vari)
  m <- glm(model, data=df, family=binomial)
  z1 <- summary(m)$coefficients[2,3]
  #pnorm(summary(m1)$coefficients[2,3], lower.tail = F)
  cl <- parallel::makeCluster(num_cores) # not overload your computer
  doParallel::registerDoParallel(cl)
  `%dopar%` <- foreach::`%dopar%`
  permuted <- foreach::foreach(i = 1:10000, .combine=c) %dopar% {
    df1 <- df
    df1[,vari] <- sample(df1[,vari])
    m1 <- glm(model, data=df1, family=binomial)
    summary(m1)$coefficients[vari, 'z value']
  }
  parallel::stopCluster(cl)
  permuted <- sort(permuted)
  permuted_p <- 1-findInterval(z1, sort(permuted))/length(permuted)
  print(permuted_p)
  newlist<- list('data'=df, 'm'=m, 'coeff'=summary(m)$coefficients, 'or'=exp(cbind(coef(m), confint(m))), 'p_permuted'=permuted_p)
  newlist
}

filename = '17k_count_hq.txt'
count <- read.csv(filename, sep="\t")
count$cnv_count <- count$rcnv_count + count$ccnv_count
count$cnv_len <- count$rcnv_len + count$ccnv_len
count$inv_count <- count$rinv_count + count$cinv_count
count$inv_len <- count$rinv_len + count$cinv_len
count$ins_count <- count$rins_count + count$cins_count
count$ins_len <- count$rins_len + count$cins_len
count$cnv_count_ <- count$rcnv_count_ + count$ccnv_count_
count$cnv_len_ <- count$rcnv_len_ + count$ccnv_len_
count$sv_count <- count$cnv_count + count$inv_count + count$ins_count
count$sv_len <- count$cnv_len + count$inv_len + count$ins_len

adspphe <- read.csv("adsp.phe", sep="", header=F)  # phenotype file
adspcov <- read.csv("adsp.cov", sep=",", header=T, stringsAsFactors=T)  # covariate file
adspcov$AD <- factor(adspphe$V3[match(adspcov$sample.id, adspphe$V2)])
df <- merge(count, adspcov, by.x="sample", by.y="sample.id")

results = list()
# analysis of burden for CNV/INS/INV
for (i in 2:dim(count)[2]) {
  name <- colnames(df)[i]
  results[[name]] <- burden(df, name, 10)
}
saveRDS(results, paste0(filename,'.RDS'))

