get_input <- function(filename) {
  count <- read.csv(filename, sep="\t")
  if ('ccnv_count' %in% colnames(count)){
    count$del_count <- count$rdel_count + count$cdel_count
    count$del_len <- count$rdel_len + count$cdel_len
    count$dup_count <- count$rdup_count + count$cdup_count
    count$dup_len <- count$rdup_len + count$cdup_len
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
  }
  adspphe <- read.csv("../covariate/adsp.phe", sep="", header=F)
  adspcov <- read.csv("../covariate/adsp.cov", sep=",", header=T, stringsAsFactors=T)
  adspcov$AD <- factor(adspphe$V3[match(adspcov$sample.id, adspphe$V2)])
  df <- merge(count, adspcov, by.x="sample", by.y="sample.id")
  df
}

power_analysis <- function(df, vari1, pp, runs=100, num_cores=6) {
  cl <- parallel::makeCluster(num_cores) # not overload your computer
  doParallel::registerDoParallel(cl)
  `%dopar%` <- foreach::`%dopar%`  #在环境中指定dopar是foreach的dopar
  permuted <- foreach::foreach(i = 1:runs, .combine=c) %dopar% {
    df1 <- df
    simuY <- sapply(pp, rbinom, n=1, size=1)
    # Then do logistic regression that regresses Y on X and PC, and check if beta is significant
    df1[, "AD"] <- simuY
    model <- sprintf("AD~%s+PC1+PC2+PC3+PC4+PC5+age+sex+Baylor+Broad+GENENTECH+Illumina+NYGC+WashU+Yes.HiSeq_2000.2500+No.HiSeq_X.Ten+Yes.HiSeq_X.Ten+No.HiSeq_2000.2500", vari1)
    m1 <- glm(model, data=df1, family=binomial)
    summary(m1)$coefficients[vari1, 'Pr(>|z|)']
  }
  parallel::stopCluster(cl)
  pow <- sum(permuted <= 0.05)/runs
  pow
}

inputs <- c("./17k_count_hq.txt", "./17k_count_s_hq.txt", "./17k_count_s_4_hq.txt", "./17k_count_h_hq.txt")
fileConn<-"burden_power.txt"
write("File,Name,OR,Power", fileConn)
for (i in inputs) {
  r<-readRDS(paste0(i, ".RDS"))
  df <- get_input(i)
  if (!grepl('s', i, fixed = TRUE)) {
    varis <- c('del_count', 'dup_count', 'ins_count', 'inv_count')
  } else {
    varis <- c('rdel_count', 'rdup_count', 'rins_count', 'rinv_count')
  }
  for (vari in varis) {
    gamma <- r[[vari]]$coeff[,'Estimate']
    df[,vari] <- (df[,vari] - mean(df[,vari]))/sd(df[,vari])
    df['intercept'] <- 1
    covs <- df[, c("intercept", vari, "PC1", "PC2", "PC3", "PC4", "PC5", "age", "sex", "Baylor", "Broad",
                   "GENENTECH", "Illumina", "NYGC", "WashU", "Yes.HiSeq_2000.2500", 
                   "No.HiSeq_X.Ten", "Yes.HiSeq_X.Ten", "No.HiSeq_2000.2500")]
    for (b in c(1.02, 1.04, 1.06, 1.08, 1.10, 1.12)) {
      beta <- log(b)
      gamma[2] <- beta
      eta <- as.matrix(covs) %*% gamma 
      pp <- exp(eta)/ (1+exp(eta))
      pow <- power_analysis(df, vari, pp, runs=100, num_cores=6)
      write(paste(i, vari, b, pow, sep=","), fileConn, append = TRUE)
    }
  }
}


# plot power
library(ggplot2)
po <- read.csv("./burden_power.txt")
po$`SV Type` <- plyr::mapvalues(po$Name, from=c(
  "del_count", "dup_count", "ins_count", "inv_count", "rdel_count", "rdup_count", "rins_count", "rinv_count"),
  to=c("Deleiton", "Duplication", "Insertion", "Inversion", "Deleiton", "Duplication", "Insertion", "Inversion"))
po$Input <- plyr::mapvalues(po$File, from=inputs,
  to=c("All", "Singletons", "Ultra-rare", "Homozygous"))
ggplot(data=po, aes(x=OR, y=Power, color=`SV Type`)) +
  geom_line() + geom_point()+
  facet_wrap(~Input) +
  theme_minimal()
ggplot(data=po[po$Input=="All",], aes(x=OR, y=Power, color=`SV Type`)) +
  geom_line() + geom_point()+
  theme_minimal()
