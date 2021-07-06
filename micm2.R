bef_aft_EI <- read.csv("original_and_new_EI_calcium2.csv",header=TRUE)

snps_wcx <- unique(bef_aft_EI[,1])

pv_vect <- NULL
st_vect <- NULL
na_counter<-0
for (i in snps_wcx) {
  
  before<-bef_aft_EI[bef_aft_EI[,1] %in% i,][,3]
  
  after<-bef_aft_EI[bef_aft_EI[,1] %in% i,][,4]
  
  wtest<-wilcox.test(before,after,paired = TRUE, exact = FALSE)$p.value
  wstats <- wilcox.test(before,after,paired = TRUE, exact = FALSE)$statistic
  
  if(is.na(wtest)) {
    wtest<-1
    na_counter<-na_counter+1
  }
  pv_vect<-c(pv_vect, wtest)
  st_vect <- c(st_vect, wstats)
  
}

wcx_results <- cbind(snps_wcx, pv_vect,st_vect)

wcx_results <- as.data.frame(wcx_results)
colnames(wcx_results) <- c("SNP", "p value","stats")
write.csv(wcx_results, file = "wcx_snp_calcium.csv", row.names = FALSE)

