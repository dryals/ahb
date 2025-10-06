#for parallel processing
args = commandArgs(trailingOnly=TRUE)
chrno = args[1]

suppressMessages(library(tidyverse))

#read data
setwd('~/scratch/ahb')
ld = read.delim('plink/preprune.ld', sep = '', header = T) %>% 
  filter(CHR_A == chrno, CHR_B == chrno)
aim = read.delim('aim/aim.ia.txt', sep = '', header = F)
  colnames(aim) = c("chr", "pos", "ia")

#attach ia for A
ld = ld %>% left_join(aim %>%  select(ia_A = ia, CHR_A = chr, BP_A = pos), by = c("CHR_A", "BP_A"))
#... for B
ld = ld %>% left_join(aim %>%  select(ia_B = ia, CHR_B = chr, BP_B = pos), by = c("CHR_B", "BP_B"))
#difference 
ld$diff = ld$ia_A - ld$ia_B
#max ia
ld$mia = (ld$diff > 0) * ld$ia_A + (ld$diff < 0) * ld$ia_B


ld2 = ld %>% arrange(desc(mia)) %>% 
  select(SNP_A, SNP_B, ia_A, ia_B, diff, mia)
toremove = vector("character", 1)

#remove site with lower ld
while(nrow(ld2)>0){
  
  #keep snp with best ia
  if(ld2$diff[1] > 0){
    maxsnp = ld2$SNP_A[1]
  } else{
    maxsnp = ld2$SNP_B[1]
  }
  #identify all snps in ld with this one
  tmpremove = c(ld2$SNP_A[ld2$SNP_B == maxsnp], ld2$SNP_B[ld2$SNP_A == maxsnp])
  #remove all pairs containing those snps
  ld2 = ld2[!((ld2$SNP_A %in% tmpremove) | (ld2$SNP_B %in% tmpremove)),]
  #add to remove list
  toremove = c(toremove, tmpremove)
  rm(tmpremove)
  
}

#write out

write.table(toremove, paste0("aim/chr", chrno, "/LDremove.txt"), quote = F, row.names = F, col.names = F)

print(paste0("chr ", chrno, ": removing ", length(toremove), " sites"))

