library(dplyr)
fname = read.delim("/scratch/bell/dryals/ahb/plink/plink_admix_filename.txt", header = F)
reffam = read.table("refData.txt", header = T, sep = "\t")
  id = read.delim(paste0("/scratch/bell/dryals/ahb/plink/", fname ,".fam"), sep = " ", header = F) %>% 
    select(1)
  colnames(id) = "id"

sup = id %>% left_join(reffam %>% select(id = SRR, lineage))
sup$lineage[is.na(sup$lineage)] = "-"

write.table(sup$lineage, file = paste0("/scratch/bell/dryals/ahb/plink/", fname ,".pop"),
            sep = "\t",
            quote = F, row.names = F, col.names = F)
