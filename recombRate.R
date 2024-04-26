#Dylan Ryals Dec 18 2023 
#Recomb-rate interpolation using Marley map
setwd("/depot/bharpur/data/projects/fikere")
load("myobjects")

library(tidyverse)

snp_data = read.delim("chr_pos_DyData.txt", header = F)
colnames(snp_data) = c('chr', 'phys')


# Calculate genetic position
#try all values of span until we get an alway-increasing fit

marey_tmp = marey_data_chr %>% filter(chr == 2)
snp_tmp = snp_data %>% filter(chr == 2)

#calc marey recombination rate
n = nrow(marey_tmp)
A = marey_tmp[-1, 4:5]
B = marey_tmp[-n, 4:5]
rrate = (A[,2] - B[,2])/((A[,1]-B[,1])/1000000)
marey_tmp$rrate = c(NA, rrate)
#fix first NA
marey_tmp$rrate[1] = marey_tmp$rrate[2]

#for (SPAN in seq(0.1, 1, 0.1)){
  SPAN = 0.4
  snp_tmp$rrate = NA
  marey_fit <- (loess(formula = rrate ~ phys, data = marey_tmp, 
                      na.action = na.exclude, span = SPAN, degree = 1, 
                      control = loess.control(surface = "direct")))
  
  snp_tmp$rrate <- predict(marey_fit, newdata = snp_tmp)
  
  # n = nrow(snp_tmp)
  # A = snp_tmp$rrate[-1]
  # B = snp_tmp$rrate[-n]
  # if(sum((A - B < 0)) == 0) break
# }
# rm(A,B)

snp_samp = snp_tmp[sample(1:nrow(snp_tmp), 100000, replace = F),]
plot(snp_samp$phys, snp_samp$rrate)
points(marey_tmp$phys, marey_tmp$rrate, col = 'red')

# calculate the shift value
shift_value <- min(snp_tmp$gen)
#shift the values
snp_tmp$gen <- snp_tmp$gen - shift_value

#create recombination rate
n = nrow(snp_tmp)
A = snp_tmp[-1, 2:3]
B = snp_tmp[-n, 2:3]

rrate = (A[,2] - B[,2])/((A[,1]-B[,1])/1000000)

snp_tmp$rrate = c(NA, rrate)

#fix first NA
snp_tmp$rrate[1] = snp_tmp$rrate[2]



snp_samp = snp_tmp[sample(1:nrow(snp_tmp), 100000, replace = F),]
plot(snp_samp$gen, snp_samp$rrate) 
points(marey_tmp$gen, marey_tmp$rrate, col = 'red')







snp_data = snp_data %>% select(pposition = phys, rrate, gposition = gen)

# create output file with genetic distance added as the third column
write.table(snp_data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = T)
