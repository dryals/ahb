#get chr number from command args
args = commandArgs(trailingOnly=TRUE)

#get N for each pop and total
setwd("/depot/bharpur/data/projects/ryals/ahb/references")
popN = read.delim("refN.txt", header = F)

#read freqs for chr
setwd(paste0("/scratch/bell/dryals/ahb/aim/chr", args[1]))
#setwd("/home/dylan/Documents/bees/harpurlab/project/diversity/aim")

print(paste(args[1], "loading..."))

freqs0 = read.delim(paste0("chr",args[1], ".popfrq"), sep = "", header = F)

#insure this is the correct order!
pops = c("A", "C", "M", "O")
colnames(freqs0) = c("chr", "pos", pops)
k = 4

#remove any repeat values (two pops with same freq) or NA loci (".")
rows = apply(freqs0, 1, function(L){
  length(unique(L[3:(k+2)])) == k &
    !(any(L == "."))
})
freqs = freqs0[rows,]

#convert to numeric
class(freqs$A) = "numeric"
class(freqs$C) = "numeric"
class(freqs$M) = "numeric"
class(freqs$O) = "numeric"


#calculate all Pij's
  #multiply all pop's freq by N, divide by total N       
freqs$p = as.matrix(freqs[3:(k+2)]) %*% popN$V1[1:k] / popN$V1[(k+1)]

#calculate constant
  #stir = abs(Stirling1(k+1,2))
  stir = 50 
  stirfack = stir/factorial(k)

#protected log function: log(0) = 0
p.log = function(x){
  ifelse(x == 0, 0, log(x))
}

print(paste(args[1], "running..."))

#calculate Ia across all sites L
Ia = apply(freqs, 1, function(L){
  
  #allele frequencies at site
  f = unlist(L[3:(k+2)])
  #matrix for combinatorial products
    #j = 0
  pd1 = sapply(1:k, function(x){prod(f[x] - f[-x])}, simplify = T)
    #since prod(j=1) = -prod(j=0)
  P = c(pd1, -1 * pd1)
  #inner summation: for j =0 and j = 1
    #for biallelic sites: p(j=0) = 1-p(j=1)
  s1 = c(f**k * p.log(f), 
                (1-f)**k * p.log(1-f)) / (k * P)
  #sum across j = 0, j = 1
  S = c(sum(s1[1:k]), sum(s1[(k+1):(2*k)]))
  #outer summation
    #pull pop-wide freq
  pj = unlist(L[k+3])
    #calculate for j=0 and j=1, add to inner summation S
    #sum these for total Ia expression
  sum(c(pj * (-1 * p.log(pj) + 1 - stirfack), 
           (1-pj) * (-1 * p.log(1-pj) + 1 - stirfack)) + S)

})

#include chr and pos
out = cbind(freqs[1:2], Ia)

#write
write.table(out, file = paste0("chr", args[1], ".ia"), row.names = F, quote = F)

print(paste(args[1], "done"))
