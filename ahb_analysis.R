library(purrr)
library(tidyverse)
library(reshape2)
library(readxl)
select = dplyr::select

setwd("/home/dylan/Documents/bees/harpurlab/project/popgen/ahb")


#choose references
#####

#reference ID's on vcf
reffam = read.table("references/refData.txt", header = T, sep = "\t") %>%
  filter(lineage != "ACER")
  #harpurPNAS samples
  reffam$harpur = grepl("^.[0-9]{3}$", reffam$DOGANTID)

#additional sample info
ref.info = read_excel("references/sciadv.abj2151_data_s1_to_s3 4.xlsx", 
                     sheet = "Data S1", range = "A2:P268")
  include = ref.info %>% filter(`Subspseices Assignment / Exclusions` != "Excluded")
  
#just keep included samples
  reffam = reffam %>% filter(DOGANTID %in% include$`Sample ID` | harpur)
  
  
reffam %>% group_by(lineage) %>% summarise(n= n())


#select ~20 references of each lineage
set.seed(2024)
balanced = reffam %>% filter(F)
rand2 = sample(c(1:sum(reffam$lineage == "O")), 20, replace = F)
rand3 = sample(c(1:sum(reffam$lineage == "A")), 20, replace = F)

balanced = rbind(reffam %>% filter(lineage == "M"),
                 reffam %>% filter(lineage == "C"),
                 (reffam %>% filter(lineage == "O"))[rand2,],
                 (reffam %>% filter(lineage == "A"))[rand3,])

#write out
# write.table(balanced$SRR, file = "references/all_refs.txt",
#             quote = F, row.names = F, col.names = F)
# linu = unique(balanced$lineage)
# #write lists
# for(i in 1:length(linu)){
#   
#   p = balanced$SRR[which(balanced$lineage == linu[i])]
#   write.table(p, file = paste0("references/", linu[i], ".txt"),
#               quote = F, col.names = F, row.names = F)
#   notp = balanced %>% filter(!(SRR %in% p)) %>% select(SRR)
#   write.table(notp$SRR, file = paste0("references/not", linu[i], ".txt"),
#               quote = F, col.names = F, row.names = F)
#   
# }

#####

#read mitotypes from first run, make calls
#####
  mito.raw = read.delim("data/mitotype3.out", header = F, sep = "")
  mito = mito.raw[1,]
  mito[1,] = NA
  mito$sampleid = NA
  mito.runs = rle(is.na(mito.raw$V3))
  
  S = 0
  for(i in 1:length(mito.runs[[1]])){
    #sum of rows so far
    S = S + mito.runs$lengths[i]
    #if this is a run of sample names
    if(mito.runs$values[i]){
      #grab rows under the run (blast data)
      tmp = mito.raw[(S+1):(S+mito.runs$lengths[(i+1)]),]
      #give them correct sampleid (last row in run)
      tmp$sampleid = mito.raw[S,1]
      #map to output
      mito = rbind(mito, tmp)
    }
  }
  mito$sampleid = gsub("sampleID[:]", "", mito$sampleid)
  mito = mito[-1,]
  
  #make calls
  usamples = unique(mito$sampleid)
  callmito = function(s){
    #pull best hits based on bitscore
    x = mito %>% filter(sampleid == s) %>% 
      group_by(V2, V12) %>% slice(1) %>%
      select(sampleid, V2, V3, V4, V11, V12) %>% arrange(V11)
    #select min E-values, settling ties by accuracy
    minE = min(x$V11)
    x = x %>% filter(V11 == minE) %>% arrange(desc(V3))
    #report possible calls in order of accuracy
    possible = paste(unique(x$V2), collapse = "/")
    #report top call
    call = x$V2[1]
    #report max bitscore
    bitscore = max(x$V12)
    return(data.frame(s, call, bitscore, possible))
  }

  mito.calls = map_dfr(usamples, callmito)
#####

  
  
#test against 3rd party mitotyping 
  
  ahb = read.csv("AHB_meta2.csv", header = T)
  ahb = ahb %>% left_join(mito.calls %>% rename(vcfid = s))
  #concordance with FL
  ahb %>% filter(loc == "FL", Amito == 1, call != "A1e")
  ahb %>% filter(loc == "FL", Amito == 0, call == "A1e")
  
  
  #who's missing
  ahb %>% filter(is.na(call)) %>% select(vcfid)
  
#read in larger sample set
  ahb = read.csv("AHB_meta3.csv", header = T)
  ahb = ahb %>% left_join(mito.calls %>% rename(vcfid = s))
  
  #fix some id issues
    #az: only one bee per colony
    ahb$id2[ahb$state == "AZ"] = 
      ahb$id[ahb$state == "AZ"]
    #pa -these ids are really messed up, hopefully this will help some
    ahb$id2[ahb$state == "PA"] = 
      gsub("-[0-9]{1,2}-[0-9]{5,}", "",ahb$id2[ahb$state == "PA"])
    #tx - count all colonies sampled over multiple years as one colony
    ahb$id2[ahb$state == "TX"] = 
      gsub("-20[0-9]{2}", "",ahb$id2[ahb$state == "TX"])
    #IN - account for SP collections
    ahb$id2[ahb$state == "IN" & !grepl("202[0-9]Q|IN23", ahb$id2)] = 
      ahb$id[ahb$state == "IN" & !grepl("202[0-9]Q|IN23", ahb$id2)]
      
  
  #verify bees from the same colony have the same mitotype
  ahb %>% group_by(id2) %>% filter(!is.na(call)) %>%
    summarise(n = n(), umito = length(unique(call))) %>% 
    filter(umito>1)
  #these may be due to low bitscores... consider filtering ...
  

#test admixture panel on references
#####
  
  refadmix = cbind( read.delim("data/reference.maf05.4.Q", header = F, sep =""),
    read.delim("data/reference.maf05.fam", sep = "", header = F) %>% select(oldid = V1)) %>%
    left_join(reffam %>% select(lineage, oldid = SRR))
  
  refadmix.melt = refadmix %>% melt(id.vars = c("oldid", "lineage"),
                                    measure.vars = c("V1", "V2", "V3", "V4"),
                                    variable.name = "fam")
  #bar chart
  ggplot(data = refadmix.melt) +
    geom_bar(aes(x = oldid, y = value, fill = fam), 
             stat='identity', width = 1) +
    facet_grid(cols = vars(lineage), scales = "free_x") +
    scale_fill_brewer(palette = "PRGn") +
    labs(x = "sample", y = "Proportion of Genome", fill = "Lineage",
         title = "maf 0.05") + 
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank())

  
  
#####
  
 
#analysis
#####
  
#set up model 
  ahb$AmitoBin = ifelse(ahb$call == "A1e", 1, 0)
  
  #output this for admixture analysis
    #possible different group of samples
  ahb5 = ahb %>%
    group_by(state, id2) %>% arrange(desc(id)) %>% slice(1) %>%
    select(vcfid, oldid, id, id2, state, A, AmitoBin, call)
 # write.csv(ahb4, file = "AHBmeta4.csv", row.names = F, quote = F)


#read in AHB4
ahb4 = read.csv("AHBmeta4.csv")
ahb4 = ahb4 %>% left_join(mito.calls %>% select(vcfid = s, call))
ahb4$AmitoBin = ifelse(ahb4$call == "A1e", 1, 0)
    
 
#connect to most-recent admix run
  fam = read.delim("data/admix.maf05.fam", sep = "", header = F) %>% select(oldid = V1)
  sum(!grepl("SRR", fam$oldid))
  admix = read.delim("data/admix.maf05.4.Q", header = F, sep ="")
  admix = cbind(admix, fam) 
  
  #verify lineage identity
  lins = admix %>% left_join(reffam %>% select(oldid = SRR, lineage)) %>% 
    filter(grepl ("SRR", oldid)) %>%
    pivot_longer(starts_with("V")) %>%
    group_by(name, lineage) %>%
    mutate(m = mean(value)) %>% ungroup %>% 
    group_by(name) %>% arrange(desc(m)) %>% slice(1) %>% select(lineage, name)
  
  lins
  
  #rename
  admix = admix %>% rename(A = V4, M = V1, C = V2, O = V3)
  
  ahb.plot =  ahb4 %>% left_join(admix %>% select(A, gencove_id = oldid))
  
  
  #establish populations
  ahb.plot$pop = ahb.plot$state
  ahb.plot$pop[grepl("NM|AZ", ahb.plot$pop)] = "AZ"
  ahb.plot$pop = factor(ahb.plot$pop,
                        levels = c("IN", "PA", "FL", "Jamaica","TX", "AZ"))
  
#create model
  
  ahb.moddat = ahb.plot %>% filter(!is.na(AmitoBin))
  
  ahbmod = glm(AmitoBin ~ pop + A, family = binomial, data = ahb.moddat)
  summary(ahbmod)
  
    #fit some values
        newdat = map_dfr(unique(ahb.plot$pop), function(x){
                  fakeA = seq(min(ahb.plot$A[ahb.plot$pop == x]), 
                        max(ahb.plot$A[ahb.plot$pop == x]), len = 100)
          return(data.frame(pop = x, A = fakeA))
        })
        newdat$AmitoBin = predict(ahbmod, newdat, type = 'response')
      
  #model fit
  with(ahbmod, pchisq(null.deviance - deviance, 
                      df.null - df.residual, lower.tail = FALSE))
  
  #anova N vs S.
    ahb.ano = ahb.plot
    ahb.ano$globe = ifelse(grepl("IN|PA", ahb.ano$pop), "N", "S")
    globelm = lm(A ~ globe, data = ahb.ano)
    anova(globelm)
    #distributions
    # ggplot(ahb.ano, aes(x = globe, y = A)) + geom_boxplot()
    # ggplot(ahb.ano, aes(x = A)) + geom_histogram() + facet_wrap(facets = vars(globe))
  
#summaries
  #n by state
  ahb.plot %>% group_by(pop) %>% summarise(n=n())
  #n by pop (filtered)
  ahb.plot %>% group_by(pop) %>% summarise(n=n())
  #n each mitotype
  ahb.plot %>% group_by(call) %>% summarise(n=n())
  #range of A calues
  ahb.m = ahb.plot %>% group_by(pop) %>% 
    summarise(r1 = range(A)[1], r2 = range(A)[2], m = mean(A))
  ahb.m
  
  #modeled threshold values of A
  thresh = newdat %>% group_by(pop) %>% mutate(t = abs(0.5 - AmitoBin)) %>%
    arrange(t) %>% slice(1) %>% filter(AmitoBin > 0.4 & AmitoBin < 0.6)
  thresh
  
  #correct calls
  sum( (ahbmod$fitted.values > 0.5 & ahb.moddat$AmitoBin == 1) |
         (ahbmod$fitted.values < 0.5 & ahb.moddat$AmitoBin == 0)) / nrow(ahb.moddat)
  #false positives
  sum( (ahbmod$fitted.values > 0.5 & ahb.moddat$AmitoBin == 0)  )
  #false negatives
  sum( (ahbmod$fitted.values < 0.5 & ahb.moddat$AmitoBin == 1)  )
  
  
#color pallette 
  plot.colors = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")
  

#plot all pops
  ggplot(ahb.plot, aes(x = A, y = AmitoBin, color = pop)) + 
    geom_point(size = 2.5, alpha = 0.5) +
    geom_line(data = newdat, color = 'black', linetype = 2) +
    facet_grid(rows = vars(pop)) + theme_bw() + 
    #plot thresholds
    #geom_vline(data = thresh, aes(xintercept = A), linetype = 2) +
    #plot means
    #geom_vline(data = ahb.m, aes(xintercept = m), linetype = 3)+
    #colors
    scale_color_manual(values = plot.colors)+
    #labels etc
    labs(x = "Estimated A-lineage Ancestry", y = "A-lineage Mitochondrion", 
         color = "State") +
    scale_y_continuous(breaks = c(0,1), 
                       labels = c("False", "True"), limits = c(-.2, 1.2))+
    theme(legend.position = "none")
  
#admixture bar chart
  
  ahb.bar = admix %>% left_join(ahb.plot %>% select(oldid = gencove_id, pop))
  
  #add bar id
  ahb.bar = ahb.bar %>% arrange(pop, A)
  ahb.bar$barid = (1:nrow(ahb.bar))
  
  admix.melt = ahb.bar %>% melt(id.vars = c("oldid", "barid", "pop"), 
                         measure.vars = c("A", "M", "C", "O"),
                         variable.name = "fam")
  
  #bar chart
  ggplot(data = admix.melt %>% filter(!is.na(pop))) +
    geom_bar(aes(x = as.factor(barid), y = value, fill = fam), 
             stat='identity', width = 1) +
    facet_grid(cols = vars(pop), scales = "free_x") +
    scale_fill_brewer(palette = "PRGn") +
    labs(x = "sample", y = "Proportion of Genome", fill = "Lineage") + 
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank())
  
#pca
  pca = read.delim("data/samps.maf05.eigenvec", header = F, sep = " ") %>%
    select(1:5)
  colnames(pca) = c("fam", "oldid", "PC1", "PC2", "PC3")
  pca = pca %>% 
    left_join(ahb.plot %>% select(oldid = gencove_id, individual_id, colony_id, state, call, pop)) %>%
    left_join(admix %>% select(oldid, A, M, C, O), by = 'oldid')
  
  pca$pop = factor(pca$pop,
                        levels = c("IN", "PA", "FL", "Jamaica","TX", "AZ"))
  
  #% explained
  eig = read.delim("data/samps.maf05.eigenval", header = F) %>% filter(V1 > 0)
  eig$x = 1:nrow(eig)
    #ggplot(eig, aes(y = V1, x = x)) + geom_point()
  
  PC1ev = round(eig$V1[1] / sum(eig$V1) * 100, 1)
  PC2ev = round(eig$V1[2] / sum(eig$V1) * 100, 1)
  
#PCA by population
  ggplot(pca, aes(x = PC1, y = PC2, color = pop)) + geom_point(alpha = 0.5, size = 3) +
    theme_bw() +
    scale_color_manual(values = plot.colors) +
    labs(color = "Population",
         x = paste0("PC1 (", PC1ev, "%)"),
         y = paste0("PC2 (", PC2ev, "%)"))
  
#PCA by lineage
  # pcapair = pca %>% select(PC1, PC2, PC3, A, M, C, O)
  # pairs(pcapair, upper.panel = NULL)
  
  
#all pca
  #load data
  #ahb4 = read.csv("AHBmeta4.csv") %>% rename(oldid = gencove_id)
  allpca = read.delim("data/all.maf05.eigenvec", header = F, sep = " ") %>%
    select(1:5)
  colnames(allpca) = c("fam", "oldid", "PC1", "PC2", "PC3")
  
  #% explained
  eig = read.delim("data/all.maf05.eigenval", header = F) %>% filter(V1 > 0)
  PCev = c(round(eig$V1[1] / sum(eig$V1) * 100, 1),
           round(eig$V1[2] / sum(eig$V1) * 100, 1))
  PCev
  
  #attache lineage names
  reflins = read.delim("/home/dylan/Documents/bees/harpurlab/project/popgen/admixResults/fullref/refData.txt")
  allpca = allpca %>% left_join(reflins %>% select(oldid = SRR, lineage, country))
  allpca$lineage[is.na(allpca$lineage)] = "admixed"
  allpca = allpca %>% left_join(ahb.plot %>% select(oldid = gencove_id, pop))
  
  allpca$lineage2 = allpca$lineage
  allpca$lineage2[allpca$lineage == "C" & allpca$country == "Italy"] = "Ci"
  allpca$lineage2[allpca$lineage == "C" & allpca$country != "Italy"] = "Cc"
  
  
  #relevel
  allpca$lineage = as.factor(allpca$lineage)
  allpca$lineage = relevel(allpca$lineage, "admixed")
  # #explicit names
  # allpca$expl = as.character(allpca$lineage)
  #   allpca$expl[allpca$expl == "A"] = "Africa"
  #   allpca$expl[allpca$expl == "M"] = "N. Europe"
  #   allpca$expl[allpca$expl == "C"] = "S. Europe"
  #   allpca$expl[allpca$expl == "O"] = "Mid. East"
  #   #relevel
  #   allpca$expl= as.factor(allpca$expl)
  #   allpca$expl = relevel(allpca$expl, "USA")
  
  #allpca2 = allpca %>% filter(! state %in% c("Jamaica", "NM"))
  
  
  ggplot(allpca, 
         aes(x = PC1, y = PC2, color = pop, shape = lineage)) + 
    geom_point(size = 3, alpha = 0.5) + 
    labs(shape = "Lineage", color = "Population",
         x = paste0("PC1 (", PCev[1], "%)"),
         y = paste0("PC2 (", PCev[2], "%)"))+
    scale_color_manual(values = plot.colors)+
    theme_bw()
  
  
  
 
  
  
#####
  
  
  
  
  
   
#old
  
#####
  
#other mitos
  
mitoadmix = admix %>% left_join(ahb.plot %>% select(oldid = vcfid, pop, call)) %>%
  filter(!is.na(call))

mitoadmix %>% group_by(call) %>% summarise(A = mean(A), M = mean(M), O = mean(O),
                                           C = mean(C))


plot(mitoadmix$A, mitoadmix$call == "A1e")
  



###

#unsupervised admix

#connect to most-recent admix run
fam = read.delim("data/admix.fam", sep = "", header = F) %>% select(oldid = V1)
sum(!grepl("SRR", fam$oldid))
admix = read.delim("data/unsupadmix.4.Q", header = F, sep ="")
admix = cbind(admix, fam)

#attach lineage
reflins = read.delim("/home/dylan/Documents/bees/harpurlab/project/popgen/admixResults/fullref/refData.txt")
ahb.bar = admix %>% left_join(ahb4 %>% select(oldid, state))
ahb.bar = ahb.bar %>% left_join(reflins %>% select(oldid = SRR, lineage))
ahb.bar$state[!is.na(ahb.bar$lineage)] = ahb.bar$lineage[!is.na(ahb.bar$lineage)]

ahb.bar = ahb.bar %>% select(-lineage)


#add bar id
ahb.bar = ahb.bar %>% arrange(state, V1)
ahb.bar$barid = (1:nrow(ahb.bar))

admix.melt = ahb.bar %>% melt(id.vars = c("oldid", "barid", "state"), 
                              measure.vars = c("V1", "V2", "V3", "V4"),
                              variable.name = "fam")
admix.decode = admix.melt
admix.decode$fam = as.character(admix.decode$fam)
admix.decode$fam[admix.decode$fam == "V1"] = "A"
admix.decode$fam[admix.decode$fam == "V2"] = "C"
admix.decode$fam[admix.decode$fam == "V3"] = "M"
admix.decode$fam[admix.decode$fam == "V4"] = "O"
admix.decode = admix.decode %>% filter(!grepl("SRR", oldid))


#bar chart
ggplot(data = admix.decode) +
  geom_bar(aes(x = as.factor(barid), y = value, fill = fam), 
           stat='identity', width = 1) +
  facet_grid(cols = vars(state), scales = "free_x") +
  scale_fill_brewer(palette = "PRGn") +
  labs(x = "sample", y = "Proportion of Genome", fill = "Lineage") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank())
  


  
##TODO: calculate heterozygosity as well, that would be interesting...
  


### old 
#-----------
  
  #just curious, PA N vs NA
  pa = ahb.plot %>% filter(state == "PA")
  pa$aggBin = grepl("NA", pa$id)
  ggplot(pa, aes(x = A, y = aggBin)) + geom_point()
  
  
  #just admixed pops
  ggplot(ahb.plot %>% filter(pop %in% c("AZ", "FL")), aes(x = A, y = AmitoBin, color = pop)) + 
    geom_point(size = 3, alpha = 0.5) +
    geom_line(data = newdat %>% filter(pop %in% c("AZ", "FL"))) +
    facet_grid(facet = vars(pop)) + theme_bw() + 
    labs(x = "Genomic A-lineage", y = "A-lineage Mitochondrion", 
         title = "Detecting AHB") +
    scale_y_continuous(breaks = c(0,1), 
                       labels = c("False", "True"), limits = c(-.2, 1.2)) +
    theme(legend.position = "none")
  
  
    #scale_y_continuous(breaks = c(-3, 0,1, 3), labels = c("-","False", "True","-")) 
  
  # ahb.filter = ahb.plot %>% filter(bitscore > 1000)
  # ggplot(ahb.filter, aes(x = A, y = AmitoBin, color = state)) + geom_point() + 
  #   facet_grid(facet = vars(pop)) + theme_bw() + 
  #   labs(x = "A admixture", y = "A mitotype detection", 
  #        title = "Mitotyping by Population")
  
  
  oldahb = read.csv("AHB_meta2.csv", header = T)
  oldahb = oldahb %>% left_join(mito.calls %>% rename(vcfid = s))
  
  oldahb %>% filter(loc == "FL", Amito == 1, call != "A1e")
  oldahb %>% filter(loc == "FL", Amito == 0, call == "A1e")
  
  #read admix
  ahb.full = ahb %>% left_join(fix2 %>% select(oldid, V1, V2, V3, V4))
  in.test = ahb.full %>% filter(state == "IN") %>% 
    select(id, call, V1, V2, V3, V4) %>% arrange(desc(V3))
  
  
  
# very old
  
  
  #read admix
  ahb = ahb %>% left_join(fix2 %>% select(oldid, A = V1))

###comparison and graphing
  pcr.lm = glm(Amito ~ A, data = ahb, family = binomial("logit"))
  with(pcr.lm, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
  
  ahb$Amt = ifelse(ahb$call == "A1e", 1, 0)
  mt.lm = glm(Amt ~ A, data = ahb, family = binomial("logit"))
  with(mt.lm, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
  
  ggplot(ahb %>% filter(Amt != Amito), aes(x= A, y = Amt, color = loc)) + 
    geom_point()
  
  
  ahb.compare = rbind(ahb %>% select(vcfid, loc, A, detection = Amito) %>% 
                        mutate(method = "pcr", detection = as.logical(detection)),
                      ahb %>% select(vcfid, loc, A, detection = Amt) %>% 
                        mutate(method = "blast", detection = as.logical(detection)))
  ahb.compare = ahb.compare %>% filter(!is.na(detection))
  
  ggplot(ahb.compare, aes(x = A, y = detection, color = method)) + 
    geom_jitter(height = 0.03, width = 0, alpha = 0.8, size = 2) + 
    labs(y = "A detection", x = "A admixture") +
    theme_bw()
  
  ggplot(ahb.compare %>% filter(loc == "AZ"), aes(x = A, y = detection, color = method)) + 
    geom_jitter(height = 0.03, width = 0, alpha = 0.8, size = 2) + 
    labs(y = "A detection", x = "A admixture") +
    theme_bw()
  
  
  fail = ahb %>% filter(is.na(call)) %>% select(vcfid)
  write.table(fail, file = "failed_samples_2.txt", row.names = F, quote = F, col.names = F)
  
  

  
  #TODO: why aren't some samples running?
  #TODO: why is are NM mitotypes so unreliable?


### more samples! : AZ, NM, PA, TX, IN, FL, Jamaica
  ahb2 = fix2 %>% filter(state %in% c('AZ', 'NM', 'PA', 'TX', 'IN', "Jamaica", "FL"))
  ahb2$vcfid = ahb2$oldid
  
  #fix oddities
  oddities = ahb$oldid[ahb$vcfid != ahb$oldid]
  for(i in 1:nrow(ahb2)){
    if(ahb2$oldid[i] %in% oddities){
      ahb2$vcfid[i] = ahb2$id[i]
    }
  }
  
  
  #now to find gencove sample names for all these, download fastqs, and call mitotypes!
  allgencove = read.delim("all_gencove_sampleNames.txt", header=F, sep = "\t") %>%
    select(vcfid = V3, gencoveid = V2, qc = V4) %>% filter(!grepl("failed", qc))
  
  
  ahb2 = ahb2 %>% left_join(allgencove)


###unarchiving samples
  library(viscomplexr)
  projid = read.delim("projIDs.txt", header=F, sep = '\t')
  missing = read.delim("missing.txt", header = F)
  ahb = read.csv("AHB_meta3.csv") %>% filter(vcfid %in% missing$V1) %>%
    left_join(projid %>% select(gencoveid = V1, projid = V2))
  
  unarch = ahb %>% group_by(projid) %>% mutate(arg = vector2String(gencoveid)) %>%
    slice(1) %>% select(projid, arg)
  projids = unarch$projid
  unarch$arg = gsub("^c[(]|[])]|[ ]", "", unarch$arg)
  unarch = str_split(unarch$arg, ",")
  names(unarch) = projids
  
  unarch.out = data.frame(projid = NA, cmd = NA)
  for (p in projids){
    x = unarch[[p]]
    y = split(x, ceiling(seq_along(x)/20))
    z = sapply(y, vector2String)
    z = gsub("^c[(]|[])]|[ ]", "", z)
    out = data.frame(projid = p, cmd = z)
    unarch.out = rbind(unarch.out, out)
  }
  unarch.out = unarch.out[-1,] %>% arrange(projid)
  
  
  write.table(unarch.out, file = "unarchive_cmd.txt", 
              col.names = F, quote = F, row.names = F, sep = "\t")
  
  write.table(ahb %>% select(projid, gencoveid), file = "unarchive_cmd.txt", 
              col.names = F, quote = F, row.names = F, sep = "\t")

#we have a bunch from the same colony! they should ALL have the same mitotype -> good test!
  
ahb2.out = ahb2 %>% 
  select(id, id2, oldid, vcfid, A, state, gencoveid, projid) %>% 
  arrange(state, A)
write.csv(ahb2.out, file = "AHB_meta3.csv", row.names = F, quote = F)

ahb2 = read.csv("AHB_meta3.csv")


#####