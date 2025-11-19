#####
# load libraries and setup
#####

  library(purrr)
  library(tidyverse)
  library(reshape2)
  library(readxl)
  library(cowplot)
  library(usmap)
  library(rnaturalearth)
  library(sf)
  library(ggspatial)
  
  select = dplyr::select
  theme_set(theme_bw())
  
  setwd("/home/dylan/Documents/bees/harpurlab/project/popgen/ahb")

#TODO:
  # full PCA figure
  
#####
# concordance figure
#####
  
  
gconc = read.csv("data/Andy bee concordance table - Sheet1.csv")
  colnames(gconc) = tolower(colnames(gconc))
  
  ggplot(gconc, aes(x = concordance)) + 
    geom_histogram(binwidth = 0.05) + 
    labs(y = "samples (count)", x = "concordance rate",
         title = "Concordance: 10x versus 4x imputed")


#####
# choose references
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
    reffam.filter = reffam %>% filter(DOGANTID %in% include$`Sample ID` | harpur)
  # # #just keep AMCO
  # #   reffam.filter = reffam.filter %>% filter(lineage %in% c('A', 'M', 'C', 'O'))
  # #   reffam.filter %>% group_by(lineage) %>% summarise(n= n())
  #   
  #   #write out
  #   write.table(reffam.filter$SRR, file = "references/allAMCO.txt",
  #               quote = F, row.names = F, col.names = F)
  #   
    
  #just keep pure samples (>90% lineage assignment)
    refadmix = read.delim("data/allRefThin.7.Q", header = F, sep ="")
    refadmix$max = apply(refadmix, 1, max)
    refadmix$id = read.delim("data/allRefThin.fam", sep = "", header = F)[,1]
    pure = refadmix$id[refadmix$max > 0.85]
    
    reffam.pure = reffam.filter %>% 
      filter(SRR %in% pure, lineage %in% c('A', 'M', 'C', 'O'))
    
    table(reffam.pure$lineage)
    
  #test
    #real observations: without imputation
    refadmix.unimpt = cbind( read.delim("data/allRefThin.7.Q", header = F, sep =""),
                             read.delim("data/allRefThin.fam", sep = "", header = F) %>% select(oldid = V1)) %>%
      left_join(reffam %>% select(lineage, oldid = SRR))
    
    refadmix.unimpt.melt = refadmix.unimpt %>% melt(id.vars = c("oldid", "lineage"),
                                                    measure.vars = c("V1", "V2", "V3", "V4", "V5", "V6", "V7"),
                                                    variable.name = "fam")
    #bar chart:all
    p.allref = refadmix.unimpt.melt %>% filter(!is.na(lineage)) %>% 
    ggplot( ) +
      geom_bar(aes(x = oldid, y = value, fill = fam), 
               stat='identity', width = 1) +
      facet_grid(cols = vars(lineage), scales = "free_x") +
      scale_fill_brewer(palette = "PRGn") +
      labs(x = "reference genomes", y = "Proportion of Genome", fill = "Lineage",
           title = "All References") + 
      theme_bw() + 
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            legend.position = "none")
    
    
  #select ~20 references of each lineage
  set.seed(2025)
  balanced = reffam.pure %>% filter(F)
  rand1 = sample(c(1:sum(reffam.pure$lineage == "A")), 20, replace = F)
  
  balanced = rbind(reffam.pure %>% filter(lineage == "M"),
                   reffam.pure %>% filter(lineage == "C"),
                   reffam.pure %>% filter(lineage == "O"),
                   (reffam.pure %>% filter(lineage == "A"))[rand1,])
  
  #bar chart:pure
  p.pureref = refadmix.unimpt.melt %>% 
    filter(oldid %in% balanced$SRR) %>% 
    ggplot() +
    geom_bar(aes(x = oldid, y = value, fill = fam), 
             stat='identity', width = 1) +
    facet_grid(cols = vars(lineage), scales = "free_x") +
    scale_fill_brewer(palette = "PRGn") +
    labs(x = "reference genomes", y = "", fill = "Lineage",
         title = "Selected References") + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          legend.position = "none")
  
  #export compound figure
    supp1 = plot_grid(p.allref, p.pureref, labels = "AUTO")
    supp1
    ggsave("manuscript/figs/supp1.pdf", plot = supp1, height = 5, width = 10.5,
           units = "in")
  
  
  # #write out
  #   write.table(balanced$SRR, file = "references/pureRefs.txt",
  #               quote = F, row.names = F, col.names = F)
  #   
  #   
  # #write lists
  # linu = unique(balanced$lineage)
  # for(i in 1:length(linu)){
  # 
  #   p = balanced$SRR[which(balanced$lineage == linu[i])]
  #   write.table(p, file = paste0("references/", linu[i], ".txt"),
  #               quote = F, col.names = F, row.names = F)
  #   # notp = balanced %>% filter(!(SRR %in% p)) %>% select(SRR)
  #   # write.table(notp$SRR, file = paste0("references/not", linu[i], ".txt"),
  #   #             quote = F, col.names = F, row.names = F)
  #   }

#####
# read mitotypes, make calls
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
  
  #test against 3rd party mitotyping 
  ahb = read.csv("ahb_metadata.csv")
  
  ahb = ahb %>% left_join(mito.calls %>% rename(vcfid = s))
  #concordance with FL
  ahb %>% filter(state == "FL", Amito_PCR == 1, call != "A1e")
  ahb %>% filter(state == "FL", Amito_PCR == 0, call == "A1e")
  
  #set up model 
  ahb$AmitoBin = ifelse(ahb$call == "A1e", 1, 0)
  
#####
# Load genotype, PCA, and phenotypic data
#####
  
  #oct25
  #unimpt
  
  loadDat = function(VERSION){
  
    #connect to most-recent admix run
    fam = read.delim(paste0("data/admix.",VERSION,".fam"), sep = "", header = F) %>% 
      select(oldid = V1)
    sum(!grepl("SRR", fam$oldid))
    admix = read.delim(paste0("data/admix.",VERSION,".4.Q"), header = F, sep ="")
    admix = cbind(admix, fam) 
    #local pca
    pca = read.delim(paste0("data/samps.",VERSION, ".eigenvec"), header = F, sep = " ") %>%
      select(1:6)
    colnames(pca) = c("fam", "oldid", "PC1", "PC2", "PC3", "PC4")
    #% explained
    eig = read.delim(paste0("data/samps.",VERSION, ".eigenval"), header = F) %>% filter(V1 > 0)
    eig$x = 1:nrow(eig)
    #ggplot(eig, aes(y = V1, x = x)) + geom_point()
    PCev = round(eig$V1 / sum(eig$V1) * 100, 1)
    #global pca
    allpca = read.delim(paste0("data/all.", VERSION,".eigenvec"), header = F, sep = " ") %>%
      select(1:6)
    colnames(allpca) = c("fam", "oldid", "PC1", "PC2", "PC3", "PC4")
    #% explained
    alleig = read.delim(paste0("data/all.", VERSION,".eigenval"), header = F) %>% filter(V1 > 0)
    allPCev = round(alleig$V1 / sum(alleig$V1) * 100, 1)
    
    #verify lineage identity
    lins = admix %>% left_join(reffam %>% select(oldid = SRR, lineage)) %>% 
      filter(grepl ("SRR", oldid)) %>%
      pivot_longer(starts_with("V")) %>%
      group_by(name, lineage) %>%
      mutate(m = mean(value)) %>% ungroup %>% 
      group_by(name) %>% arrange(desc(m)) %>% slice(1) %>% select(lineage, name)
    
    #rename
    admix = admix %>% rename(A = V4, M = V1, C = V2, O = V3)
    
    ahb.plot =  ahb %>% select(-A) %>%
                left_join(admix %>% select(A, oldid))
  
    #establish populations
    ahb.plot$pop = ahb.plot$state
    ahb.plot$pop[ahb.plot$pop == "Jamaica"] = "CAR"
    ahb.plot$pop = factor(ahb.plot$pop,
                          levels = c("IN", "PA", "FL","TX", "AZ", "CAR"))
    
    #admixture bar chart
    ahb.bar = admix %>% left_join(ahb.plot %>% select(oldid, pop))
    
      #add bar id
      ahb.bar = ahb.bar %>% arrange(pop, A)
      ahb.bar$barid = (1:nrow(ahb.bar))
      #melt --- long format
      admix.melt = ahb.bar %>% melt(id.vars = c("oldid", "barid", "pop"), 
                                    measure.vars = c("A", "M", "C", "O"),
                                    variable.name = "fam")
    
  #create logistic model for mitotypes
    ahb.moddat = ahb.plot %>% filter(!is.na(AmitoBin))
    
    ahbmod = glm(AmitoBin ~ pop * A, family = binomial, data = ahb.moddat)
    
      #fit values to plot regression line
          newdat = map_dfr(unique(ahb.plot$pop), function(x){
                    fakeA = seq(min(ahb.plot$A[ahb.plot$pop == x]), 
                          max(ahb.plot$A[ahb.plot$pop == x]), len = 100)
            return(data.frame(pop = x, A = fakeA))
          })
          newdat$AmitoBin = predict(ahbmod, newdat, type = 'response')
          
    #load phenotypic data
    pheno = read_excel("/home/dylan/Documents/bees/harpurlab/project/popgen/samples/beekData/AZ Colony Sample ID.xlsx",
                       skip = 1) %>% 
      select(id = 2, def = Aggression) %>% 
      mutate(id = toupper(id))
    
    #fix pheno IDs to match genetic IDs
    pheno$id[pheno$id == "E4 OG"] = "E4OG"
    pheno$id[pheno$id == "D6"] = "D6-OG"
    pheno$id[pheno$id == "G10"] = "z-G10"
    pheno$id[pheno$id == "A21B"] = "A21b"
    
    pheno = pheno %>% 
      left_join(ahb.plot %>% select(id, A, AmitoBin)) %>% 
      mutate(AmitoShape = as.character(AmitoBin))
    
    pca = pca %>% 
      left_join(ahb.plot %>% select(oldid, id, id2, state, call, pop)) %>%
      left_join(admix %>% select(oldid, A, M, C, O), by = 'oldid')
    
    pca$pop = factor(pca$pop,
                     levels = c("IN", "PA", "FL","TX", "AZ", "CAR"))
    
    #attache lineage names to PCA
    reflins = read.delim("/home/dylan/Documents/bees/harpurlab/project/popgen/admixResults/fullref/refData.txt")
    allpca = allpca %>% left_join(reflins %>% select(oldid = SRR, lineage, country))
    allpca$lineage[is.na(allpca$lineage)] = "admixed"
    allpca = allpca %>% left_join(ahb.plot %>% select(oldid, pop))
    
    allpca$lineage2 = allpca$lineage
    allpca$lineage2[allpca$lineage == "C" & allpca$country == "Italy"] = "Ci"
    allpca$lineage2[allpca$lineage == "C" & allpca$country != "Italy"] = "Cc"
    
    
    #relevel PCA
    allpca$lineage = as.factor(allpca$lineage)
    allpca$lineage = relevel(allpca$lineage, "admixed")
    

    out = dplyr::lst(admix, pca, PCev, allpca, allPCev, ahb.plot, ahb.bar,
                     newdat, pheno, ahbmod, admix.melt, ahb.moddat)
  }
  
  
  AHBDAT = loadDat("oct25")
  list2env(AHBDAT, .GlobalEnv)
          
#####
# Analysis
#####
          
    #mean admixture statistics
      meanse = function(x){
        m = mean(x, na.rm = T) %>% round(4)
        se = sd(x, na.rm = T) /  sqrt(sum(!is.na(x)))
        se = round(se, 4)
        return(paste0(m,", ",se))
      }
    admix %>% filter(!grepl("SRR", oldid)) %>%
      select(-oldid) %>%
      summarise_all(meanse)
    
    #admixture anova N vs S.
    ahb.ano = ahb.plot
    ahb.ano$globe = ifelse(grepl("IN|PA", ahb.ano$pop), "N", "S")
    globelm = lm(A ~ globe, data = ahb.ano)
    anova(globelm)
    #distributions
    # ggplot(ahb.ano, aes(x = globe, y = A)) + geom_boxplot()
    # ggplot(ahb.ano, aes(x = A)) + geom_histogram() + facet_wrap(facets = vars(globe))
         
    #logistic mitotype model fit
    with(ahbmod, pchisq(null.deviance - deviance, 
                        df.null - df.residual, lower.tail = FALSE))
    
  #summaries
    #n by state
    ahb.plot %>% group_by(pop) %>% summarise(n=n())
    #n by pop (filtered)
    ahb.plot %>% group_by(pop) %>% summarise(n=n())
    #n each mitotype
    ahb.plot %>% group_by(call) %>% summarise(n=n())
    #range of A values
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
    
  
  #logistic regression of mitotypes
    p.logimito = ggplot(ahb.plot, aes(x = A, y = AmitoBin, color = pop)) + 
      geom_point(size = 2, alpha = 0.5) +
      geom_line(data = newdat, color = 'black', linetype = 2) +
      facet_grid(rows = vars(pop)) + theme_bw() + 
      #plot thresholds
      #geom_vline(data = thresh, aes(xintercept = A), linetype = 2) +
      #plot means
      #geom_vline(data = ahb.m, aes(xintercept = m), linetype = 3)+
      #colors
      scale_color_manual(values = plot.colors)+
      #labels etc
      labs(x = "A-lineage Admixture", y = "A-lineage Mitotype", 
           color = "State") +
      scale_y_continuous(breaks = c(0,1), 
                         labels = c("False", "True"), limits = c(-.2, 1.2))+
      theme(legend.position = "none",
            axis.title = element_text(size = 9),
            strip.text.y = element_text(size = 7),
            axis.text = element_text(size = 7))
    
    p.logimito
    
    #admixture bar chart
    p.sampleadmix = ggplot(data = admix.melt %>% filter(!is.na(pop))) +
      geom_bar(aes(x = as.factor(barid), y = value, fill = fam), 
               stat='identity', width = 1) +
      facet_grid(cols = vars(pop), scales = "free_x") +
      scale_fill_brewer(palette = "PRGn") +
      labs(x = "sample", y = "Proportion of Genome", fill = "Lineage") + 
      theme_bw() + 
      theme(axis.text.y = element_text(size = 7),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            legend.title=element_text(size=7),
            legend.text=element_text(size = 6),
            axis.title = element_text(size = 9),
            legend.key.size = unit(0.8,"line"),
            strip.text.x = element_text(size = 7),
            legend.margin = margin(l = -8)) + 
      guides(fill = guide_legend(override.aes = list(linewidth = 0.8)))
      
    p.sampleadmix
    
    #phenotypic data
    
    range(pheno$def)

    #plot
    qdmse = summary(lm(def ~ A, data = pheno))$residuals **2 %>%
      mean() %>% 
      round(.,3)
    qdmse = paste0("MSE = ", qdmse)
    
    p.quantdef = ggplot(pheno, aes(x = A, y = def)) + 
      geom_jitter(color = "#66a61e", size = 2, alpha = 0.5,
                  height = 0.05, width = 0) + 
      geom_smooth(method = 'lm', linetype = 2, color = 'black',
                  linewidth = 0.5) + 
      annotate(geom = "text", x = 0.2, y = 2, label = qdmse, size = 6) +
      labs(x = "A-lineage Admixture",
           y = "Defensive Score") +
      theme(axis.title = element_text(size = 9),
            axis.text = element_text(size = 7)) 
    p.quantdef
    
    
    bdmse = summary(lm(def ~ AmitoBin, data = pheno))$residuals **2 %>%
      mean() %>% 
      round(.,3)
    bdmse = paste0("MSE = ", bdmse)
    
    p.bindef = pheno %>% 
      filter(!is.na(AmitoBin)) %>% 
      #mutate(AmitoBin = as.factor(AmitoBin)) %>% 
      ggplot(aes(x = AmitoBin, y = def)) + 
      geom_jitter(color = "#66a61e", size = 2, alpha = 0.5,
                  height = 0.05, width = 0.05) + 
      geom_smooth(method = 'lm', linetype = 2, color = 'black',
                  linewidth = 0.5) +
     # geom_hline(aes(yintercept = 2), linetype = 2, color = 'red') +
      annotate(geom = "text", x = 0.4, y = 2, label = bdmse, size = 6) +
      scale_x_continuous(breaks = c(0,1), labels = c("False","True"),
                         limits = c(-0.2, 1.2)) +
      labs(x = "A-lineage Mitotype",
           y = "") +
      theme(axis.title = element_text(size = 9),
            axis.text = element_text(size = 7))
    p.bindef
  #   
  # #leave-one-out
  #   #quant
  #   pheno.quant = pheno %>% filter(!is.na(A))
  #   LOO = 
  #   errs = rep(NA, nrow(pheno.quant))
  #   for(i in 1:nrow(pheno.quant)){
  #     xi = pheno.quant[i,]
  #     xm1 = pheno.quant[-i,]
  # 
  #     mod = lm(def ~ A, data = xm1)
  # 
  #     errs[i] = (xi$def - predict(mod, xi))**2
  #   }
  #   mean(errs)

 # #plot false and correct assignments
 #      
 #      pheno.gtest = pheno
 #      pheno.gtest$assign = "correct"
 #      pheno.gtest$assign[pheno.gtest$def < 2 & pheno.gtest$A > 0.5] = 'FalsePositive'
 #      pheno.gtest$assign[pheno.gtest$def >= 2 & pheno.gtest$A < 0.5] = 'FalseNegative'
 #      
 #        
 #      ggplot(pheno.gtest, aes(x = A, y = def, color = assign)) + 
 #        geom_jitter(size = 2, alpha = 0.5,
 #                    height = 0, width = 0) + 
 #        geom_smooth(method = 'lm', linetype = 2, color = 'black',
 #                    linewidth = 0.5) +
 #        geom_hline(aes(yintercept = 2), linetype = 2, color = 'red') +
 #        geom_vline(aes(xintercept = 0.5), linetype = 2, color = 'red') +
 #        labs(x = "A-lineage Admixture",
 #             y = "Defensive Score") +
 #        theme(axis.title = element_text(size = 9),
 #              axis.text = element_text(size = 7))
 #      
 #      table(pheno.gtest$assign)
 #      (17+2)/48
 #      
 #      
 #      pheno.mtest = pheno %>% 
 #        filter(!is.na(AmitoBin))
 #      pheno.mtest$assign = "Correct"
 #      pheno.mtest$assign[pheno.mtest$def < 2 & pheno.mtest$AmitoBin == TRUE] = 'FalsePositive'
 #      pheno.mtest$assign[pheno.mtest$def >= 2 & pheno.mtest$AmitoBin == FALSE] = 'FalseNegative'
 #      
 #      ggplot(pheno.mtest,aes( x = AmitoBin, y = def, color = assign)) + 
 #      geom_jitter(size = 2, alpha = 0.5,
 #                  height = 0, width = 0.05) + 
 #      geom_smooth(method = 'lm', linetype = 2, color = 'black',
 #                  linewidth = 0.5) +
 #      geom_hline(aes(yintercept = 2), linetype = 2, color = 'red') +
 #      scale_x_continuous(breaks = c(0,1), labels = c("False","True"),
 #                         limits = c(-0.2, 1.2)) +
 #      labs(x = "A-lineage Mitotype",
 #           y = "") +
 #      theme(axis.title = element_text(size = 9),
 #            axis.text = element_text(size = 7))
 #      
 #      table(pheno.mtest$assign)
 #    
 #    #TODO: add sig to fig
 #    
 #    summary(lm(def ~ A, data = pheno))
 #    summary(lm(def ~ AmitoBin, data = pheno))
    
    #PCA


  #PCA by population
    ggplot(pca, aes(x = PC1, y = PC2, color = pop)) + geom_point(alpha = 0.5, size = 3) +
      theme_bw() +
      scale_color_manual(values = plot.colors) +
      labs(color = "Population",
           x = paste0("PC1 (", PCev[1], "%)"),
           y = paste0("PC2 (", PCev[2], "%)"))
    
  #PCA by lineage
    pcapair = pca %>% select(PC1, PC2, PC3, PC4, A, M, C, O)
    pairs(pcapair, upper.panel = NULL)
    
    
  #all pca
    #load data
    #ahb4 = read.csv("AHBmeta4.csv") %>% rename(oldid = gencove_id)
    
    

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
    
    
    p.globpca = ggplot(allpca, 
           aes(x = PC1, y = PC2, color = pop, shape = lineage)) + 
      geom_point(size = 2, alpha = 0.5) + 
      labs(shape = "Lineage", color = "Population",
           x = paste0("PC1 (", PCev[1], "%)"),
           y = paste0("PC2 (", PCev[2], "%)"))+
      scale_color_manual(values = plot.colors)+
      theme_bw() +
      theme(legend.title=element_text(size=7),
            legend.text=element_text(size = 6),
            legend.spacing.y = unit(1, 'line'),
            axis.title = element_text(size = 9),
            axis.text = element_text(size = 7),
            legend.key.size = unit(0.6,"line"),
            legend.margin = margin(l = -8)) +
      guides(color = guide_legend(byrow = TRUE),
             shape = guide_legend(byrow = TRUE))
    p.globpca
    
    # ggplot(allpca, 
    #        aes(x = PC3, y = PC4, color = pop, shape = lineage)) + 
    #   geom_point(size = 3, alpha = 0.5) + 
    #   labs(shape = "Lineage", color = "Population",
    #        x = paste0("PC3 (", PCev[3], "%)"),
    #        y = paste0("PC4 (", PCev[4], "%)"))+
    #   scale_color_manual(values = plot.colors)+
    #   theme_bw()
  
    
  ### PC's vs Admixture 
    pca.admix = allpca %>% left_join(admix, by = 'oldid') %>%
      filter(!grepl("SRR", oldid))
    
    summary(lm(A~PC2, data = pca.admix))
    summary(lm(C~PC2, data = pca.admix))
    summary(lm(M~PC1, data = pca.admix))
    summary(lm(O~PC1, data = pca.admix))
    
  
#####
# Compound Figures
#####
    #main: ancestry results
    sub.a = plot_grid(p.sampleadmix, p.globpca, labels = c("A", "B"), nrow = 1,
                      rel_widths = c(0.9,1))
    sub.b = plot_grid(p.logimito, NULL, labels = c("C", ""), nrow = 1, 
                      rel_widths = c(1, 0.2))
    sub.c = plot_grid(p.quantdef, p.bindef, labels = c("D", "E"), nrow = 1)
    
    p.ancfig = plot_grid(sub.a, sub.b, sub.c, nrow = 3, 
                         rel_heights = c(0.9,1,0.7))
    
    ggsave("manuscript/figs/fig3.pdf", plot = p.ancfig, height = 9, width = 7,
           units = "in")
    
    #extended PCA
    
    
  # ### Compound figure: two smaller figs
  #   ggsave("manuscript/figs/fig3a.pdf", plot = sub.a, height = 3, width = 7,
  #          units = "in")
  #   
  #   sub.c2 = plot_grid(p.quantdef, p.bindef, labels = c("B", "C"), nrow = 1)
  #   
  #   p.ancfig.b = plot_grid(p.logimito, sub.c2, nrow = 2, 
  #                        rel_heights = c(1,0.6), labels = c("A",NULL))
  #   
  #   ggsave("manuscript/figs/fig3b.pdf", plot = p.ancfig.b, height = 8, width = 7,
  #          units = "in")
  #  
  # 
##### 
#compare imputed and unimptued data
##### 
    
    #load unimputed data
    AHBDAT.unimpt = loadDat("unimpt")
    list2env(AHBDAT.unimpt, .GlobalEnv)
    
  
  #PCA
    p.u.globpca = ggplot(allpca, 
                       aes(x = PC1, y = PC2, color = pop, shape = lineage)) + 
      geom_point(size = 2, alpha = 0.5) + 
      labs(shape = "Lineage", color = "Population",
           x = paste0("PC1 (", PCev[1], "%)"),
           y = paste0("PC2 (", PCev[2], "%)"))+
      scale_color_manual(values = plot.colors)+
      theme_bw() +
      theme(legend.title=element_text(size=7),
            legend.text=element_text(size = 6),
            legend.spacing.y = unit(1, 'line'),
            axis.title = element_text(size = 9),
            axis.text = element_text(size = 7),
            legend.key.size = unit(0.6,"line"),
            legend.margin = margin(l = -8),
            legend.position = "none") +
      guides(color = guide_legend(byrow = TRUE),
             shape = guide_legend(byrow = TRUE))
    p.u.globpca
    
  #admixture
    p.u.sampleadmix = ggplot(data = admix.melt %>% filter(!is.na(pop))) +
      geom_bar(aes(x = as.factor(barid), y = value, fill = fam), 
               stat='identity', width = 1) +
      facet_grid(cols = vars(pop), scales = "free_x") +
      scale_fill_brewer(palette = "PRGn") +
      labs(x = "sample", y = NULL, fill = "Lineage") + 
      theme_bw() + 
      theme(axis.text.y = element_text(size = 7),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            legend.title=element_text(size=7),
            legend.text=element_text(size = 6),
            axis.title = element_text(size = 9),
            legend.key.size = unit(0.8,"line"),
            strip.text.x = element_text(size = 7),
            legend.margin = margin(l = -8),
            legend.position = "none") + 
      guides(fill = guide_legend(override.aes = list(linewidth = 0.8)))
    
    p.u.sampleadmix
  #mitotypes
    #logistic regression of mitotypes
    p.u.logimito = ggplot(ahb.plot, aes(x = A, y = AmitoBin, color = pop)) + 
      geom_point(size = 2, alpha = 0.5) +
      geom_line(data = newdat, color = 'black', linetype = 2) +
      facet_grid(rows = vars(pop)) + theme_bw() + 
      scale_color_manual(values = plot.colors)+
      #labels etc
      labs(x = "A-lineage Admixture", y = NULL, 
           color = "State") +
      scale_y_continuous(breaks = c(0,1), 
                         labels = c("False", "True"), limits = c(-.2, 1.2))+
      theme(legend.position = "none",
            axis.title = element_text(size = 9),
            strip.text.y = element_text(size = 7),
            axis.text = element_text(size = 7),
            axis.text.y = element_blank()) 
      
    p.u.logimito
    
      #correct calls
      sum( (ahbmod$fitted.values > 0.5 & ahb.moddat$AmitoBin == 1) |
             (ahbmod$fitted.values < 0.5 & ahb.moddat$AmitoBin == 0)) / nrow(ahb.moddat)
      #false positives
      sum( (ahbmod$fitted.values > 0.5 & ahb.moddat$AmitoBin == 0)  )
      #false negatives
      sum( (ahbmod$fitted.values < 0.5 & ahb.moddat$AmitoBin == 1)  )
    
  #phenotypes
      qdmse = summary(lm(def ~ A, data = pheno))$residuals **2 %>%
        mean() %>% 
        round(.,3)
      qdmse = paste0("MSE = ", qdmse)
      
      p.u.quantdef = ggplot(pheno, aes(x = A, y = def)) + 
        geom_jitter(color = "#66a61e", size = 2, alpha = 0.5,
                    height = 0.05, width = 0) + 
        geom_smooth(method = 'lm', linetype = 2, color = 'black',
                    linewidth = 0.5) + 
        annotate(geom = "text", x = 0.4, y = 2, label = qdmse, size = 6) +
        #ylim(0.5, 3) +
        labs(x = "A-lineage Admixture",
             y = "Defensive Score") +
        theme(axis.title = element_text(size = 9),
              axis.text = element_text(size = 7)) 
      p.u.quantdef
      

#####
# Compound Figure
#####   
      
  impt.comp = plot_grid(NULL, NULL, p.sampleadmix, p.u.sampleadmix,
            p.globpca, p.u.globpca,
            p.logimito, p.u.logimito,
            p.quantdef, p.u.quantdef,
            nrow = 5,
            rel_heights = c(0.15, 01,1,1.5,1),
            rel_widths = c(1.1, 1),
            labels = c("Imputed data", "Unimputed data", 
                       NA, NA, NA, NA, NA, NA, NA, NA))
      
  ggsave("manuscript/figs/fig4.pdf", plot = impt.comp, height = 11, width = 7,
         units = "in")
  
  
#####
# test admixture panel on references
#####

  #real observations: without imputation
  refadmix.unimpt = cbind( read.delim("data/reference.unimpt.4.Q", header = F, sep =""),
                           read.delim("data/reference.unimpt.fam", sep = "", header = F) %>% select(oldid = V1)) %>%
    left_join(reffam %>% select(lineage, oldid = SRR))
  
  refadmix.unimpt.melt = refadmix.unimpt %>% melt(id.vars = c("oldid", "lineage"),
                                                  measure.vars = c("V1", "V2", "V3", "V4"),
                                                  variable.name = "fam")
  #bar chart
  plot.unimpt = ggplot(data = refadmix.unimpt.melt) +
    geom_bar(aes(x = oldid, y = value, fill = fam), 
             stat='identity', width = 1) +
    facet_grid(cols = vars(lineage), scales = "free_x") +
    scale_fill_brewer(palette = "PRGn") +
    labs(x = "reference genomes", y = "Proportion of Genome", fill = "Lineage",
         title = "Unimputed Sites n=95") + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          legend.position = "none")
  
  
  #with imputation
  refadmix = cbind( read.delim("data/reference.oct25.4.Q", header = F, sep =""),
                    read.delim("data/reference.oct25.fam", sep = "", header = F) %>% select(oldid = V1)) %>%
    left_join(reffam %>% select(lineage, oldid = SRR))
  
  refadmix.melt = refadmix %>% melt(id.vars = c("oldid", "lineage"),
                                    measure.vars = c("V1", "V2", "V3", "V4"),
                                    variable.name = "fam")
  
  plot.impt = ggplot(data = refadmix.melt) +
    geom_bar(aes(x = oldid, y = value, fill = fam), 
             stat='identity', width = 1) +
    facet_grid(cols = vars(lineage), scales = "free_x") +
    scale_fill_brewer(palette = "PRGn") +
    labs(x = "reference genomes", y = NULL, fill = "Ancestry\nComponents",
         title = "Imputed Sites n=30199") + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank())
  
  #supplementary
  supp.admixcomp = plot_grid(plot.unimpt, plot.impt)
  supp.admixcomp


  
  
##### 
# mapping sampling locations
##### 

#load gps positions of samples
ahb = read.csv("ahb_metadata.csv")
samp.gps = read_excel("../old_ahb/sra/biosample2.xlsx",) %>%
  select(1,4,10)
  samp.gps[ samp.gps == "NA"] = NA
  
  #format 
  samp.ll = str_split(samp.gps$Lat_Lon, ",", simplify = T)
    colnames(samp.ll) = c("lat", "lon")
    
  samp.gps = cbind(samp.gps[,-3], samp.ll) %>%
    filter(!is.na(lat)) %>%
    group_by(lat) %>%
    slice(1) %>%
    mutate(lat = as.numeric(lat),
           lon = as.numeric(lon)) %>%
    rename(id = `Sample Name`) %>%
    left_join(ahb %>% select(id, state))
  
    samp.gps$state[samp.gps$state == "Jamaica"] = "CAR"
    samp.gps$state = factor(samp.gps$state,
                     levels = c("IN", "PA", "FL","TX", "AZ", "CAR"))
    
  #anonymize CAR sample locations
    caranonctr = c(-80, 18)
    caranonpts = cbind(rnorm(5, caranonctr[1], 0.5), rnorm(5, caranonctr[2], 0.5))
    samp.gps[samp.gps$state == "CAR", c("lon", "lat")] = caranonpts
    
    
  
  #build map limits using max and min values from data (plus a bit of padding)
  map.limits = rbind(range(samp.gps$lon), range(samp.gps$lat))
    map.limits[,1] = map.limits[,1] - 1.2
    map.limits[,2] = map.limits[,2] + 1.2
    #a bit more padding south
    map.limits[2,1] = map.limits[2,1] - 2
  
  
  #region <- ne_states()
  region = ne_download( scale = 10L, type = "states", category = "cultural")
  world = ne_countries(scale='medium')
  
  
  #convert to sf object
  converted <- st_as_sf( samp.gps %>% select(lon, lat)
                     , coords = c("lon", "lat"), crs = st_crs(region), 
                     agr = "constant")
  samp.gps = cbind(samp.gps, converted)
  
  #color pallette 
  plot.colors = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")
  
  
  
  #full map (world)
  cutout = ggplot(data = world) +
      geom_sf() +
      #add rectangle for cutout
      geom_rect(xmin = map.limits[1,1], xmax = map.limits[1,2], 
                ymin = map.limits[2,1], ymax = map.limits[2,2], 
                fill = NA, colour = "red", size = 0.8) +
      theme(panel.background = element_rect(fill = "azure"),
            panel.border = element_rect(fill = NA))
  
  #just sampling range
  sampling = ggplot(data = region) +
      geom_sf() +
      geom_sf(data = samp.gps, aes(fill = state, geometry = geometry),
              size = 2, shape = 23) +
      # annotate(geom = "text", x = -90, y = 26, label = "Gulf of Mexico", 
      #          fontface = "italic", color = "grey22", size = 6) +
      coord_sf(xlim = map.limits[1,], 
               ylim = map.limits[2,], expand = FALSE) +
    scale_fill_manual(values = plot.colors) +
    labs(fill = "Population") +
    #scale bars
    annotation_scale(
      location = "tl",
      bar_cols = c("grey60", "white")) +
    #north arrow
    annotation_north_arrow(
      location = "tl", which_north = "true",
      pad_x = unit(0.28, "in"), pad_y = unit(0.45, "in"),
      style = ggspatial::north_arrow_nautical(
        fill = c("grey40", "white"),
        line_col = "grey20")) +
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank(), panel.background = element_rect(fill = "azure"), 
            panel.border = element_rect(fill = NA)) +
    #annotation for CAR samples
    annotate("point",x = -80 , y = 18, size = 80, alpha = 0.2,
             color = '#e6ab02')
  #sampling
  
  #plot both together
  full.map = ggdraw(sampling) + 
    draw_plot(cutout, width = 0.28, height = 0.32, 
              x = 0.1, y = 0.1)
  full.map
  #save
  
  pdf("manuscript/figs/fig2.pdf", width = 7, height = 7)

    full.map

  dev.off()
 
  
  
# Mulu Missing
# found = read.delim("data/muluFound2.txt", header = F)
# meta = read.delim("../old_ahb/sra/sra_meta.tsv", sep = "\t")
# 
# mulusent = read.csv("../old_ahb/sra/muluexport.csv")
# 
#  sum(mulusent[,1] %in% found[,1])
#  sum(found[,1] %in% mulusent[,1])

# nrow(mulusent)
# nrow(found)/2


projids = read.delim("sample_lists/projIDs.txt") %>% 
  rename(projid = 2)

found = found %>% filter(grepl("R1", V1)) %>% 
  rename(filename = V1) %>% 
  left_join(meta %>% select(sample_name, filename)) %>% 
  mutate(gencoveid = gsub("_.*", "", filename))
  left_join()
   
  
  
#Metadata for SRA  
#####
  
# sra = ahb %>% select(`Sample Name` = id, state,
#                       `mitochondrial haplotype` = call)
# sra$state[sra$state == "NM"] = "AZ"
#   #long-form state names
#   long.state =
#     data.frame(state= unique(sra$state),
#                long = c("Indiana", "Pennsylvania", "Florida",
#                         "Texas", "Arizona", "Jamaica"))
# 
#   sra = sra %>% left_join(long.state)
#   sra$`geographic location` = paste0("USA:", sra$long)
# 
# #collection data
#   fl.dates = sra %>% filter(state == "FL") %>% select(`Sample Name`)
# 
#     fl.sheet = read_excel("/home/dylan/Documents/bees/harpurlab/project/popgen/samples/beekData/FL_AHB_combined.xlsx") %>%
#       select(2, 5, 7, 8)
#       colnames(fl.sheet) = c("Sample Name", "date", "lat", "lon")
# 
#       fl.sheet$lat = fl.sheet$lat %>% as.numeric() %>% round(6)
#       fl.sheet$lon = fl.sheet$lon %>% as.numeric() %>% round(6)
# 
#       fl.sheet$`Sample Name` = gsub("[ |.].*", "",fl.sheet$`Sample Name`)
#       fl.sheet$gps = paste0(fl.sheet$lat, ", ", fl.sheet$lon)
#       fl.sheet$gps[fl.sheet$gps == "NA, NA"] = NA
#       fl.sheet = fl.sheet %>% select(-lat, -lon)
# 
#     fl.dates = fl.dates %>% left_join(fl.sheet, multiple = 'first')
#     #use ID for dates that do not match
#     fl.match = fl.dates %>% filter(is.na(date))
#     fl.dates = fl.dates %>% filter(!is.na(date))
# 
#     fl.match$date = gsub("^.*([0-9]{8}).*$", "\\1", fl.match$`Sample Name`)
#     fl.match$date = gsub('^([0-9]{2})([0-9]{2})([0-9]{4})$', '\\3-\\1-\\2',
#                          fl.match$date) %>% as.Date()
#     fl.dates = rbind(fl.dates, fl.match)
#     #other metadata
#     fl.dates$collected_by = "Florida Department of Agriculture and Consumer Services"
# 
#   az.dates = sra %>% filter(state == "AZ") %>% select(`Sample Name`)
#     az.sheet = read_excel("/home/dylan/Documents/bees/harpurlab/project/popgen/samples/beekData/AZ Colony Sample ID.xlsx",
#                           skip = 1) %>%
#       rename(id = 2) %>%
#       select(id, gps = `GPS coordinates`) %>%
#       mutate(id = toupper(id)) %>%
#       mutate(id = gsub(" ", "", id))
#       az.sheet$id[az.sheet$id == "G10"] = "z-G10"
#       az.sheet$id[az.sheet$id == "D6"] = "D6-OG"
# 
#     az.dates$date = as.Date("2021-08-17")
#     az.dates = az.dates %>% left_join(az.sheet %>% select('Sample Name' = id, gps))
#     #other meta data
#     az.dates$collected_by = "Ethel Villalobos"
# 
#   tx.dates = sra %>% filter(state == "TX") %>% select(`Sample Name`)
#     tx.dates$date = paste0(gsub("^.*[-]", "", tx.dates$`Sample Name`),
#                            "-08-01") %>% as.Date()
#     #other meta
#     tx.sheet = read_excel("/home/dylan/Documents/bees/harpurlab/project/popgen/samples/beekData/Sample_Info_Dickey_Rangel.xlsx") %>%
#       select(`Sample Name` = 1, lat = 5, lon = 6) %>%
#       mutate(lat = gsub(" N", "", lat),
#              lon = gsub(" W", "", lon)) %>%
#       mutate(lat = round(as.numeric(lat), 6),
#              lon = round(as.numeric(lon), 6)) %>%
#       mutate(gps = paste0(lat, ", ", lon)) %>%
#       select(-lat, -lon)
#     tx.dates = tx.dates %>% left_join(tx.sheet)
#     tx.dates$collected_by = "Myra Dickey"
# 
#   pa.dates = sra %>% filter(state == "PA") %>% select(`Sample Name`)
#     pa.dates$date = as.Date("2022-08-01")
#     #other meta
#     pa.sheet = read_excel("/home/dylan/Documents/bees/harpurlab/project/popgen/samples/beekData/dkredit_CDean_Gencove_Metadata_Corrected.xlsx",
#                           range = "A1:H79") %>%
#       select(id = 2, lat = 7, lon = 8) %>%
#       filter(!is.na(lat)) %>%
#       mutate(lat = round(as.numeric(lat), 6),
#              lon = round(as.numeric(lon), 6)) %>%
#       mutate(gps = paste0(lat, ", ", lon)) %>%
#       select(-lat, -lon) %>%
#       mutate(site = toupper(gsub("^(.*)[ |-][N|A].*", "\\1", id))) %>%
#       distinct(gps, site)
# 
#     pa.dates = pa.dates %>%
#       mutate(site2 = toupper(`Sample Name`))
#       pa.dates$gps = NA
# 
#       for(i in 1:nrow(pa.sheet)){
#         pa.dates$gps[grepl(pa.sheet$site[i], pa.dates$site2)] = pa.sheet$gps[i]
#       }
# 
#     pa.dates = pa.dates %>% select(-site2)
#     pa.dates$collected_by = "Charles C. Dean"
# 
# 
#   in.dates = sra %>% filter(state == "IN") %>% select(`Sample Name`)
#     in.dates$date = as.Date("2020-06-01")
#     #KF
#     in.dates$date[grepl("IN23KF", in.dates$`Sample Name`)] = as.Date("2023-06-01")
#     #KG
#     in.dates$date[grepl("202.Q", in.dates$`Sample Name`)] =
#       paste0(gsub("^.*[_]([0-9]{4}).*", "\\1",
#                   in.dates$`Sample Name`[grepl("202.Q", in.dates$`Sample Name`)]),
#              "-06-01")
#     #other meta
#     in.dates$gps = NA
# 
#     maddie.sheet = read_excel("/home/dylan/Documents/bees/harpurlab/project/popgen/samples/beekData/carpenter_gpsfix.xlsx") %>%
#       select(site = 1, gps = 4)
#     for(i in 1:nrow(maddie.sheet)){
#       in.dates$gps[grepl(maddie.sheet$site[i], in.dates$`Sample Name`)] = maddie.sheet$gps[i]
#     }
# 
#     in.dates$collected_by = "Madeline Carpenter"
#     in.dates$collected_by[grepl("IN23KF", in.dates$`Sample Name`)] = "Ken Foster"
#     in.dates$gps[grepl("IN23KF", in.dates$`Sample Name`)] = "40.46469, -86.7221"
# 
#     in.dates$collected_by[grepl("202.Q", in.dates$`Sample Name`)] = "Krispn Given"
#     in.dates$gps[grepl("202.Q", in.dates$`Sample Name`)] = "40.42857, -86.94861"
# 
#   all.dates = rbind(fl.dates, az.dates, tx.dates, pa.dates, in.dates)
#   sra = sra %>% left_join(all.dates)
# 
# 
# 
# #format for sra
#   sra.out = sra %>%
#     mutate(Organism = "Apis mellifera",
#            tissue = "full body",
#            `isolation source` = "Managed colony",
#            Sex = "female",
#            Dev_stage = "adult") %>%
#     select(`Sample Name`, Organism, collection_date = date, `geographic location`,
#            tissue, `isolation source`, collected_by, Sex, Dev_stage, Lat_Lon = gps)
# 
# 
#   #remove Jamaica (already in SRA)
#   #sra.out = sra.out %>% filter(!grepl("Jamaica", `geographic location`))
# 
# 
#   #fix date format
#   sra.out$collection_date = as.character(sra.out$collection_date)
#   sra.out$sample_identifier = sra.out$`Sample Name`
#   sra.out$ecotype = "Admixed"
# 
#     write_tsv(sra.out, file = "../old_ahb/sra/biosample.tsv")
#     
#     
# #sra metadata
#     ahb3 = read.csv("AHB_meta3.csv")
#     nrow(ahb3)
#     nrow(ahb4)
#     ahb4 = read.csv("AHBmeta4.csv")
#     
#     meta = ahb4 %>% left_join(ahb3 %>% select(oldid, gencoveid)) %>%
#       select(id, gencoveid) %>%
#       mutate(sample_name = id,
#              library_ID = id,
#              title = "full genome sequence of apis mellifera from North America",
#              library_strategy = "WGS",
#              library_source = "GENOMIC",
#              library_selection = "other",
#              library_layout = "paired",
#              platform = "ILLUMINA",
#              instrument_model = "Illumina NovaSeq 6000",
#              design_description = "Bead-Linked Transposome using Illumina Nextera DNA Flex Library Preparation kit",
#              filetype = "fastq",
#              filename = paste0(gencoveid, "_R1.fastq.gz"),
#              filename2 = paste0(gencoveid, "_R2.fastq.gz")) %>%
#       select(-id, -gencoveid)
#     
#     #remove jamaica (already submitted)
#     meta = meta %>% filter(!grepl("C2021", sample_name))
    
    #write_tsv(meta, file = "../old_ahb/sra/sra_meta.tsv")

#working with gencove samples
#####  
    

# #unarchiving samples
#   library(viscomplexr)
#   #projid = read.delim("projIDs.txt", header=F, sep = '\t')
#   missing = read.delim("outputs/failed.out", header = F)
#   ahb3 = read.csv("AHB_meta3.csv") %>% filter(vcfid %in% missing$V1)
#   
#   unarch = ahb3 %>% group_by(projid) %>% mutate(arg = vector2String(gencoveid)) %>%
#     slice(1) %>% select(projid, arg)
#   projids = unarch$projid
#   unarch$arg = gsub("^c[(]|[])]|[ ]", "", unarch$arg)
#   unarch = str_split(unarch$arg, ",")
#   names(unarch) = projids
#   
#   unarch.out = data.frame(projid = NA, cmd = NA)
#   for (p in projids){
#     x = unarch[[p]]
#     y = split(x, ceiling(seq_along(x)/20))
#     z = sapply(y, vector2String)
#     z = gsub("^c[(]|[])]|[ ]", "", z)
#     out = data.frame(projid = p, cmd = z)
#     unarch.out = rbind(unarch.out, out)
#   }
#   unarch.out = unarch.out[-1,] %>% arrange(projid)
#     
#   # write.table(unarch.out, file = "unarchive_cmd2.txt", 
#   #             col.names = F, quote = F, row.names = F, sep = "\t")
#   
#   
#   # write.table(ahb %>% select(projid, gencoveid), file = "unarchive_cmd.txt", 
#   #             col.names = F, quote = F, row.names = F, sep = "\t")
# 
#       
# # downloading fastq
#   dl.list = ahb3$gencoveid
#   write.table(dl.list, "outputs/dl-list.txt", quote = F, row.names = F, col.names = F)
#   