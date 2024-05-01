library(purrr)
library(tidyverse)
library(reshape2)
library(readxl)
select = dplyr::select

setwd("/home/dylan/Documents/bees/harpurlab/project/popgen/ahb")


##TODO:
  #fix multiple ahb csv's
  #change AZ/NM all to AZ


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
  
  ahb.plot =  ahb4 %>% select(-A) %>%
              left_join(admix %>% select(A, oldid))
  
  
  #establish populations
  ahb.plot$pop = ahb.plot$state
  ahb.plot$pop[grepl("NM|AZ", ahb.plot$pop)] = "AZ"
  ahb.plot$pop = factor(ahb.plot$pop,
                        levels = c("IN", "PA", "FL", "Jamaica","TX", "AZ"))
  
#create model
  
  ahb.moddat = ahb.plot %>% filter(!is.na(AmitoBin))
  
  ahbmod = glm(AmitoBin ~ pop + A, family = binomial, data = ahb.moddat)
  summary(ahbmod)$aic
  
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
  
  ahb.bar = admix %>% left_join(ahb.plot %>% select(oldid, pop))
  
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
    left_join(ahb.plot %>% select(oldid, id, id2, state, call, pop)) %>%
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
  allpca = allpca %>% left_join(ahb.plot %>% select(oldid, pop))
  
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

  
#Metadata for SRA  
#####
  
  #TODO: 
  # lat/lon
  # isolation source (wild or managed)
  
sra = ahb4 %>% select(`Sample Name` = id, state, 
                      `mitochondrial haplotype` = call)
sra$state[sra$state == "NM"] = "AZ"
  #long-form state names
  long.state = 
    data.frame(state= unique(sra$state),
               long = c("Indiana", "Pennsylvania", "Florida", 
                        "Texas", "Arizona", "Jamaica"))
  
  sra = sra %>% left_join(long.state)
  sra$`geographic location` = paste0("USA:", sra$long)
    
#collection data
  fl.dates = sra %>% filter(state == "FL") %>% select(`Sample Name`)

    fl.sheet = read_excel("/home/dylan/Documents/bees/harpurlab/project/popgen/samples/beekData/FL_AHB_combined.xlsx") %>%
      select(2, 5, 7, 8)
      colnames(fl.sheet) = c("Sample Name", "date", "lat", "lon")
      
      fl.sheet$lat = fl.sheet$lat %>% as.numeric() %>% round(6)
      fl.sheet$lon = fl.sheet$lon %>% as.numeric() %>% round(6)

      fl.sheet$`Sample Name` = gsub("[ |.].*", "",fl.sheet$`Sample Name`)
      fl.sheet$gps = paste0(fl.sheet$lat, ", ", fl.sheet$lon)
      fl.sheet$gps[fl.sheet$gps == "NA, NA"] = NA
      fl.sheet = fl.sheet %>% select(-lat, -lon)
      
    fl.dates = fl.dates %>% left_join(fl.sheet, multiple = 'first')
    #use ID for dates that do not match
    fl.match = fl.dates %>% filter(is.na(date))
    fl.dates = fl.dates %>% filter(!is.na(date))
    
    fl.match$date = gsub("^.*([0-9]{8}).*$", "\\1", fl.match$`Sample Name`)
    fl.match$date = gsub('^([0-9]{2})([0-9]{2})([0-9]{4})$', '\\3-\\1-\\2', 
                         fl.match$date) %>% as.Date()
    fl.dates = rbind(fl.dates, fl.match)
    #other metadata
    fl.dates$collected_by = "Florida Department of Agriculture and Consumer Services"
    
  az.dates = sra %>% filter(state == "AZ") %>% select(`Sample Name`)
    az.sheet = read_excel("/home/dylan/Documents/bees/harpurlab/project/popgen/samples/beekData/AZ Colony Sample ID.xlsx",
                          skip = 1) %>%
      rename(id = 2) %>%
      select(id, gps = `GPS coordinates`) %>%
      mutate(id = toupper(id)) %>%
      mutate(id = gsub(" ", "", id))
      az.sheet$id[az.sheet$id == "G10"] = "z-G10"
      az.sheet$id[az.sheet$id == "D6"] = "D6-OG"
    
    az.dates$date = as.Date("2021-08-17")
    az.dates = az.dates %>% left_join(az.sheet %>% select('Sample Name' = id, gps))
    #other meta data
    az.dates$collected_by = "Ethel Villalobos"
    
  tx.dates = sra %>% filter(state == "TX") %>% select(`Sample Name`)
    tx.dates$date = paste0(gsub("^.*[-]", "", tx.dates$`Sample Name`),
                           "-08-01") %>% as.Date()
    #other meta
    tx.sheet = read_excel("/home/dylan/Documents/bees/harpurlab/project/popgen/samples/beekData/Sample_Info_Dickey_Rangel.xlsx") %>%
      select(`Sample Name` = 1, lat = 5, lon = 6) %>%
      mutate(lat = gsub(" N", "", lat),
             lon = gsub(" W", "", lon)) %>%
      mutate(lat = round(as.numeric(lat), 6),
             lon = round(as.numeric(lon), 6)) %>%
      mutate(gps = paste0(lat, ", ", lon)) %>%
      select(-lat, -lon)
    tx.dates = tx.dates %>% left_join(tx.sheet)
    tx.dates$collected_by = "Myra Dickey"
  
  pa.dates = sra %>% filter(state == "PA") %>% select(`Sample Name`)
    pa.dates$date = as.Date("2022-08-01")
    #other meta
    pa.sheet = read_excel("/home/dylan/Documents/bees/harpurlab/project/popgen/samples/beekData/dkredit_CDean_Gencove_Metadata_Corrected.xlsx",
                          range = "A1:H79") %>%
      select(id = 2, lat = 7, lon = 8) %>%
      filter(!is.na(lat)) %>%
      mutate(lat = round(as.numeric(lat), 6),
             lon = round(as.numeric(lon), 6)) %>%
      mutate(gps = paste0(lat, ", ", lon)) %>%
      select(-lat, -lon) %>%
      mutate(site = toupper(gsub("^(.*)[ |-][N|A].*", "\\1", id))) %>%
      distinct(gps, site)
    
    pa.dates = pa.dates %>% 
      mutate(site2 = toupper(`Sample Name`))
      pa.dates$gps = NA
      
      for(i in 1:nrow(pa.sheet)){
        pa.dates$gps[grepl(pa.sheet$site[i], pa.dates$site2)] = pa.sheet$gps[i]
      }
    
    pa.dates = pa.dates %>% select(-site2)
    pa.dates$collected_by = "Charles C. Dean"
    

  in.dates = sra %>% filter(state == "IN") %>% select(`Sample Name`)
    in.dates$date = as.Date("2020-06-01")
    #KF
    in.dates$date[grepl("IN23KF", in.dates$`Sample Name`)] = as.Date("2023-06-01")
    #KG
    in.dates$date[grepl("202.Q", in.dates$`Sample Name`)] =
      paste0(gsub("^.*[_]([0-9]{4}).*", "\\1", 
                  in.dates$`Sample Name`[grepl("202.Q", in.dates$`Sample Name`)]),
             "-06-01")
    #other meta
    in.dates$gps = NA
    
    maddie.sheet = read_excel("/home/dylan/Documents/bees/harpurlab/project/popgen/samples/beekData/carpenter_gpsfix.xlsx") %>%
      select(site = 1, gps = 4)
    for(i in 1:nrow(maddie.sheet)){
      in.dates$gps[grepl(maddie.sheet$site[i], in.dates$`Sample Name`)] = maddie.sheet$gps[i]
    }
    
    in.dates$collected_by = "Madeline Carpenter"
    in.dates$collected_by[grepl("IN23KF", in.dates$`Sample Name`)] = "Ken Foster"
    in.dates$gps[grepl("IN23KF", in.dates$`Sample Name`)] = "40.46469, -86.7221"
    
    in.dates$collected_by[grepl("202.Q", in.dates$`Sample Name`)] = "Krispn Given"
    in.dates$gps[grepl("202.Q", in.dates$`Sample Name`)] = "40.42857, -86.94861"
    
  all.dates = rbind(fl.dates, az.dates, tx.dates, pa.dates, in.dates)
  sra = sra %>% left_join(all.dates)
  
    
    
#format for sra
  sra.out = sra %>%
    mutate(Organism = "Apis mellifera",
           tissue = "full body",
           `isolation source` = "Managed colony",
           Sex = "female",
           Dev_stage = "adult") %>%
    select(`Sample Name`, Organism, collection_date = date, `geographic location`,
           tissue, `isolation source`, collected_by, Sex, Dev_stage, Lat_Lon = gps)
    
  
  #remove Jamaica (already in SRA)
  sra.out = sra.out %>% filter(!grepl("Jamaica", `geographic location`))
    # jam = read_excel("/home/dylan/Documents/bees/harpurlab/project/popgen/jamaica/biosample.xlsx",
    #                  sheet = 'Sheet1') %>%
    #   select(-breed)
    # sra.out = rbind(sra.out, jam)
    
  #fix date format
  sra.out$collection_date = as.character(sra.out$collection_date)
  sra.out$sample_identifier = sra.out$`Sample Name`
  sra.out$ecotype = "Admixed"
    
    write_tsv(sra.out, file = "../old_ahb/sra/biosample.tsv")
    
  
  
  
  