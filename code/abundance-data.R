# load packages
source("code/load_packages.R")

# get data
source("code/get_cleaned_data.R")

theme_set(theme_bw())
# pal = "Sequential"
# scale_colour_discrete <-  function(palname=pal, ...){
#   scale_colour_brewer(palette=palname, ...)
# }
# scale_fill_discrete <-  function(palname=pal, ...){
#   scale_fill_brewer(palette=palname, ...)
# }


# transform to relative abundances
phylo.data.nci.umd <- transform_sample_counts(phylo.data.nci.umd, function(x){x / sum(x)})
phylo.data.tcga.RNAseq <- transform_sample_counts(phylo.data.tcga.RNAseq, function(x){x / sum(x)})
phylo.data.tcga.WGS <- transform_sample_counts(phylo.data.tcga.WGS, function(x){x / sum(x)})

# melt data down for use
dat.16s <- psmelt(phylo.data.nci.umd)
dat.rna <- psmelt(phylo.data.tcga.RNAseq)
dat.wgs <- psmelt(phylo.data.tcga.WGS)

# fix otu formatting
dat.rna$otu2 <- "a"
dat.wgs$otu2 <- "a"
i <- 1
for(i in 1:nrow(dat.rna)){
  dat.rna$otu2[i] <- str_split(dat.rna$OTU[i], ";")[[1]][7]
}
for(i in 1:nrow(dat.wgs)){
  dat.wgs$otu2[i] <- str_split(dat.wgs$OTU[i], ";")[[1]][7]
}

# set up variables
dat.16s$sample_type <- 0
dat.16s$sample_type[dat.16s$tissue=="T" &
                      dat.16s$Histology=="ADC" &
                      dat.16s$Barretts. == "Y"] <- "EAC tissues w/ Barretts History"
dat.16s$sample_type[dat.16s$tissue=="T" &
                      dat.16s$Histology=="ADC" &
                      dat.16s$Barretts. == "N"] <- "EAC tissues w/o Barretts History"
dat.16s$sample_type[dat.16s$tissue=="N" &
                      dat.16s$Histology=="ADC" &
                      dat.16s$Barretts. == "Y"] <- "EAC-adjacent tissue w/ Barretts History"
dat.16s$sample_type[dat.16s$tissue=="BO" &
                      dat.16s$Histology=="Barrets only"&
                      dat.16s$Barretts. == "Y"] <-"Barretts Only"


dat.rna$sample_type <- 0
dat.rna$sample_type[(dat.rna$morphology=="8140/3" |
                       dat.rna$morphology=="8480/3") &
                      dat.rna$SampleType_Level2=="Tumor" &
                      dat.rna$Barrett.s.Esophagus.Reported=="Yes"] <- "EAC tissues w/ Barretts History"
dat.rna$sample_type[(dat.rna$morphology=="8140/3" | 
                       dat.rna$morphology=="8480/3") &
                      dat.rna$SampleType_Level2=="Normal"&
                      dat.rna$Barrett.s.Esophagus.Reported=="Yes"] <- "EAC-adjacent tissue w/ Barretts History"
dat.rna$sample_type[(dat.rna$morphology=="8140/3" |
                       dat.rna$morphology=="8480/3") &
                      dat.rna$SampleType_Level2=="Tumor" &
                      dat.rna$Barrett.s.Esophagus.Reported=="No"] <- "EAC tissues w/o Barretts History"


dat.rna$EACcomp <- 0
dat.rna$EACcomp[(dat.rna$morphology=="8140/3" |
                   dat.rna$morphology=="8480/3") &
                  dat.rna$SampleType_Level2=="Tumor" &
                  dat.rna$Barrett.s.Esophagus.Reported=="Yes"] <- "EAC tissues w/ Barretts History"
dat.rna$EACcomp[(dat.rna$morphology=="8140/3" |
                   dat.rna$morphology=="8480/3") &
                  dat.rna$SampleType_Level2=="Tumor" &
                  dat.rna$Barrett.s.Esophagus.Reported=="No"] <- "EAC tissues w/o Barretts History"

dat.rna$EACcomp[(dat.rna$morphology=="8140/3" |
                   dat.rna$morphology=="8480/3") &
                  dat.rna$SampleType_Level2=="Tumor" &
                  dat.rna$Barrett.s.Esophagus.Reported=="No"] <- "EAC tissues w/o Barretts History"


dat.wgs$sample_type <- 0
dat.wgs$sample_type[(dat.wgs$morphology=="8140/3" |
                       dat.wgs$morphology=="8480/3") &
                      dat.wgs$SampleType_Level2=="Tumor" &
                      dat.wgs$Barrett.s.Esophagus.Reported=="Yes"] <- "EAC tissues w/ Barretts History"
dat.wgs$sample_type[(dat.wgs$morphology=="8140/3" |
                       dat.wgs$morphology=="8480/3") &
                      dat.wgs$SampleType_Level2=="Tumor" &
                      dat.wgs$Barrett.s.Esophagus.Reported=="No"] <- "EAC tissues w/o Barretts History"
dat.wgs$sample_type[(dat.wgs$morphology=="8140/3" | 
                       dat.wgs$morphology=="8480/3") &
                      dat.wgs$SampleType_Level2=="Normal"&
                      dat.wgs$Barrett.s.Esophagus.Reported=="Yes"] <- "EAC-adjacent tissue w/ Barretts History"

dat.wgs$EACcomp <- 0
dat.wgs$EACcomp[(dat.wgs$morphology=="8140/3" |
                   dat.wgs$morphology=="8480/3") &
                  dat.wgs$SampleType_Level2=="Tumor" &
                  dat.wgs$Barrett.s.Esophagus.Reported=="Yes"] <- "EAC tissues w/ Barretts History"
dat.wgs$EACcomp[(dat.wgs$morphology=="8140/3" |
                   dat.wgs$morphology=="8480/3") &
                  dat.wgs$SampleType_Level2=="Tumor" &
                  dat.wgs$Barrett.s.Esophagus.Reported=="No"] <- "EAC tissues w/o Barretts History"
dat.wgs$EACcomp[(dat.wgs$morphology=="8140/3" |
                   dat.wgs$morphology=="8480/3") &
                  dat.wgs$SampleType_Level2=="Tumor" &
                  dat.wgs$Barrett.s.Esophagus.Reported=="No"] <- "EAC tissues w/o Barretts History"


# make tumor vs normal variable
dat.16s$tumor.cat <- factor(dat.16s$tissue, levels=c("BO", "N", "T"), labels = c("Non-Tumor", "Non-Tumor", "Tumor"))
dat.rna$tumor.cat <- factor(dat.rna$SampleType_Level2, levels=c("Normal", "Tumor"), labels = c("Non-Tumor", "Tumor"))
dat.wgs$tumor.cat <- factor(dat.wgs$SampleType_Level2, levels=c("Normal", "Tumor"), labels = c("Non-Tumor", "Tumor"))

# dataset id
dat.16s$source <- "16s"
dat.rna$source <- "rna"
dat.wgs$source <- "wgs"

# plotting ids
dat.16s$X <- paste0(dat.16s$source, "-", dat.16s$tumor.cat)
dat.rna$X <- paste0(dat.rna$source, "-", dat.rna$tumor.cat)
dat.wgs$X <- paste0(dat.wgs$source, "-", dat.wgs$tumor.cat)

# relabel as (0/1) for analysis
dat.16s$tumor <- as.numeric(factor(dat.16s$tissue, levels=c("BO", "N", "T"), labels = c("Non-Tumor", "Non-Tumor", "Tumor"))) - 1
dat.rna$tumor <- as.numeric(factor(dat.rna$SampleType_Level2, levels=c("Normal", "Tumor"), labels = c("Non-Tumor", "Tumor"))) - 1
dat.wgs$tumor <- as.numeric(factor(dat.wgs$SampleType_Level2, levels=c("Normal", "Tumor"), labels = c("Non-Tumor", "Tumor"))) - 1

# presence- absence
dat.16s$pres <- ifelse(dat.16s$Abundance > 0, 1, 0)
dat.16s$pres[is.na(dat.16s$pres)] <- 0
dat.rna$pres <- ifelse(dat.rna$Abundance > 0, 1, 0)
dat.rna$pres[is.na(dat.rna$pres)] <- 0
dat.wgs$pres <- ifelse(dat.wgs$Abundance > 0, 1, 0)
dat.wgs$pres[is.na(dat.wgs$pres)] <- 0


# subset to fuso. nuc. only
# Streptococcus sanguinis 
# Campylobacter concisus
# Prevotella spp.

dat.16s.s <- filter(
  dat.16s,
  OTU %in% c(
    "Fusobacterium_nucleatum",
    unique(dat.16s$OTU[dat.16s$OTU %like% "Streptococcus_"]),
    unique(dat.16s$OTU[dat.16s$OTU %like% "Campylobacter_"]),
    "Prevotella_melaninogenica")
)
dat.rna.s <- filter(
  dat.rna,
  otu2 %in% c(
    "Fusobacterium nucleatum",
    "Streptococcus sanguinis",
    "Streptococcus oralis",
    "Streptococcus mitis", 
    "Streptococcus pneumoniae", 
    "Streptococcus parasanguinis", 
    "Streptococcus salivarius",
    "Campylobacter concisus",
    "Prevotella melaninogenica")
)
dat.wgs.s <- filter(
  dat.wgs,
  otu2 %in% c(
    "Fusobacterium nucleatum",
    "Streptococcus sanguinis",
    "Streptococcus oralis",
    "Streptococcus mitis", 
    "Streptococcus pneumoniae", 
    "Streptococcus parasanguinis", 
    "Streptococcus salivarius",
    "Campylobacter concisus",
    "Prevotella melaninogenica")
)

# new names
dat.16s.s$OTU1 <- factor(
  dat.16s.s$OTU,
  levels = c(
    "Fusobacterium_nucleatum",
    "Streptococcus_dentisani:Streptococcus_infantis:Streptococcus_mitis:Streptococcus_oligofermentans:Streptococcus_oralis:Streptococcus_pneumoniae:Streptococcus_pseudopneumoniae:Streptococcus_sanguinis",
    "Campylobacter_rectus:Campylobacter_showae",
    "Prevotella_melaninogenica"
  ),
  labels = c(
    "Fusobacterium nucleatum",
    "Streptococcus sanguinis",
    "Campylobacter concisus",
    "Prevotella melaninogenica"
  )
)
dat.rna.s$OTU1 <- factor(
  dat.rna.s$otu2,
  levels = c(
    "Fusobacterium nucleatum",
    "Streptococcus sanguinis",
    "Campylobacter concisus",
    "Prevotella melaninogenica"),
  labels = c(
    "Fusobacterium nucleatum",
    "Streptococcus sanguinis",
    "Campylobacter concisus",
    "Prevotella melaninogenica")
)

dat.wgs.s$OTU1 <- factor(
  dat.wgs.s$otu2,
  levels = c(
    "Fusobacterium nucleatum",
    "Streptococcus sanguinis",
    "Campylobacter concisus",
    "Prevotella melaninogenica"),
  labels = c(
    "Fusobacterium nucleatum",
    "Streptococcus sanguinis",
    "Campylobacter concisus",
    "Prevotella melaninogenica")
)

# rename bacteria
dat.16s.s$OTU <- factor(
  dat.16s.s$OTU,
  levels = c(
    "Fusobacterium_nucleatum",
    "Streptococcus_dentisani:Streptococcus_infantis:Streptococcus_mitis:Streptococcus_oligofermentans:Streptococcus_oralis:Streptococcus_pneumoniae:Streptococcus_pseudopneumoniae:Streptococcus_sanguinis",
    "Campylobacter_rectus:Campylobacter_showae",
    "Prevotella_melaninogenica"
  ),
  labels = c(
    "Fusobacterium nucleatum",
    "Streptococcus spp.",
    "Campylobacter concisus",
    "Prevotella melaninogenica"
  )
)
dat.rna.s$OTU <- factor(
  dat.rna.s$otu2,
  levels = c(
    "Fusobacterium nucleatum",
    "Streptococcus sanguinis",
    "Streptococcus oralis",
    "Streptococcus mitis", 
    "Streptococcus pneumoniae", 
    "Streptococcus parasanguinis", 
    "Streptococcus salivarius",
    "Campylobacter concisus",
    "Prevotella melaninogenica"),
  labels = c(
    "Fusobacterium nucleatum",
    "Streptococcus spp.",
    "Streptococcus spp.",
    "Streptococcus spp.", 
    "Streptococcus spp.", 
    "Streptococcus spp.", 
    "Streptococcus spp.",
    "Campylobacter concisus",
    "Prevotella melaninogenica")
)

dat.wgs.s$OTU <- factor(
  dat.wgs.s$otu2,
  levels = c(
    "Fusobacterium nucleatum",
    "Streptococcus sanguinis",
    "Streptococcus oralis",
    "Streptococcus mitis", 
    "Streptococcus pneumoniae", 
    "Streptococcus parasanguinis", 
    "Streptococcus salivarius",
    "Campylobacter concisus",
    "Prevotella melaninogenica"),
  labels = c(
    "Fusobacterium nucleatum",
    "Streptococcus spp.",
    "Streptococcus spp.",
    "Streptococcus spp.", 
    "Streptococcus spp.", 
    "Streptococcus spp.", 
    "Streptococcus spp.",
    "Campylobacter concisus",
    "Prevotella melaninogenica")
)


t16s <- dat.16s.s %>% filter(sample_type != "0") %>%
  dplyr::group_by(OTU, sample_type) %>%
  dplyr::summarise(
    MA_16s = mean(Abundance)*100
  ) %>%
  filter(
    sample_type %in% c("EAC tissues w/ Barretts History",
                       "EAC tissues w/o Barretts History")
  )

trna <- dat.rna.s %>% filter(sample_type != "0") %>%
  dplyr::group_by(OTU, sample_type) %>%
  dplyr::summarise(
    MA_rna = mean(Abundance, na.rm=T)*100
  ) %>%
  filter(
    sample_type %in% c("EAC tissues w/ Barretts History",
                       "EAC tissues w/o Barretts History")
  )

twgs <- dat.wgs.s %>% filter(sample_type != "0") %>%
  dplyr::group_by(OTU, sample_type) %>%
  dplyr::summarise(
    MA_wgs = mean(Abundance, na.rm=T)*100
  ) %>%
  filter(
    sample_type %in% c("EAC tissues w/ Barretts History",
                       "EAC tissues w/o Barretts History")
  )
tout <- left_join(t16s, trna)
tout <- left_join(tout, twgs)
write.table(tout, "output/mean-abundance-table-2021-01-19.csv",
            row.names = F, sep=",")
