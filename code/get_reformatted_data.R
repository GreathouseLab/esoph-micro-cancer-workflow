# Reformat data

# transform to relative abundances 
# in scale of 0-100
phylo.data.nci.umd <- transform_sample_counts(phylo.data.nci.umd, function(x){100*(x / sum(x))})
phylo.data.tcga.RNAseq <- transform_sample_counts(phylo.data.tcga.RNAseq, function(x){100*(x / sum(x))})
phylo.data.tcga.WGS <- transform_sample_counts(phylo.data.tcga.WGS, function(x){100*(x / sum(x))})

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
dat.16s <- dat.16s %>%
  group_by(OTU) %>%
  mutate(
    BarrettsHist = ifelse(Barretts. == "Y", 1, 0),
    bar.c = BarrettsHist - mean(BarrettsHist),
    female = ifelse(gender=="F", 1, 0),
    female.c = female-mean(female),
    age.c = age - mean(age),
    BMI.n = ifelse(BMI == "U", NA, BMI),
    BMI.n = as.numeric(BMI.n),
    bmi.c = BMI.n - mean(BMI.n, na.rm=T),
    bmi.c = ifelse(is.na(bmi.c), 0, bmi.c)
  )

dat.rna <- dat.rna %>%
  #group_by(OTU) %>%
  mutate(
    BarrettsHist = ifelse(Barrett.s.Esophagus.Reported == "Yes", 1, 0),
    BarrettsHist = ifelse(BarrettsHist == "Not Available", NA, BarrettsHist),
    bar.c = BarrettsHist - mean(BarrettsHist),
    female = ifelse(Gender=="female", 1, 0),
    female.c = female-mean(female),
    year_of_death = ifelse(is.na(year_of_death), 2021, year_of_death),
    age = year_of_death-year_of_birth,
    age.c = age - mean(age, na.rm=T),
    bmi.c = bmi - mean(bmi, na.rm=T),
    bmi.c = ifelse(is.na(bmi.c), 0, bmi.c)
  )

dat.wgs <- dat.wgs %>%
  #group_by(OTU) %>%
  mutate(
    BarrettsHist = ifelse(Barrett.s.Esophagus.Reported == "Yes", 1, 0),
    BarrettsHist = ifelse(BarrettsHist == "Not Available", NA, BarrettsHist),
    bar.c = BarrettsHist - mean(BarrettsHist, na.rm=T),
    female = ifelse(Gender=="female", 1, 0),
    female.c = female-mean(female, na.rm=T),
    year_of_death = ifelse(is.na(year_of_death), 2021, year_of_death),
    age = year_of_death-year_of_birth,
    age.c = age - mean(age, na.rm=T),
    bmi.c = bmi - mean(bmi, na.rm=T),
    bmi.c = ifelse(is.na(bmi.c), 0, bmi.c)
  )

dat.16s$sample_type <- 0
dat.16s$sample_type[dat.16s$tissue=="T" &
                      dat.16s$Histology=="ADC" &
                      dat.16s$Barretts. == "Y"] <- "EAC tissues w/ Barretts History"
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

dat.rna$EACcomp <- 0
dat.rna$EACcomp[(dat.rna$morphology=="8140/3" |
                   dat.rna$morphology=="8480/3") &
                  dat.rna$SampleType_Level2=="Tumor" &
                  dat.rna$Barrett.s.Esophagus.Reported=="Yes"] <- "EAC tissues w/ Barretts History"
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

# tumor stage variable
dat.16s$tumor.stage <- factor(
  dat.16s$pTNM,
  levels = c("0", "1", "I", "IIA", "IIB", "III", "IV"),
  labels = c("0", "1", "I", "II", "II", "III", "IV"), ordered=T
)

dat.rna[dat.rna == "not reported"] <- NA
dat.rna$tumor.stage <- factor(
  dat.rna$tumor_stage,
  levels = c("stage i", "stage ia", "stage ib", "stage iib", "stage iia", "stage iii", "stage iiia", "stage iiib", "stage iiic", "stage iv", "stage iva"),
  labels = c("I", "I", "I", "II", "II", "III", "III", "III", "III", "IV", "IV"), ordered=T
)

dat.wgs[dat.wgs == "not reported"] <- NA
dat.wgs$tumor.stage <- factor(
  dat.wgs$tumor_stage,
  levels = c("stage i", "stage ia", "stage ib", "stage iib", "stage iia", "stage iii", "stage iiia", "stage iiib", "stage iiic", "stage iv", "stage iva"),
  labels = c("I", "I", "I", "II", "II", "III", "III", "III", "III", "IV", "IV"), ordered=T
)


# subset to fuso. nuc. only
# Streptococcus sanguinis 
# Campylobacter concisus
# Prevotella spp.

# Not uniquely identified in 16s
dat.16s.s <- filter(
  dat.16s,
  OTU %in% c(
    "Fusobacterium_nucleatum",
    unique(dat.16s$OTU[dat.16s$OTU %like% "Streptococcus_sanguinis"]),
    unique(dat.16s$OTU[dat.16s$OTU %like% "Campylobacter_"]),
    unique(dat.16s$OTU[dat.16s$OTU %like% "Prevotella"]))
)
dat.rna.s <- filter(
  dat.rna,
  otu2 %in% c(
    "Fusobacterium nucleatum",
    "Streptococcus sanguinis",
    "Campylobacter concisus",
    unique(dat.rna$otu2[dat.rna$otu2 %like% "Prevotella"]))
)
dat.wgs.s <- filter(
  dat.wgs,
  otu2 %in% c(
    "Fusobacterium nucleatum",
    "Streptococcus sanguinis",
    "Campylobacter concisus",
    unique(dat.wgs$otu2[dat.wgs$otu2 %like% "Prevotella"]))
)

# new names -> due to lack of uniqueness
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
    "Streptococcus spp.*",
    "Campylobacter spp.*",
    "Prevotella spp."
  )
)
dat.rna.s$OTU1 <- factor(
  dat.rna.s$otu2,
  levels = c(
    "Fusobacterium nucleatum",
    "Streptococcus sanguinis",
    "Campylobacter concisus",
    "Prevotella melaninogenica",
    "Prevotella denticola",
    "Prevotella intermedia",
    "Prevotella ruminicola"),
  labels = c(
    "Fusobacterium nucleatum",
    "Streptococcus sanguinis",
    "Campylobacter concisus",
    "Prevotella spp.", 
    "Prevotella spp.", 
    "Prevotella spp.", 
    "Prevotella spp.")
)

dat.wgs.s$OTU1 <- factor(
  dat.wgs.s$otu2,
  levels = c(
    "Fusobacterium nucleatum",
    "Streptococcus sanguinis",
    "Campylobacter concisus",
    "Prevotella melaninogenica",
    "Prevotella denticola",
    "Prevotella intermedia",
    "Prevotella ruminicola"),
  labels = c(
    "Fusobacterium nucleatum",
    "Streptococcus sanguinis",
    "Campylobacter concisus",
    "Prevotella spp.", 
    "Prevotella spp.", 
    "Prevotella spp.", 
    "Prevotella spp.")
)

# rename bacteria/OTU variables
dat.16s.s$OTU <- dat.16s.s$OTU1
dat.rna.s$OTU <- dat.rna.s$OTU1
dat.wgs.s$OTU <- dat.wgs.s$OTU1
dat.rna$OTU <- dat.rna$otu2
dat.wgs$OTU <- dat.wgs$otu2

# fix abundance/OTU in rna and WGS data
dat.rna.s <- dat.rna.s %>%
  dplyr::group_by(OTU, ID)%>%
  mutate(Abundance = mean(Abundance))%>%
  distinct(OTU, ID, .keep_all = T)
dat.wgs.s <- dat.wgs.s %>%
  dplyr::group_by(OTU, ID)%>%
  mutate(Abundance = mean(Abundance))%>%
  distinct(OTU, ID, .keep_all = T)