# Figure 2 violin plots

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

#knitr::opts_chunk$set(out.width = "225%", out.height = 10)


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
dat.16s$sample_type <- NA
dat.16s$sample_type[dat.16s$tissue=="T" &
                      dat.16s$Histology=="ADC" &
                      dat.16s$Barretts. == "N"] <- "16SrRNA Tumor"
dat.16s$sample_type[dat.16s$tissue=="N" &
                      dat.16s$Histology=="ADC" &
                      dat.16s$Barretts. == "N"] <- "16SrRNA Non-tumor"

dat.rna$sample_type <- NA
dat.rna$sample_type[(dat.rna$morphology=="8140/3" |
                       dat.rna$morphology=="8480/3") &
                      dat.rna$SampleType_Level2=="Tumor" &
                      dat.rna$Barrett.s.Esophagus.Reported=="No"] <- "RNA-seq Tumor"
dat.rna$sample_type[(dat.rna$morphology=="8140/3" | 
                       dat.rna$morphology=="8480/3") &
                      dat.rna$SampleType_Level2=="Normal"&
                      dat.rna$Barrett.s.Esophagus.Reported=="No"] <- "RNA-seq Non-tumor"

dat.wgs$sample_type <- NA
dat.wgs$sample_type[(dat.wgs$morphology=="8140/3" |
                       dat.wgs$morphology=="8480/3") &
                      dat.wgs$SampleType_Level2=="Tumor" &
                      dat.wgs$Barrett.s.Esophagus.Reported=="No"] <- "WGS Tumor"
dat.wgs$sample_type[(dat.wgs$morphology=="8140/3" | 
                       dat.wgs$morphology=="8480/3") &
                      dat.wgs$SampleType_Level2=="Normal"&
                      dat.wgs$Barrett.s.Esophagus.Reported=="No"] <- "WGS Non-tumor"


# make tumor vs normal variable
dat.16s$tumor.cat <- factor(dat.16s$tissue, levels=c("N", "T"), labels = c("Non-Tumor", "Tumor"))
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
dat.16s$tumor <- as.numeric(factor(dat.16s$tissue, levels=c("N", "T"), labels = c("Non-Tumor", "Tumor"))) - 1
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
  labels = c("0", "1", "I", "II", "II", "III", "IV")
)

dat.rna[dat.rna == "not reported"] <- NA
dat.rna$tumor.stage <- factor(
  dat.rna$tumor_stage,
  levels = c("stage i", "stage ia", "stage ib", "stage iib", "stage iia", "stage iii", "stage iiia", "stage iiib", "stage iiic", "stage iv", "stage iva"),
  labels = c("I", "I", "I", "II", "II", "III", "III", "III", "III", "IV", "IV")
)

dat.wgs[dat.wgs == "not reported"] <- NA
dat.wgs$tumor.stage <- factor(
  dat.wgs$tumor_stage,
  levels = c("stage i", "stage ia", "stage ib", "stage iib", "stage iia", "stage iii", "stage iiia", "stage iiib", "stage iiic", "stage iv", "stage iva"),
  labels = c("I", "I", "I", "II", "II", "III", "III", "III", "III", "IV", "IV")
)


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
    "Streptococcus spp. (not uniquely identified)",
    "Campylobacter spp. (not uniquely identified)",
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
dat.rna.s$OTU <- dat.rna.s$otu2
dat.wgs.s$OTU <- dat.wgs.s$otu2
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

#```




# Violin Plot Fuso

#```{r violin-fuso, warning=F, error=F, message=F, fig.dim=c(7,5), out.height="100%", out.width="100%"}

# merge datasets by subsetting to specific variables then merging
analysis.dat <- dat.16s.s %>% 
  dplyr::mutate(ID = as.factor(accession.number)) %>%
  dplyr::select(OTU, sample_type, Abundance, ID, source)

dat <- dat.rna.s %>% 
  dplyr::select(OTU, sample_type, Abundance, ID, source)
analysis.dat <- full_join(analysis.dat, dat)

dat <- dat.wgs.s %>% 
  dplyr::select(OTU, sample_type, Abundance, ID, source)

analysis.dat <- full_join(analysis.dat, dat) %>%
  mutate(pres = ifelse(Abundance > 0, 1, 0)) %>%
  filter(OTU=="Fusobacterium nucleatum")# create a presence/absences variable


tb <- analysis.dat %>%
  filter(is.na(sample_type)==F)%>%
  group_by(sample_type, OTU) %>%
  summarise(
    N=n(),
    p = sum(pres, na.rm=T),
    percent = p/N*100
  )

kable(tb, format="html")%>%
  kable_styling(full_width = T)

analysis.dat <- analysis.dat %>%
  filter(is.na(sample_type)==F)%>%
  mutate(
    Abund = Abundance*100
  )

#root function
root<-function(x){
  x <- ifelse(x < 0, 0, x)
  x**(0.2)
}
#inverse root function
invroot<-function(x){
  x**(5)
}

# colors
cols <- c("16SrRNA Tumor" = "#C0392B",
          "16SrRNA Non-tumor"= "#336699",
          "RNA-seq Tumor"= "#C0392B",
          "RNA-seq Non-tumor"= "#336699",
          "WGS Tumor" = "#C0392B",
          "WGS Non-tumor"= "#336699")
# "#2C3E50", "#2980B9", "#8E44AD", "#C0392B", "#D35400", "#F39C12", "#F1C40F", "#27AE60", "#196f3e"
#   mycolors = c()
#   mycolors$Tissue["T"]   = "#C0392B"
#   mycolors$Tissue["N"]   = "#336699"


p <- ggplot(analysis.dat, aes(sample_type, Abund, color=sample_type)) +
  geom_violin(scale="width")+
  geom_point(alpha=0.75)+
  scale_y_continuous(
    trans=scales::trans_new("root", root, invroot),
    breaks=c(0, 0.001,0.01, 0.1, 1,10,50, 100),
    labels = c(0, 0.001,0.01, 0.1, 1,10,50, 100),
    limits = c(0, 110)
  )+
  scale_color_manual(values=cols)+
  annotate(
    "text", x=c(1:6), y=c(rep(110, 6)),
    label=c(paste0(round(tb[1,5], 0),"% (",tb[1,4],"/",tb[1,3],")"), 
            paste0(round(tb[2,5], 0),"% (",tb[2,4],"/",tb[2,3],")"),
            paste0(round(tb[3,5], 0),"% (",tb[3,4],"/",tb[3,3],")"),
            paste0(round(tb[4,5], 0),"% (",tb[4,4],"/",tb[4,3],")"),
            paste0(round(tb[5,5], 0),"% (",tb[5,4],"/",tb[5,3],")"),
            paste0(round(tb[6,5], 0),"% (",tb[6,4],"/",tb[6,3],")"))
  )+
  labs(x=NULL, y="% Abundance")+
  theme(
    #axis.text.x = element_text(angle=30, hjust=0.95, vjust=0.95),
    legend.position = "none"
  )
p
ggsave("figure2/figure2_violin_fuso.pdf", p, units = "in", width = 8, height = 4)
ggsave("figure2/figure2_violin_fuso.png", p, units = "in", width = 8, height = 4)

# Violin Plot Strepto
#```{r violin-strepto, warning=F, error=F, message=F, fig.dim=c(7,5), out.height="400px", out.width="100%"}

# merge datasets by subsetting to specific variables then merging
analysis.dat <- dat.16s.s %>% 
  dplyr::mutate(ID = as.factor(accession.number)) %>%
  dplyr::select(OTU, sample_type, Abundance, ID, source)

dat <- dat.rna.s %>% 
  dplyr::select(OTU, sample_type, Abundance, ID, source)
analysis.dat <- full_join(analysis.dat, dat)

dat <- dat.wgs.s %>% 
  dplyr::select(OTU, sample_type, Abundance, ID, source)

analysis.dat <- full_join(analysis.dat, dat) %>%
  mutate(pres = ifelse(Abundance > 0, 1, 0)) %>%
  filter(OTU%like%"Streptococcus")# create a presence/absences variable


tb <- analysis.dat %>%
  filter(is.na(sample_type)==F)%>%
  group_by(sample_type, OTU) %>%
  summarise(
    N=n(),
    p = sum(pres, na.rm=T),
    percent = p/N*100
  )

kable(tb, format="html")%>%
  kable_styling(full_width = T)

analysis.dat <- analysis.dat %>%
  filter(is.na(sample_type)==F)%>%
  mutate(
    Abund = Abundance*100
  )

#root function
root<-function(x){
  x <- ifelse(x < 0, 0, x)
  x**(0.2)
}
#inverse root function
invroot<-function(x){
  x**(5)
}

p <- ggplot(analysis.dat, aes(sample_type, Abund, color=sample_type)) +
  geom_violin(scale="width")+
  geom_point(alpha=0.75)+
  scale_y_continuous(
    trans=scales::trans_new("root", root, invroot),
    breaks=c(0, 0.001,0.01, 0.1, 1,10,50, 100),
    labels = c(0, 0.001,0.01, 0.1, 1,10,50, 100),
    limits = c(0, 110)
  )+
  scale_color_manual(values=cols)+
  annotate(
    "text", x=c(1:6), y=c(rep(110, 6)),
    label=c(paste0(round(tb[1,5], 0),"% (",tb[1,4],"/",tb[1,3],")"), 
            paste0(round(tb[2,5], 0),"% (",tb[2,4],"/",tb[2,3],")"),
            paste0(round(tb[3,5], 0),"% (",tb[3,4],"/",tb[3,3],")"),
            paste0(round(tb[4,5], 0),"% (",tb[4,4],"/",tb[4,3],")"),
            paste0(round(tb[5,5], 0),"% (",tb[5,4],"/",tb[5,3],")"),
            paste0(round(tb[6,5], 0),"% (",tb[6,4],"/",tb[6,3],")"))
  )+
  labs(x=NULL, y="% Abundance")+
  theme(
    legend.position = "none"
  )
p
ggsave("figure2/figure2_violin_strepto.pdf", p, units = "in", width = 8, height = 4)
ggsave("figure2/figure2_violin_strepto.png", p, units = "in", width = 8, height = 4)
#```

# Violin Plot Campy
#```{r violin-campy, warning=F, error=F, message=F, fig.dim=c(7,5), out.height="400px", out.width="100%"}

# merge datasets by subsetting to specific variables then merging
analysis.dat <- dat.16s.s %>% 
  dplyr::mutate(ID = as.factor(accession.number)) %>%
  dplyr::select(OTU, sample_type, Abundance, ID, source)

dat <- dat.rna.s %>% 
  dplyr::select(OTU, sample_type, Abundance, ID, source)
analysis.dat <- full_join(analysis.dat, dat)

dat <- dat.wgs.s %>% 
  dplyr::select(OTU, sample_type, Abundance, ID, source)

analysis.dat <- full_join(analysis.dat, dat) %>%
  mutate(pres = ifelse(Abundance > 0, 1, 0)) %>%
  filter(OTU %like% "Campylobacter")# create a presence/absences variable


tb <- analysis.dat %>%
  filter(is.na(sample_type)==F)%>%
  group_by(sample_type, OTU) %>%
  summarise(
    N=n(),
    p = sum(pres, na.rm=T),
    percent = p/N*100
  )

kable(tb, format="html")%>%
  kable_styling(full_width = T)

analysis.dat <- analysis.dat %>%
  filter(is.na(sample_type)==F)%>%
  mutate(
    Abund = Abundance*100
  )

#root function
root<-function(x){
  x <- ifelse(x < 0, 0, x)
  x**(0.2)
}
#inverse root function
invroot<-function(x){
  x**(5)
}

p <- ggplot(analysis.dat, aes(sample_type, Abund, color=sample_type)) +
  geom_violin(scale="width")+
  geom_point(alpha=0.75)+
  scale_y_continuous(
    trans=scales::trans_new("root", root, invroot),
    breaks=c(0, 0.001,0.01, 0.1, 1,10,50, 100),
    labels = c(0, 0.001,0.01, 0.1, 1,10,50, 100),
    limits = c(0, 110)
  )+
  scale_color_manual(values=cols)+
  annotate(
    "text", x=c(1:6), y=c(rep(110, 6)),
    label=c(paste0(round(tb[1,5], 0),"% (",tb[1,4],"/",tb[1,3],")"), 
            paste0(round(tb[2,5], 0),"% (",tb[2,4],"/",tb[2,3],")"),
            paste0(round(tb[3,5], 0),"% (",tb[3,4],"/",tb[3,3],")"),
            paste0(round(tb[4,5], 0),"% (",tb[4,4],"/",tb[4,3],")"),
            paste0(round(tb[5,5], 0),"% (",tb[5,4],"/",tb[5,3],")"),
            paste0(round(tb[6,5], 0),"% (",tb[6,4],"/",tb[6,3],")"))
  )+
  labs(x=NULL, y="% Abundance")+
  theme(
    legend.position = "none"
  )
p
ggsave("figure2/figure2_violin_campy.pdf", p, units = "in", width = 8, height = 4)
ggsave("figure2/figure2_violin_campy.png", p, units = "in", width = 8, height = 4)
#```

# Violin Plot Prevo
#```{r violin-prevo, warning=F, error=F, message=F, fig.dim=c(7,5), out.height="400px", out.width="100%"}

# merge datasets by subsetting to specific variables then merging
analysis.dat <- dat.16s.s %>% 
  dplyr::mutate(ID = as.factor(accession.number)) %>%
  dplyr::select(OTU, sample_type, Abundance, ID, source)

dat <- dat.rna.s %>% 
  dplyr::select(OTU, sample_type, Abundance, ID, source)
analysis.dat <- full_join(analysis.dat, dat)

dat <- dat.wgs.s %>% 
  dplyr::select(OTU, sample_type, Abundance, ID, source)

analysis.dat <- full_join(analysis.dat, dat) %>%
  mutate(pres = ifelse(Abundance > 0, 1, 0)) %>%
  filter(OTU %like% "Prevotella melaninogenica" | OTU == "Prevotella spp.")# create a presence/absences variable


tb <- analysis.dat %>%
  filter(is.na(sample_type)==F)%>%
  group_by(sample_type, OTU) %>%
  summarise(
    N=n(),
    p = sum(pres, na.rm=T),
    percent = p/N*100
  )

kable(tb, format="html")%>%
  kable_styling(full_width = T)

analysis.dat <- analysis.dat %>%
  filter(is.na(sample_type)==F)%>%
  mutate(
    Abund = Abundance*100
  )

#root function
root<-function(x){
  x <- ifelse(x < 0, 0, x)
  x**(0.2)
}
#inverse root function
invroot<-function(x){
  x**(5)
}

p <- ggplot(analysis.dat, aes(sample_type, Abund, color=sample_type)) +
  geom_violin(scale="width")+
  geom_point(alpha=0.75)+
  scale_y_continuous(
    trans=scales::trans_new("root", root, invroot),
    breaks=c(0, 0.001,0.01, 0.1, 1,10,50, 100),
    labels = c(0, 0.001,0.01, 0.1, 1,10,50, 100),
    limits = c(0, 110)
  )+
  scale_color_manual(values=cols)+
  annotate(
    "text", x=c(1:6), y=c(rep(110, 6)),
    label=c(paste0(round(tb[1,5], 0),"% (",tb[1,4],"/",tb[1,3],")"), 
            paste0(round(tb[2,5], 0),"% (",tb[2,4],"/",tb[2,3],")"),
            paste0(round(tb[3,5], 0),"% (",tb[3,4],"/",tb[3,3],")"),
            paste0(round(tb[4,5], 0),"% (",tb[4,4],"/",tb[4,3],")"),
            paste0(round(tb[5,5], 0),"% (",tb[5,4],"/",tb[5,3],")"),
            paste0(round(tb[6,5], 0),"% (",tb[6,4],"/",tb[6,3],")"))
  )+
  labs(x=NULL, y="% Abundance")+
  theme(
    legend.position = "none"
  )
p
ggsave("figure2/figure2_violin_prevo.pdf", p, units = "in", width = 8, height = 4)
ggsave("figure2/figure2_violin_prevo.png", p, units = "in", width = 8, height = 4)