# Updated heatmaps


# load packages
source("code/load_packages.R")

# Navy, gold, blue
#myColor <- c("#002454", "#FAC01A", "#0186e9")
myColor <- rev(colorRampPalette(c("#416DD4", "#565E6A", "#FFBC33"))(100))
# get data
source("code/get_cleaned_data.R")

theme_set(theme_bw())

# transform to relative abundances
phylo.data.nci.umd <- transform_sample_counts(phylo.data.nci.umd, function(x){x / sum(x)})
phylo.data.tcga.RNAseq <- transform_sample_counts(phylo.data.tcga.RNAseq, function(x){x / sum(x)})
phylo.data.tcga.WGS <- transform_sample_counts(phylo.data.tcga.WGS, function(x){x / sum(x)})

# end here and redo from here down for ESCAh plot. 
# May be excessive but I want to make sure it works before cutting out any code.

# master controls of saving plots
save.plots <- T
save.Date <- Sys.Date()

## =================================================== ##
# NCI 16S Data

# melt data down for use
dat.16s <- psmelt(phylo.data.nci.umd)

# set up variables
dat.16s$sample_type <- 0
dat.16s$sample_type[dat.16s$tissue=="T" &
                      dat.16s$Histology=="ADC" &
                      dat.16s$Barretts. == "Y"] <- "ESCA tissues w/ Barretts History"
dat.16s$sample_type[dat.16s$tissue=="T" &
                      dat.16s$Histology=="ADC" &
                      dat.16s$Barretts. == "N"] <- "ESCA tissues w/o Barretts History"
dat.16s$sample_type[dat.16s$tissue=="N" &
                      dat.16s$Histology=="ADC" &
                      dat.16s$Barretts. == "Y"] <- "ESCA-adjacent tissue w/ Barretts History"
dat.16s$sample_type[dat.16s$tissue=="BO" &
                      dat.16s$Histology=="Barrets only"&
                      dat.16s$Barretts. == "Y"] <-NA

# make tumor vs normal variable
dat.16s$tumor.cat <- factor(dat.16s$tissue, levels=c("N", "T"), labels = c("Non-Tumor", "Tumor"))

# dataset id
dat.16s$source <- "16s"

# plotting ids
dat.16s$X <- paste0(dat.16s$source, "-", dat.16s$tumor.cat)

# relabel as (0/1) for analysis
dat.16s$tumor <- as.numeric(factor(dat.16s$tissue, levels=c("N", "T"), labels = c( "Non-Tumor", "Tumor"))) - 1

# presence- absence
dat.16s$pres <- ifelse(dat.16s$Abundance > 0, 1, 0)
dat.16s$pres[is.na(dat.16s$pres)] <- 0

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


#get counts

# Barretts, blood normal, non-tumor, tumor, tumor w/ barretts

nci_sample_data <- dat.16s.s %>%
  filter(OTU == "Fusobacterium nucleatum")

table(nci_sample_data$tumor)

table(nci_sample_data$tissue)

table(nci_sample_data$Barretts.)

table(nci_sample_data$tumor, nci_sample_data$Barretts.)

nrow(nci_sample_data)
#```


### Relative Abudance Cutoff: 0.01

#```{r ht-slide5-001-plot-w-barretts, warning=F, error=F, message=F}

# melt data down for use
dat.16s <- psmelt(phylo.data.nci.umd)

# set up variables
dat.16s$sample_type <- 0
dat.16s$sample_type[dat.16s$tissue=="T" &
                      dat.16s$Histology=="ADC" &
                      dat.16s$Barretts. == "N"] <- "ESCA tissues"
dat.16s$sample_type[dat.16s$tissue=="N" &
                      dat.16s$Histology=="ADC" &
                      dat.16s$Barretts. == "N"] <- "ESCA-adjacent tissue"

# make tumor vs normal variable
dat.16s$tumor.cat <- factor(dat.16s$tissue, levels=c("N", "T"), labels = c("Non-Tumor", "Tumor"))

# dataset id
dat.16s$source <- "16s"

# plotting ids
dat.16s$X <- paste0(dat.16s$source, "-", dat.16s$tumor.cat)

# relabel as (0/1) for analysis
dat.16s$tumor <- as.numeric(factor(dat.16s$tissue, levels=c("N", "T"), labels = c("Non-Tumor", "Tumor"))) - 1

# presence- absence
dat.16s$pres <- ifelse(dat.16s$Abundance > 0, 1, 0)
dat.16s$pres[is.na(dat.16s$pres)] <- 0

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

analysis.dat <- dat.16s # insert dataset to be used in analysis
avgRelAbundCutoff <- 0.01 # minimum average relative abundance for OTUs

otu.dat <- analysis.dat %>% filter(sample_type != "0") %>%
  dplyr::group_by(OTU) %>%
  dplyr::summarise(AverageRelativeAbundance=mean(Abundance))%>%
  dplyr::filter(AverageRelativeAbundance>=avgRelAbundCutoff) %>%
  dplyr::arrange(desc(AverageRelativeAbundance))

kable(otu.dat[,c(2,1)], format="html", digits=3) %>%
  kable_styling(full_width = T)%>%
  scroll_box(width="100%", height="400px")

plot.dat <- analysis.dat %>% 
  filter(sample_type != "0") %>%
  dplyr::group_by(OTU) %>%
  dplyr::mutate(aveAbund=mean(Abundance)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(aveAbund>=avgRelAbundCutoff)  %>%
  dplyr::mutate(ID = as.factor(accession.number),
                Genus = substr(Genus, 4, 1000),
                Phylum = substr(Phylum, 4, 1000)) %>%
  dplyr::select(sample_type, Phylum, Genus, ID, Abundance, aveAbund, Sample.ID) %>%
  dplyr::filter(Genus != "unassigned")
length(unique(plot.dat$Sample.ID))

# save data
write.csv(plot.dat, file="figure1_heatmap_nci_raw_otu_data.csv")

# widen plot.dat for dendro
dat.wide <- plot.dat %>%
  dplyr::mutate(
    ID = paste0(ID, "_",sample_type)
  ) %>%
  dplyr::select(ID, Genus, Abundance) %>%
  dplyr::group_by(ID, Genus) %>%
  dplyr::summarise(
    Abundance = mean(Abundance)
  ) %>%
  tidyr::pivot_wider(
    id_cols = Genus,
    names_from = ID,
    values_from = Abundance,
    values_fill = 0
  )
rn <- dat.wide$Genus
mat <- as.matrix(dat.wide[,-1])
rownames(mat) <- rn

sample_names <- colnames(mat)

# Obtain the dendrogram
dend <- as.dendrogram(hclust(dist(mat)))
dend_data <- dendro_data(dend)

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data <- with(
  segment(dend_data), 
  data.frame(x = y, y = x, xend = yend, yend = xend))
# Use the dendrogram label data to position the gene labels
gene_pos_table <- with(
  dend_data$labels, 
  data.frame(y_center = x, gene = as.character(label), height = 1))

# Table to position the samples
sample_pos_table <- data.frame(sample = sample_names) %>%
  dplyr::mutate(x_center = (1:n()), 
                width = 1)

# Neglecting the gap parameters
heatmap_data <- mat %>% 
  reshape2::melt(value.name = "expr", varnames = c("gene", "sample")) %>%
  left_join(gene_pos_table) %>%
  left_join(sample_pos_table)

# extract and rejoin sample IDs and sample_type names for plotting
# first for the heatmap data.frame
A <- str_split(heatmap_data$sample, "_")
heatmap_data$ID <- heatmap_data$sample_type <- "0"
for(i in 1:nrow(heatmap_data)){
  heatmap_data$ID[i] <- A[[i]][1]
  heatmap_data$sample_type[i] <- A[[i]][2]
}
# second for the sample position dataframe (dendo)
A <- str_split(sample_pos_table$sample, "_")
sample_pos_table$ID <- sample_pos_table$sample_type <- "0"
for(i in 1:nrow(sample_pos_table)){
  sample_pos_table$ID[i] <- A[[i]][1]
  sample_pos_table$sample_type[i] <- A[[i]][2]
}

# Limits for the vertical axes
gene_axis_limits <- with(
  gene_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
) + 
  0.1 * c(-1, 1) # extra spacing: 0.1

## Build Heatmap Pieces
htmap.size <- numeric(3)
# Part 1: "ESCA-adjacent tissue
hmd <- filter(heatmap_data, sample_type == "ESCA-adjacent tissue")
hmd$x_center <- as.numeric(as.factor(hmd$x_center))
spd <- filter(sample_pos_table, sample_type == "ESCA-adjacent tissue")
spd$x_center <- as.numeric(as.factor(spd$x_center))
htmap.size[3] <- length(spd$x_center)

plt_hmap1 <- ggplot(hmd, 
                    aes(x = x_center, y = y_center, fill = expr, 
                        height = height, width = width)) + 
  geom_tile() +
  #facet_wrap(.~sample_type)+
  scale_fill_gradient2("Rel. Abund.",trans="sqrt", high=myColor[1], mid = myColor[50], low=myColor[100], midpoint = 0.40, breaks=c(0, 0.1, 0.10, 0.30, 0.50, 0.80)) +
  scale_x_continuous(breaks = spd$x_center, 
                     labels = spd$ID, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  scale_y_continuous(breaks = gene_pos_table[, "y_center"], 
                     labels = rep("", nrow(gene_pos_table)),
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "ESCA-adj.", y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank(),#element_text(size = rel(1), hjust = 0.5,vjust=0.5, angle = 90), 
        axis.ticks.y = element_blank(),
        # margin: top, right, bottom, and left
        plot.margin = unit(c(1, 0.01, 0.01, -0.7), "cm"), 
        panel.grid.minor = element_blank(),
        legend.position = "none")

# Part 3: "ESCA tissues"
hmd <- filter(heatmap_data, sample_type == "ESCA tissues")
hmd$x_center <- as.numeric(as.factor(hmd$x_center))
spd <- filter(sample_pos_table, sample_type == "ESCA tissues")
spd$x_center <- as.numeric(as.factor(spd$x_center))
htmap.size[2] <- length(spd$x_center)

plt_hmap2 <- ggplot(hmd, 
                    aes(x = x_center, y = y_center, fill = expr, 
                        height = height, width = width)) + 
  geom_tile() +
  #facet_wrap(.~sample_type)+
  scale_fill_gradient2("Rel. Abund.",trans="sqrt", high=myColor[1], mid = myColor[50], low=myColor[100], midpoint = 0.40, breaks=c(0, 0.10, 0.30, 0.50, 0.80)) +
  scale_x_continuous(breaks = spd$x_center, 
                     labels = spd$ID, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  scale_y_continuous(breaks = gene_pos_table[, "y_center"], 
                     labels = rep("", nrow(gene_pos_table)),
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "ESCA", y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank(),#element_text(size = rel(1), hjust = 0.75, vjust=0.5, angle = 90), 
        axis.ticks.y = element_blank(),
        # margin: top, right, bottom, and left
        plot.margin = unit(c(1, 0.01, 0.01, -0.7), "cm"), 
        panel.grid.minor = element_blank())


# Dendrogram plot
plt_dendr <- ggplot(segment_data) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  scale_x_reverse(expand = c(0, 0.5)) + 
  scale_y_continuous(breaks = gene_pos_table$y_center, 
                     labels = gene_pos_table$gene, 
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "", y = "", colour = "", size = "") +
  theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(1, 0.01, 0.01, -0.7), "cm"))

prntRelAbund <- avgRelAbundCutoff*100


#```{r ht-slide5-001-plot-w-barretts-make, warning=F, error=F, message=F, fig.height=10, out.width="225%", out.height="750px"}
htmap.size <- htmap.size/max(htmap.size)
htmap.size[1] <- 0.25
p <- plt_dendr+plt_hmap2+plt_hmap1+ 
  plot_layout(
    nrow=1, widths = htmap.size,
    guides="collect"
  ) +
  plot_annotation(
    title="NCI-16s Data showing average relative abundance of genera by individual",
    subtitle=paste0("Subset to OTU average relative abundance > ",prntRelAbund,"%")
  )
p
if(save.plots == T){
  ggsave(paste0("figure1_heatmap_nci_",save.Date,".pdf"), plot=p, units="in", width=7, height=6)
  ggsave(paste0("figure1_heatmap_nci_",save.Date,".png"), plot=p, units="in", width=7, height=6)
}


## 
# TCGH Data

#```{r ht-6-data, include=FALSE, echo=FALSE, warning=F, error=F, message=F}

# melt data down for use
dat.rna <- psmelt(phylo.data.tcga.RNAseq)

# fix otu formatting
dat.rna$otu2 <- "a"
i <- 1
for(i in 1:nrow(dat.rna)){
  dat.rna$otu2[i] <- str_split(dat.rna$OTU[i], ";")[[1]][7]
}

# set up variables
dat.rna$sample_type <- 0
dat.rna$sample_type[(dat.rna$morphology=="8140/3" | 
                       dat.rna$morphology=="8480/3") &
                      dat.rna$SampleType_Level2=="Normal"&
                      dat.rna$Barrett.s.Esophagus.Reported=="No"] <- "ESCA-adj."
dat.rna$sample_type[(dat.rna$morphology=="8140/3" |
                       dat.rna$morphology=="8480/3") &
                      dat.rna$SampleType_Level2=="Tumor" &
                      dat.rna$Barrett.s.Esophagus.Reported=="No"] <- "ESCA"


dat.rna$ESCAcomp <- 0
dat.rna$ESCAcomp[(dat.rna$morphology=="8140/3" |
                    dat.rna$morphology=="8480/3") &
                   dat.rna$SampleType_Level2=="Tumor" &
                   dat.rna$Barrett.s.Esophagus.Reported=="No"] <- "ESCA"

dat.rna$ESCAcomp[(dat.rna$morphology=="8140/3" |
                    dat.rna$morphology=="8480/3") &
                   dat.rna$SampleType_Level2=="Normal" &
                   dat.rna$Barrett.s.Esophagus.Reported=="No"] <- "ESCA-adj."


# make tumor vs normal variable
dat.rna$tumor.cat <- factor(dat.rna$SampleType_Level2, levels=c("Normal", "Tumor"), labels = c("Non-Tumor", "Tumor"))

# dataset id
dat.rna$source <- "rna"

# plotting ids
dat.rna$X <- paste0(dat.rna$source, "-", dat.rna$tumor.cat)

# relabel as (0/1) for analysis
dat.rna$tumor <- as.numeric(factor(dat.rna$SampleType_Level2, levels=c("Normal", "Tumor"), labels = c("Non-Tumor", "Tumor"))) - 1

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
# new names
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

# rename bacteria
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


#```{r ht-slide6-001-plot, warning=F, error=F, message=F}

analysis.dat <- dat.rna # insert dataset to be used in analysis
avgRelAbundCutoff <- 0.001 # minimum average relative abundance for OTUs


otu.dat <- analysis.dat %>% filter(sample_type != "0") %>%
  dplyr::group_by(otu2) %>%
  dplyr::summarise(AverageRelativeAbundance=mean(Abundance, na.rm=T))%>%
  dplyr::filter(AverageRelativeAbundance>=avgRelAbundCutoff) %>%
  dplyr::arrange(desc(AverageRelativeAbundance))

kable(otu.dat[,c(2,1)], format="html", digits=3) %>%
  kable_styling(full_width = T)%>%
  scroll_box(width="100%", height="100%")

plot.dat <- analysis.dat %>% filter(sample_type != "0") %>%
  dplyr::group_by(OTU) %>%
  dplyr::mutate(aveAbund=mean(Abundance, na.rm=T)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(aveAbund>=avgRelAbundCutoff) %>%
  dplyr::mutate(ID = as.factor(Patient_ID),
                Abundance = ifelse(is.na(Abundance), 0, Abundance)) %>%
  dplyr::select(sample_type, Phylum, Genus, ID, Abundance, aveAbund)


# widen plot.dat for dendro
dat.wide <- plot.dat %>%
  dplyr::mutate(
    ID = paste0(ID, "_",sample_type)
  ) %>%
  dplyr::select(ID, Genus, Abundance) %>%
  dplyr::group_by(ID, Genus) %>%
  dplyr::summarise(
    Abundance = mean(Abundance)
  ) %>%
  tidyr::pivot_wider(
    id_cols = Genus,
    names_from = ID,
    values_from = Abundance,
    values_fill = 0
  )
rn <- dat.wide$Genus
mat <- as.matrix(dat.wide[,-1])
rownames(mat) <- rn

sample_names <- colnames(mat)

# Obtain the dendrogram
dend <- as.dendrogram(hclust(dist(mat)))
dend_data <- dendro_data(dend)

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data <- with(
  segment(dend_data), 
  data.frame(x = y, y = x, xend = yend, yend = xend))
# Use the dendrogram label data to position the gene labels
gene_pos_table <- with(
  dend_data$labels, 
  data.frame(y_center = x, gene = as.character(label), height = 1))

# Table to position the samples
sample_pos_table <- data.frame(sample = sample_names) %>%
  dplyr::mutate(x_center = (1:n()), 
                width = 1)

# Neglecting the gap parameters
heatmap_data <- mat %>% 
  reshape2::melt(value.name = "expr", varnames = c("gene", "sample")) %>%
  left_join(gene_pos_table) %>%
  left_join(sample_pos_table)

# extract and rejoin sample IDs and sample_type names for plotting
# first for the heatmap data.frame
A <- str_split(heatmap_data$sample, "_")
heatmap_data$ID <- heatmap_data$sample_type <- "0"
for(i in 1:nrow(heatmap_data)){
  heatmap_data$ID[i] <- A[[i]][1]
  heatmap_data$sample_type[i] <- A[[i]][2]
}
# second for the sample position dataframe (dendo)
A <- str_split(sample_pos_table$sample, "_")
sample_pos_table$ID <- sample_pos_table$sample_type <- "0"
for(i in 1:nrow(sample_pos_table)){
  sample_pos_table$ID[i] <- A[[i]][1]
  sample_pos_table$sample_type[i] <- A[[i]][2]
}

# Limits for the vertical axes
gene_axis_limits <- with(
  gene_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
) + 
  0.1 * c(-1, 1) # extra spacing: 0.1

## Build Heatmap Pieces
# object to hold relative prop.
htmap.size <- numeric(3)
# ESCA
hmd <- filter(heatmap_data, sample_type == "ESCA")
hmd$x_center <- as.numeric(as.factor(hmd$x_center))
spd <- filter(sample_pos_table, sample_type == "ESCA")
spd$x_center <- as.numeric(as.factor(spd$x_center))
htmap.size[2] <- length(spd$x_center)

plt_hmap1 <- ggplot(hmd, 
                    aes(x = x_center, y = y_center, fill = expr, 
                        height = height, width = width)) + 
  geom_tile() +
  #facet_wrap(.~sample_type)+
  scale_fill_gradient2("Rel. Abund.",trans="sqrt", high=myColor[1], mid = myColor[50], low=myColor[100], midpoint = 0.10, breaks=c(0, 0.05, 0.10, 0.15, 0.20)) +
  scale_x_continuous(breaks = spd$x_center, 
                     labels = spd$ID, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  scale_y_continuous(breaks = gene_pos_table[, "y_center"], 
                     labels = rep("", nrow(gene_pos_table)),
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "ESCA", y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank(),#element_text(size = rel(1), hjust = 0.5,vjust=0.5, angle = 90), 
        # margin: top, right, bottom, and left
        #axis.ticks.y = element_blank(),
        #axis.title.x = element_text(angle=90, hjust=0.75, vjust=0.5),
        plot.margin = unit(c(1, 0.01, 0.01, -0.7), "cm"), 
        panel.grid = element_blank(),
        legend.position = "none")

# Part 2: "ESCA-adjacent"
hmd <- filter(heatmap_data, sample_type == "ESCA-adj.")
hmd$x_center <- as.numeric(as.factor(hmd$x_center))
spd <- filter(sample_pos_table, sample_type == "ESCA-adj.")
spd$x_center <- as.numeric(as.factor(spd$x_center))
htmap.size[3] <- length(spd$x_center)

plt_hmap2 <- ggplot(hmd, 
                    aes(x = x_center, y = y_center, fill = expr, 
                        height = height, width = width)) + 
  geom_tile() +
  #facet_wrap(.~sample_type)+
  scale_fill_gradient2("Rel. Abund.",trans="sqrt", high=myColor[1], mid = myColor[50], low=myColor[100], midpoint = 0.10, breaks=c(0, 0.05, 0.10, 0.15, 0.20)) +
  scale_x_continuous(breaks = spd$x_center, 
                     labels = spd$ID, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  scale_y_continuous(breaks = gene_pos_table[, "y_center"], 
                     labels = rep("", nrow(gene_pos_table)),
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "ESCA-adj.", y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank(),#element_text(size = rel(1), hjust = 0.5,vjust=0.5, angle = 90), 
        axis.ticks.y = element_blank(),
        #axis.title.x = element_text(angle=90, hjust=0.75, vjust=0.5),
        # margin: top, right, bottom, and left
        plot.margin = unit(c(1, 0.01, 0.01, -0.7), "cm"), 
        panel.grid.minor = element_blank())


# Dendrogram plot
plt_dendr <- ggplot(segment_data) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  scale_x_reverse(expand = c(0, 0.5)) + 
  scale_y_continuous(breaks = gene_pos_table$y_center, 
                     labels = gene_pos_table$gene, 
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "", y = "", colour = "", size = "") +
  theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(1, 0.01, 0.01, -0.7), "cm"))

prntRelAbund <- avgRelAbundCutoff*100


#```{r ht-slide6-0001-plot-make, warning=F, error=F, message=F, fig.height=5,out.width="100%", out.height="500px"}
htmap.size <- htmap.size/max(htmap.size)
htmap.size[1] <- 0.25
p <- plt_dendr+plt_hmap1+plt_hmap2+ 
  plot_layout(
    nrow=1, widths = htmap.size,
    guides="collect"
  ) +
  plot_annotation(
    title="TCGA RNAseq Data showing average relative abundance of genera by individual",
    subtitle=paste0("Subset to OTU average relative abundance > ",prntRelAbund,"%")
  )
p

if(save.plots == T){
  ggsave(paste0("figure1_heatmap_tcga_rnaseq_",save.Date,".pdf"), plot=p, units="in", width=7, height=5)
  ggsave(paste0("figure1_heatmap_tcga_rnaseq_",save.Date,".png"), plot=p, units="in", width=7, height=5)
}


## ==== 
# TCGA WGS
# melt data down for use
dat.wgs <- psmelt(phylo.data.tcga.WGS)

# fix otu formatting
dat.wgs$otu2 <- "a"
for(i in 1:nrow(dat.wgs)){
  dat.wgs$otu2[i] <- str_split(dat.wgs$OTU[i], ";")[[1]][7]
}

# set up variables
dat.wgs$sample_type <- 0
dat.wgs$sample_type[(dat.wgs$morphology=="8140/3" |
                       dat.wgs$morphology=="8480/3") &
                      dat.wgs$SampleType_Level2=="Tumor" &
                      dat.wgs$Barrett.s.Esophagus.Reported=="No"] <- "ESCA"

dat.wgs$sample_type[(dat.wgs$morphology=="8140/3" | 
                       dat.wgs$morphology=="8480/3") &
                      dat.wgs$SampleType_Level2=="Normal"&
                      dat.wgs$Barrett.s.Esophagus.Reported=="No"] <- "ESCA-adj."

dat.wgs$ESCAcomp <- 0
dat.wgs$ESCAcomp[(dat.wgs$morphology=="8140/3" |
                    dat.wgs$morphology=="8480/3") &
                   dat.wgs$SampleType_Level2=="Tumor" &
                   dat.wgs$Barrett.s.Esophagus.Reported=="No"] <- "ESCA"
dat.wgs$ESCAcomp[(dat.wgs$morphology=="8140/3" |
                    dat.wgs$morphology=="8480/3") &
                   dat.wgs$SampleType_Level2=="Normal" &
                   dat.wgs$Barrett.s.Esophagus.Reported=="No"] <- "ESCA-adj."


# make tumor vs normal variable
dat.wgs$tumor.cat <- factor(dat.wgs$SampleType_Level2, levels=c("Normal", "Tumor"), labels = c("Non-Tumor", "Tumor"))

# dataset id
dat.wgs$source <- "wgs"

# plotting ids
dat.wgs$X <- paste0(dat.wgs$source, "-", dat.wgs$tumor.cat)

# relabel as (0/1) for analysis
dat.wgs$tumor <- as.numeric(factor(dat.wgs$SampleType_Level2, levels=c("Normal", "Tumor"), labels = c("Non-Tumor", "Tumor"))) - 1

# subset to fuso. nuc. only
# Streptococcus sanguinis 
# Campylobacter concisus
# Prevotella spp.
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

#``{r ht-slide7-001-plot, warning=F, error=F, message=F}

analysis.dat <- dat.wgs # insert dataset to be used in analysis
avgRelAbundCutoff <- 0.01 # minimum average relative abundance for OTUs


otu.dat <- analysis.dat %>% filter(sample_type != "0") %>%
  dplyr::group_by(otu2) %>%
  dplyr::summarise(AverageRelativeAbundance=mean(Abundance, na.rm=T))%>%
  dplyr::filter(AverageRelativeAbundance>=avgRelAbundCutoff) %>%
  dplyr::arrange(desc(AverageRelativeAbundance))

kable(otu.dat[,c(2,1)], format="html", digits=3) %>%
  kable_styling(full_width = T)%>%
  scroll_box(width="100%", height="100%")

plot.dat <- analysis.dat %>% filter(sample_type != "0") %>%
  dplyr::group_by(OTU) %>%
  dplyr::mutate(aveAbund=mean(Abundance, na.rm=T)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(aveAbund>=avgRelAbundCutoff) %>%
  dplyr::mutate(ID = as.factor(Patient_ID),
                Abundance = ifelse(is.na(Abundance), 0, Abundance)) %>%
  dplyr::select(sample_type, Phylum, Genus, ID, Abundance, aveAbund) %>%
  dplyr::filter(Genus != "Unclassified")


# widen plot.dat for dendro
dat.wide <- plot.dat %>%
  dplyr::mutate(
    ID = paste0(ID, "_",sample_type)
  ) %>%
  dplyr::select(ID, Genus, Abundance) %>%
  dplyr::group_by(ID, Genus) %>%
  dplyr::summarise(
    Abundance = mean(Abundance)
  ) %>%
  tidyr::pivot_wider(
    id_cols = Genus,
    names_from = ID,
    values_from = Abundance,
    values_fill = 0
  )
rn <- dat.wide$Genus
mat <- as.matrix(dat.wide[,-1])
rownames(mat) <- rn

sample_names <- colnames(mat)

# Obtain the dendrogram
dend <- as.dendrogram(hclust(dist(mat)))
dend_data <- dendro_data(dend)

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data <- with(
  segment(dend_data), 
  data.frame(x = y, y = x, xend = yend, yend = xend))
# Use the dendrogram label data to position the gene labels
gene_pos_table <- with(
  dend_data$labels, 
  data.frame(y_center = x, gene = as.character(label), height = 1))

# Table to position the samples
sample_pos_table <- data.frame(sample = sample_names) %>%
  dplyr::mutate(x_center = (1:n()), 
                width = 1)

# Neglecting the gap parameters
heatmap_data <- mat %>% 
  reshape2::melt(value.name = "expr", varnames = c("gene", "sample")) %>%
  left_join(gene_pos_table) %>%
  left_join(sample_pos_table)

# extract and rejoin sample IDs and sample_type names for plotting
# first for the heatmap data.frame
A <- str_split(heatmap_data$sample, "_")
heatmap_data$ID <- heatmap_data$sample_type <- "0"
for(i in 1:nrow(heatmap_data)){
  heatmap_data$ID[i] <- A[[i]][1]
  heatmap_data$sample_type[i] <- A[[i]][2]
}
# second for the sample position dataframe (dendo)
A <- str_split(sample_pos_table$sample, "_")
sample_pos_table$ID <- sample_pos_table$sample_type <- "0"
for(i in 1:nrow(sample_pos_table)){
  sample_pos_table$ID[i] <- A[[i]][1]
  sample_pos_table$sample_type[i] <- A[[i]][2]
}

# Limits for the vertical axes
gene_axis_limits <- with(
  gene_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
) + 
  0.1 * c(-1, 1) # extra spacing: 0.1

## Build Heatmap Pieces
htmap.size <- numeric(3)
# ESCA
hmd <- filter(heatmap_data, sample_type == "ESCA")
hmd$x_center <- as.numeric(as.factor(hmd$x_center))
spd <- filter(sample_pos_table, sample_type == "ESCA")
spd$x_center <- as.numeric(as.factor(spd$x_center))
htmap.size[2] <- length(spd$x_center)

plt_hmap1 <- ggplot(hmd, 
                    aes(x = x_center, y = y_center, fill = expr, 
                        height = height, width = width)) + 
  geom_tile() +
  #facet_wrap(.~sample_type)+
  scale_fill_gradient2("Rel. Abund.",trans="sqrt", high=myColor[1], mid = myColor[50], low=myColor[100], midpoint = 0.40, breaks=c(0, 0.10, 0.30, 0.50, 0.80)) +
  scale_x_continuous(breaks = spd$x_center, 
                     labels = spd$ID, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  scale_y_continuous(breaks = gene_pos_table[, "y_center"], 
                     labels = rep("", nrow(gene_pos_table)),
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "ESCA", y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank(), #element_text(size = rel(1), hjust = 0.5,vjust=0.5, angle = 90), 
        # margin: top, right, bottom, and left
        #axis.ticks.y = element_blank(),
        #axis.title.x = element_text(angle=90, hjust=0.75, vjust=0.5),
        plot.margin = unit(c(1, 0.01, 0.01, -0.7), "cm"), 
        panel.grid = element_blank(),
        legend.position = "none")

# Part 2: "ESCA-adjacent"
hmd <- filter(heatmap_data, sample_type == "ESCA-adj.")
hmd$x_center <- as.numeric(as.factor(hmd$x_center))
spd <- filter(sample_pos_table, sample_type == "ESCA-adj.")
spd$x_center <- as.numeric(as.factor(spd$x_center))
htmap.size[3] <- length(spd$x_center)

plt_hmap2 <- ggplot(hmd, 
                    aes(x = x_center, y = y_center, fill = expr, 
                        height = height, width = width)) + 
  geom_tile() +
  #facet_wrap(.~sample_type)+
  scale_fill_gradient2("Rel. Abund.",trans="sqrt", high=myColor[1], mid = myColor[50], low=myColor[100], midpoint = 0.40, breaks=c(0, 0.1, 0.10, 0.30, 0.50, 0.80)) +
  scale_x_continuous(breaks = spd$x_center, 
                     labels = spd$ID, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  scale_y_continuous(breaks = gene_pos_table[, "y_center"], 
                     labels = rep("", nrow(gene_pos_table)),
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "ESCA-adj.", y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank(), #element_text(size = rel(1), hjust = 0.5,vjust=0.5, angle = 90), 
        axis.ticks.y = element_blank(),
        #axis.title.x = element_text(angle=90, hjust=0.75, vjust=0.5),
        # margin: top, right, bottom, and left
        plot.margin = unit(c(1, 0.01, 0.01, -0.7), "cm"), 
        panel.grid.minor = element_blank())

# Dendrogram plot
plt_dendr <- ggplot(segment_data) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  scale_x_reverse(expand = c(0, 0.5)) + 
  scale_y_continuous(breaks = gene_pos_table$y_center, 
                     labels = gene_pos_table$gene, 
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "", y = "", colour = "", size = "") +
  theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(1, 0.01, 0.01, -0.7), "cm"))


prntRelAbund <- avgRelAbundCutoff*100

htmap.size<-htmap.size/max(htmap.size)
htmap.size[1] <- 0.25

p <- plt_dendr+plt_hmap1+plt_hmap2+ 
  plot_layout(
    nrow=1, widths = htmap.size,
    guides="collect"
  ) +
  plot_annotation(
    title="TCGA WGS Data showing average relative abundance of genera by individual",
    subtitle=paste0("Subset to OTU average relative abundance > ",prntRelAbund,"%")
  )
p

if(save.plots == T){
  ggsave(paste0("figure1_heatmap_tcga_wgs_",save.Date,".pdf"), plot=p, units="in", width=7, height=5)
  ggsave(paste0("figure1_heatmap_tcga_wgs_",save.Date,".png"), plot=p, units="in", width=7, height=5)
}


