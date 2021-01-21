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
# For the full dendrogram
library(plyr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggdendro)
library(gridExtra)
library(dendextend)
library(cowplot)

plot.dat <- dat.wgs %>% filter(sample_type != "0", Abundance > 0.01) %>%
  mutate(ID = as.factor(Patient_ID)) %>%
  select(sample_type, Phylum, Genus, ID, Abundance)
plot.dat <- na.omit(plot.dat) 

a <- plot.dat %>%
  expand(ID, Genus)

a$Abundance <- 0

plot.dat <- full_join(plot.dat, a)



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

A <- str_split(heatmap_data$sample, "_")

heatmap_data$ID <- heatmap_data$sample_type <- "0"

for(i in 1:nrow(heatmap_data)){
  heatmap_data$ID[i] <- A[[i]][1]
  heatmap_data$sample_type[i] <- A[[i]][2]
}

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

# Heatmap plot
idlen <- numeric(3)
# by parts
hmd <- filter(heatmap_data, sample_type == "EAC tissues w/ Barretts History")
idlen[1] <- length(unique(hmd$ID))
hmd$x_center <- as.numeric(as.factor(hmd$x_center))
spd <- filter(sample_pos_table, sample_type == "EAC tissues w/ Barretts History")
spd$x_center <- as.numeric(as.factor(spd$x_center))

plt_hmap1 <- ggplot(hmd, 
                   aes(x = x_center, y = y_center, fill = expr, 
                       height = height, width = width)) + 
  geom_tile() +
  #facet_wrap(.~sample_type)+
  scale_fill_gradient2("Abundance",trans="sqrt", high = "darkred", low = "darkblue", breaks=c(0, 0.1, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = spd$x_center, 
                     labels = spd$ID, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  scale_y_continuous(breaks = gene_pos_table[, "y_center"], 
                     labels = rep("", nrow(gene_pos_table)),
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "EAC tissues", y = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = rel(1), hjust = 0.5,vjust=0.5, angle = 90), 
        # margin: top, right, bottom, and left
        #axis.ticks.y = element_blank(),
        #axis.title.x = element_text(angle=90, hjust=0.9, vjust=0.30),
        plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"), 
        panel.grid = element_blank(),
        legend.position = "none")

# Part 2: "EAC-adjacent tissue w/ Barretts History"
hmd <- filter(heatmap_data, sample_type == "EAC-adjacent tissue w/ Barretts History")
idlen[2] <- length(unique(hmd$ID))
hmd$x_center <- as.numeric(as.factor(hmd$x_center))
spd <- filter(sample_pos_table, sample_type == "EAC-adjacent tissue w/ Barretts History")
spd$x_center <- as.numeric(as.factor(spd$x_center))

plt_hmap2 <- ggplot(hmd, 
                    aes(x = x_center, y = y_center, fill = expr, 
                        height = height, width = width)) + 
  geom_tile() +
  #facet_wrap(.~sample_type)+
  scale_fill_gradient2("expr",trans="sqrt", high = "darkred", low = "darkblue", breaks=c(0, 0.1, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = spd$x_center, 
                     labels = spd$ID, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  scale_y_continuous(breaks = gene_pos_table[, "y_center"], 
                     labels = rep("", nrow(gene_pos_table)),
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "EAC-adjacent", y = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = rel(1), hjust = 0.5,vjust=0.5, angle = 90), 
        axis.ticks.y = element_blank(),
        #axis.title.x = element_text(angle=80),
        # margin: top, right, bottom, and left
        plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"), 
        panel.grid.minor = element_blank(),
        legend.position = "none")

# Part 3: "EAC tissues w/o Barretts History"
hmd <- filter(heatmap_data, sample_type == "EAC tissues w/o Barretts History")
idlen[3] <- length(unique(hmd$ID))
hmd$x_center <- as.numeric(as.factor(hmd$x_center))
spd <- filter(sample_pos_table, sample_type == "EAC tissues w/o Barretts History")
spd$x_center <- as.numeric(as.factor(spd$x_center))

plt_hmap3 <- ggplot(hmd, 
                    aes(x = x_center, y = y_center, fill = expr, 
                        height = height, width = width)) + 
  geom_tile() +
  #facet_wrap(.~sample_type)+
  scale_fill_gradient2("Abundance",trans="sqrt", high = "darkred", low = "darkblue", breaks=c(0, 0.1, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = spd$x_center, 
                     labels = spd$ID, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  scale_y_continuous(breaks = gene_pos_table[, "y_center"], 
                     labels = rep("", nrow(gene_pos_table)),
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "No Barretts History", y = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = rel(1), hjust = 0.75, vjust=0.5, angle = 90), 
        axis.ticks.y = element_blank(),
        #axis.title.x = element_text(angle=90, hjust=0.25, vjust=0.30),
        # margin: top, right, bottom, and left
        plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"), 
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
        axis.ticks.x = element_blank())


p <- plot_grid(plt_dendr, plt_hmap1, plt_hmap2,plt_hmap3, align = 'h', rel_widths = c(0.5, 0.3, 0.3, 1),ncol = 4)
p
ggsave("output/slide-7-heatmap-01.pdf", plot=p, units="in", width=15, height=8)
ggsave("output/slide-7-heatmap-01.png", plot=p, units="in", width=15, height=8)
