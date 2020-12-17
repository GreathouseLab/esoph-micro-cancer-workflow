# ============================================= #
# script: get_data.R
# Project: Esoph-Microbiome-Study
# Author(s): K. L. Greathouse et al.
# ============================================= #
# Data Created: 2020-09-17
# Date Modified: 2020-09-17
# By: R. Noah Padgett
# ============================================= #
# ============================================= #
# Purpose:
# This R script is for loading and formating
#   the data file for use in analyses.
#
# ============================================= #

# need to make sure the microbiome code is read in
source("code/microbiome_statistics_and_functions.R")

# ================================= #
# read in NCI-UMD Data
# ================================= #
# meta-data
meta.data <- read_excel(
  "data/NCI-UMD/UMD Esoph dataset from EB_2019_08_06_AV edits.xlsx", 
  sheet = "FOR STATA"
)
# subset to unique "sample ids
meta.data <- meta.data %>% distinct(`Sample ID`, .keep_all = T)
meta.data$sampleid <- meta.data$`Sample ID`
  # change "_" in sampleid to "." to match .biome file
  meta.data$ID <- stringr::str_replace_all(meta.data$sampleid, "_", ".")

# get microbiome data
biom.file  <- import_biom("data/NCI-UMD/otu_table_even500.biom")
tree.file  <- read_tree("data/NCI-UMD/reps_even500.tre")

# create phyloseq object
meta <- sample_data(meta.data)
sample_names(meta) <- meta.data$ID
phylo.data0 <- merge_phyloseq(biom.file, tree.file, meta)
colnames(phylo.data0@tax_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
# phylo_data <- rarefy_even_depth(phylo_data, sample.size = 7004, rngseed=20191210)
# phylo object ready for phyloseq related analyses (alphea beta, etc..)
phylo_data <- phylo.data0

# now, extract the information from .biom/phyloseq for other analyses
otus <- psmelt(biom.file) # uses the .biom file to just get abundance data
otus <- otus[ otus$Sample %in% meta.data$ID, ]
otus <- reshape(otus, idvar = "OTU", timevar = "Sample", direction = "wide")
OTU <- otus$OTU
Kingdom <- otus[, paste0('Rank1.', meta.data$ID[1])]
Phylum <- otus[, paste0('Rank2.', meta.data$ID[1])]
Class <- otus[, paste0('Rank3.', meta.data$ID[1])]
Order <- otus[, paste0('Rank4.', meta.data$ID[1])]
Family <- otus[, paste0('Rank5.', meta.data$ID[1])]
Genus <- otus[, paste0('Rank6.', meta.data$ID[1])]
# Next get just the observed counts of the OTUs
observed.counts <- otus[, grep("Abundance",colnames(otus)) ]
i <- 1
for(i in 1:length(colnames(observed.counts))){
  colnames(observed.counts)[i] <- substring(colnames(observed.counts)[i], 11)
}
# OTU Names object
otu.name <- cbind(Kingdom,Phylum,Class,Order,Family,Genus)
rownames(otu.name) <- OTU
# Abundance List Objects
abund.list <- list( cbind(Kingdom,observed.counts),
                    cbind(Phylum,observed.counts),
                    cbind(Class,observed.counts),
                    cbind(Order,observed.counts),
                    cbind(Family,observed.counts),
                    cbind(Genus,observed.counts),
                    cbind(OTU,observed.counts))
names(abund.list) <- c("Kingdom","Phylum","Class","Order","Family","Genus","OTU")
i <- 1
for(i in 1:7){
  abund.list[[i]] <- abundance_list_create(abund.list[[i]],abund.list[[i]][,1])
}


microbiome.data <- list(otu.tab = abund.list$OTU,
                        otu.name,abund.list,
                        meta.dat = meta.data,
                        tree = tree.file)
names(microbiome.data) <- c("otu.tab", "otu.name","abund.list",
                                  "meta.dat","tree")

microbiome.data.nci.umd <- microbiome.data
phylo.data.nci.umd <- phylo.data0

# ================================== #
# Read in TCGA Data
# ================================== #
# need to merge 3 files of meta-data
meta.data1 <- read_xlsx("data/TCGA/tcga_clincal_metadata.xlsx")
meta.data2 <- read_xlsx("data/TCGA/tcga_sample_summary.xlsx")
meta.data3 <- read_xlsx("data/TCGA/tcga_samples_and_metadata.xlsx")

meta.data <- full_join(meta.data1, meta.data2, by="Patient_ID")
meta.data <- full_join(meta.data, meta.data3, by="Patient_ID")

meta.data$ID <- meta.data$UniqueID
meta.data$ID <- stringr::str_replace_all(meta.data$ID, "-", ".")

meta.data.RNAseq <- filter(meta.data, ID %like% "RNAseq")
meta.data.WGS <- filter(meta.data, ID %like% "WGS")

# otu tables for RNAseq and WGS
otus <-as.data.frame(read_xlsx("data/TCGA/tcga_otu_counts_species.xlsx"))
rownames(otus) <- otus[,1]
otus <- otus[,-1]

otus.RNAseq <- otus[, colnames(otus) %like% "RNAseq"]
otus.WGS <- otus[, colnames(otus) %like% "WGS"]

# # next "filter out" 0 columns
# # ncol(otus.RNAseq) = 173
# otus.RNAseq[nrow(otus.RNAseq) + 1,] <- colSums(otus.RNAseq)
# otus.RNAseq[nrow(otus.RNAseq),][otus.RNAseq[nrow(otus.RNAseq),] == 0] <- NA
# otus.RNAseq <- otus.RNAseq %>%  select_if(~ !any(is.na(.)))
# otus.RNAseq <- otus.RNAseq[-780, ]
# # ncol(otus.RNAseq) = 66
# 
# # ncol(otus.WGS) = 139
# otus.WGS[nrow(otus.WGS) + 1,] <- colSums(otus.WGS)
# otus.WGS[nrow(otus.WGS),][otus.WGS[nrow(otus.WGS),] == 0] <- NA
# otus.WGS <- otus.WGS %>%  select_if(~ !any(is.na(.)))
# otus.WGS <- otus.WGS[-780, ]
# # ncol(otus.WGS) = 123

# rename columns
colnames(otus.RNAseq) <- stringr::str_replace_all(colnames(otus.RNAseq), "-", ".")
colnames(otus.WGS) <- stringr::str_replace_all(colnames(otus.WGS), "-", ".")

# ============================== #
# Create data objects for RNAseq
# subset sample to only those above
meta.data <- filter(meta.data.RNAseq, ID %in% colnames(otus.RNAseq))
meta <- sample_data(meta.data)
sample_names(meta) <- meta.data$ID
otu.tab <- otu_table(otus.RNAseq, taxa_are_rows = T) # OTU table
phylo.data0 <- merge_phyloseq(otu.tab, meta)

otus <- psmelt(phylo.data0)
otus <- otus[, c("Sample", "OTU", "Abundance")]
otus <- otus[ otus$Sample %in% meta.data.RNAseq$ID, ]
# get individual taxonomy levels (kingdom, phylum, etc.)c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
otus$Rank1 <- NA
otus$Rank2 <- NA
otus$Rank3 <- NA
otus$Rank4 <- NA
otus$Rank5 <- NA
otus$Rank6 <- NA
otus$Rank7 <- NA

otus[,4:10] <- stringr::str_split(otus$OTU, ";", simplify = T)
otus <- reshape(otus, idvar = "OTU", timevar = "Sample", direction = "wide")

OTU <- otus$OTU
Kingdom <- otus[, paste0('Rank1.', meta.data$ID[1])]
Phylum <- otus[, paste0('Rank2.', meta.data$ID[1])]
Class <- otus[, paste0('Rank3.', meta.data$ID[1])]
Order <- otus[, paste0('Rank4.', meta.data$ID[1])]
Family <- otus[, paste0('Rank5.', meta.data$ID[1])]
Genus <- otus[, paste0('Rank6.', meta.data$ID[1])]
# Next get just the observed counts of the OTUs
observed.counts <- otus[, grep("Abundance",colnames(otus)) ]
i <- 1
for(i in 1:length(colnames(observed.counts))){
  colnames(observed.counts)[i] <- substring(colnames(observed.counts)[i], 11)
}
# OTU Names object
otu.name <- cbind(Kingdom,Phylum,Class,Order,Family,Genus)
rownames(otu.name) <- OTU
# Abundance List Objects
abund.list <- list( cbind(Kingdom,observed.counts),
                    cbind(Phylum,observed.counts),
                    cbind(Class,observed.counts),
                    cbind(Order,observed.counts),
                    cbind(Family,observed.counts),
                    cbind(Genus,observed.counts),
                    cbind(OTU,observed.counts))
names(abund.list) <- c("Kingdom","Phylum","Class","Order","Family","Genus","OTU")
i <- 1
for(i in 1:7){
  abund.list[[i]] <- abundance_list_create(abund.list[[i]],abund.list[[i]][,1])
}


microbiome.data <- list(otu.tab = abund.list$OTU,
                        otu.name, abund.list,
                        meta.dat = meta.data)
names(microbiome.data) <- c("otu.tab", "otu.name","abund.list","meta.dat")

# add taxa_table (otu.name object)
otu.name <- tax_table(otu.name)
phylo.data0 <- merge_phyloseq(otu.tab, otu.name, meta)

microbiome.data.tcga.RNAseq <- microbiome.data
phylo.data.tcga.RNAseq <- phylo.data0

# ============================== #
# Create data objects for WGS
# subset sample to only those above
meta.data <- filter(meta.data.WGS, ID %in% colnames(otus.WGS))
meta <- sample_data(meta.data)
sample_names(meta) <- meta.data$ID
otu.tab <- otu_table(otus.WGS, taxa_are_rows = T) # OTU table
phylo.data0 <- merge_phyloseq(otu.tab, meta)

otus <- psmelt(phylo.data0)
otus <- otus[, c("Sample", "OTU", "Abundance")]
otus <- otus[ otus$Sample %in% meta.data.WGS$ID, ]
# get individual taxonomy levels (kingdom, phylum, etc.)c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
otus$Rank1 <- NA
otus$Rank2 <- NA
otus$Rank3 <- NA
otus$Rank4 <- NA
otus$Rank5 <- NA
otus$Rank6 <- NA
otus$Rank7 <- NA

otus[,4:10] <- stringr::str_split(otus$OTU, ";", simplify = T)
otus <- reshape(otus, idvar = "OTU", timevar = "Sample", direction = "wide")

OTU <- otus$OTU
Kingdom <- otus[, paste0('Rank1.', meta.data$ID[1])]
Phylum <- otus[, paste0('Rank2.', meta.data$ID[1])]
Class <- otus[, paste0('Rank3.', meta.data$ID[1])]
Order <- otus[, paste0('Rank4.', meta.data$ID[1])]
Family <- otus[, paste0('Rank5.', meta.data$ID[1])]
Genus <- otus[, paste0('Rank6.', meta.data$ID[1])]
# Next get just the observed counts of the OTUs
observed.counts <- otus[, grep("Abundance",colnames(otus)) ]
i <- 1
for(i in 1:length(colnames(observed.counts))){
  colnames(observed.counts)[i] <- substring(colnames(observed.counts)[i], 11)
}
# OTU Names object
otu.name <- cbind(Kingdom,Phylum,Class,Order,Family,Genus)
rownames(otu.name) <- OTU
# Abundance List Objects
abund.list <- list( cbind(Kingdom,observed.counts),
                    cbind(Phylum,observed.counts),
                    cbind(Class,observed.counts),
                    cbind(Order,observed.counts),
                    cbind(Family,observed.counts),
                    cbind(Genus,observed.counts),
                    cbind(OTU,observed.counts))
names(abund.list) <- c("Kingdom","Phylum","Class","Order","Family","Genus","OTU")
i <- 1
for(i in 1:7){
  abund.list[[i]] <- abundance_list_create(abund.list[[i]],abund.list[[i]][,1])
}


microbiome.data <- list(otu.tab = abund.list$OTU,
                        otu.name, abund.list,
                        meta.dat = meta.data)
names(microbiome.data) <- c("otu.tab", "otu.name","abund.list","meta.dat")

# add taxa_table (otu.name object)
otu.name <- tax_table(otu.name)
phylo.data0 <- merge_phyloseq(otu.tab, otu.name, meta)

microbiome.data.tcga.WGS <- microbiome.data
phylo.data.tcga.WGS <- phylo.data0


# remove unnecessary items
remove(abund.list, biom.file, meta, meta.data, observed.counts, otu.name, otus, tree.file, Class, Family, Genus, i, Kingdom, new.packages, Order, OTU,packages, Phylum, meta.data1, meta.data2, meta.data3, otu.tab, phylo.data0, microbiome.data, otus.RNAseq, otus.WGS, meta.data.RNAseq, meta.data.WGS)
