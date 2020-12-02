# ============================================= #
# script: get_data.R
# Project: ECA
# Author(s): R.N. Padgett
# ============================================= #
# Data Created: 2020-01-30
# Date Modified: 2020-10-22
# By: R. Noah Padgett
# ============================================= #
# ============================================= #
# Purpose:
# This R script is for loading and formating
#   the data file for use in analyses.

# need to make sure the microbiome code is read in
source("code/microbiome_statistics_and_functions.R")

# ================================= #
# read in NCI-UMD Data
# ================================= #
# meta-data
meta.data <- read_xlsx("data/NCI-UMD/NCI_UMD_metadata_2020_09_17.xlsx")
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

# CLEAN MICROBIOME DATA
# rarify
# rarified to an even depth of
phylo_data0 <- rarefy_even_depth(phylo.data0, replace = T, rngseed = 20200101)

# compute prevalence of each feature
prevdf <- apply(X=otu_table(phylo_data0),
                MARGIN= ifelse(taxa_are_rows(phylo_data0), yes=1, no=2),
                FUN=function(x){sum(x>0)})
# store as data.frame with labels
prevdf <- data.frame(Prevalence=prevdf,
                     TotalAbundance=taxa_sums(phylo_data0),
                     tax_table(phylo_data0))

# get totals
totals <- plyr::ddply(prevdf, "Phylum",
                      function(df1){
                        A <- cbind(mean(df1$Prevalence), sum(df1$Prevalence))
                        colnames(A) <- c("Average", "Total")
                        A
                      }
) # end


filterPhyla <- totals$Phylum[totals$Total <= 3, drop=T] # drop allows some of the attributes to be removed
phylo_data1 <- subset_taxa(phylo_data0, !Phylum %in% filterPhyla)

# Remove phylum
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(phylo_data1, "Phylum"))

# prevalence threshold
prevalenceThreshold <- 0.001*nsamples(phylo_data1)

# execute the filtering to this level
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
phylo_data2 <- prune_taxa(keepTaxa, phylo_data1)

# remove contaminated genera
a <- taxa_names(phylo_data2)
conTaxa <- c("Ralstonia", "Delftia", "Agrobacterium", "Janthinobacterium", "Halomonas", "Methylobacterium", "Aquamicrobium", "Diaphorobacter", "Herbaspirillum", "Variovorax")
i <- 1
K <- 0
for(i in 1:length(conTaxa)){
  kT <- a[a %like% conTaxa[i]]
  K <- c(K, kT)
}
b <- !a %in% K
phylo_data2 <- phyloseq::prune_taxa(b, phylo_data2)

genusNames <- get_taxa_unique(phylo_data2, "Genus")

# Merged data
phylo_data = tax_glom(phylo_data2, "Genus", NArm = TRUE)

# compute alpha diversity metrics
alpha_div <- estimate_richness(phylo_data, split = TRUE,
                               measures = c("Observed", "Shannon",
                                            "Simpson", "InvSimpson"))
alpha_div$IDi <- rownames(alpha_div) %>%
  as.factor()
meta.data <- sample_data(phylo_data) %>%
  unclass() %>%
  data.frame() %>%
  mutate(
    ID = ID,
    IDi = paste0("X",ID)
  ) %>%
  left_join(alpha_div, by = "IDi")

meta.data$caseID3 <- as.factor(meta.data$IDi)
# now, extract the information from .biom/phyloseq for other analyses
otus <- psmelt(phylo_data) # uses the CLEAN PHYLO_DATA file to just get abundance data
otus <- otus[ otus$Sample %in% meta.data$ID, ]
otus <- reshape(otus, idvar = "OTU", timevar = "Sample", direction = "wide")
OTU <- otus$OTU
Kingdom <- otus[, paste0('Kingdom.', meta.data$ID[1])]
Phylum <- otus[, paste0('Phylum.', meta.data$ID[1])]
Class <- otus[, paste0('Class.', meta.data$ID[1])]
Order <- otus[, paste0('Order.', meta.data$ID[1])]
Family <- otus[, paste0('Family.', meta.data$ID[1])]
Genus <- otus[, paste0('Genus.', meta.data$ID[1])]
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
phylo.data.nci.umd <- phylo_data

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

# next "filter out" 0 columns
# ncol(otus.RNAseq) = 173
otus.RNAseq[nrow(otus.RNAseq) + 1,] <- colSums(otus.RNAseq)
otus.RNAseq[nrow(otus.RNAseq),][otus.RNAseq[nrow(otus.RNAseq),] == 0] <- NA
otus.RNAseq <- otus.RNAseq %>%  select_if(~ !any(is.na(.)))
otus.RNAseq <- otus.RNAseq[-780, ]
# ncol(otus.RNAseq) = 66

# ncol(otus.WGS) = 139
otus.WGS[nrow(otus.WGS) + 1,] <- colSums(otus.WGS)
otus.WGS[nrow(otus.WGS),][otus.WGS[nrow(otus.WGS),] == 0] <- NA
otus.WGS <- otus.WGS %>%  select_if(~ !any(is.na(.)))
otus.WGS <- otus.WGS[-780, ]
# ncol(otus.WGS) = 123

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
phylo.data <- merge_phyloseq(otu_table(otu.tab), meta)
tax_table(phylo.data) <- otu.name

# CLEAN MICROBIOME DATA
# rarify
# rarified to an even depth of
phylo_data0 <- rarefy_even_depth(phylo.data, replace = T, rngseed = 20200101)

# compute prevalence of each feature
prevdf <- apply(X=otu_table(phylo_data0),
                MARGIN= ifelse(taxa_are_rows(phylo_data0), yes=1, no=2),
                FUN=function(x){sum(x>0)})
# store as data.frame with labels
prevdf <- data.frame(Prevalence=prevdf,
                     TotalAbundance=taxa_sums(phylo_data0),
                     tax_table(phylo_data0))

# get totals
totals <- plyr::ddply(prevdf, "Phylum",
                      function(df1){
                        A <- cbind(mean(df1$Prevalence), sum(df1$Prevalence))
                        colnames(A) <- c("Average", "Total")
                        A
                      }
) # end


filterPhyla <- totals$Phylum[totals$Total <= 3, drop=T] # drop allows some of the attributes to be removed
phylo_data1 <- subset_taxa(phylo_data0, !Phylum %in% filterPhyla)

# Remove phylum
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(phylo_data1, "Phylum"))

# prevalence threshold
prevalenceThreshold <- 0.001*nsamples(phylo_data1)

# execute the filtering to this level
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
phylo_data2 <- prune_taxa(keepTaxa, phylo_data1)

# remove contaminated genera
a <- taxa_names(phylo_data2)
conTaxa <- c("Pseudomonadales", "Comamonadaceae", "Rhizobiales", "Burkholderiales", "Paenibacillaceae", "Staphylococcus epidermidis", "Propionibacterium acnes", "Escherichia", "Bacillaceae")
i <- 1
K <- 0
for(i in 1:length(conTaxa)){
  kT <- a[a %like% conTaxa[i]]
  K <- c(K, kT)
}
b <- !a %in% K
phylo_data2 <- phyloseq::prune_taxa(b, phylo_data2)


genusNames <- get_taxa_unique(phylo_data2, "Genus")

# Merged data
phylo_data = tax_glom(phylo_data2, "Genus", NArm = TRUE)

# compute alpha diversity metrics
alpha_div <- estimate_richness(phylo_data, split = TRUE,
                               measures = c("Observed","Shannon",
                                            "Simpson", "InvSimpson"))
alpha_div$ID <- rownames(alpha_div) %>%
  as.factor()
meta_data0 <- sample_data(phylo_data) %>%
  unclass() %>%
  data.frame() %>%
  left_join(alpha_div, by = "ID")

meta_data0$SubjectID3 <- as.factor(meta_data0$ID)
# now, extract the information from .biom/phyloseq for other analyses
otus <- psmelt(phylo_data) # uses the CLEAN PHYLO_DATA file to just get abundance data
otus <- otus[ otus$Sample %in% meta.data$ID, ]
otus <- reshape(otus, idvar = "OTU", timevar = "Sample", direction = "wide")
OTU <- otus$OTU
Kingdom <- otus[, paste0('Kingdom.', meta.data$ID[1])]
Phylum <- otus[, paste0('Phylum.', meta.data$ID[1])]
Class <- otus[, paste0('Class.', meta.data$ID[1])]
Order <- otus[, paste0('Order.', meta.data$ID[1])]
Family <- otus[, paste0('Family.', meta.data$ID[1])]
Genus <- otus[, paste0('Genus.', meta.data$ID[1])]
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
                        meta.dat = meta_data0)
names(microbiome.data) <- c("otu.tab", "otu.name","abund.list",
                            "meta.dat")

microbiome.data.tcga.RNAseq <- microbiome.data
phylo.data.tcga.RNAseq <- phylo.data

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
phylo.data <- merge_phyloseq(otu_table(otu.tab), meta)
tax_table(phylo.data) <- otu.name

# CLEAN MICROBIOME DATA
# rarify
# rarified to an even depth of
phylo_data0 <- rarefy_even_depth(phylo.data, replace = T, rngseed = 20200101)

# compute prevalence of each feature
prevdf <- apply(X=otu_table(phylo_data0),
                MARGIN= ifelse(taxa_are_rows(phylo_data0), yes=1, no=2),
                FUN=function(x){sum(x>0)})
# store as data.frame with labels
prevdf <- data.frame(Prevalence=prevdf,
                     TotalAbundance=taxa_sums(phylo_data0),
                     tax_table(phylo_data0))

# get totals
totals <- plyr::ddply(prevdf, "Phylum",
                      function(df1){
                        A <- cbind(mean(df1$Prevalence), sum(df1$Prevalence))
                        colnames(A) <- c("Average", "Total")
                        A
                      }
) # end


filterPhyla <- totals$Phylum[totals$Total <= 3, drop=T] # drop allows some of the attributes to be removed
phylo_data1 <- subset_taxa(phylo_data0, !Phylum %in% filterPhyla)

# Remove phylum
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(phylo_data1, "Phylum"))

# prevalence threshold
prevalenceThreshold <- 0.001*nsamples(phylo_data1)

# execute the filtering to this level
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
phylo_data2 <- prune_taxa(keepTaxa, phylo_data1)

# remove contaminated genera
a <- taxa_names(phylo_data2)
conTaxa <- c("Pseudomonadales", "Comamonadaceae", "Rhizobiales", "Burkholderiales", "Paenibacillaceae", "Staphylococcus epidermidis", "Propionibacterium acnes", "Escherichia", "Bacillaceae")
i <- 1
K <- 0
for(i in 1:length(conTaxa)){
  kT <- a[a %like% conTaxa[i]]
  K <- c(K, kT)
}
b <- !a %in% K
phylo_data2 <- phyloseq::prune_taxa(b, phylo_data2)


genusNames <- get_taxa_unique(phylo_data2, "Genus")

# Merged data
phylo_data = tax_glom(phylo_data2, "Genus", NArm = TRUE)

# compute alpha diversity metrics
alpha_div <- estimate_richness(phylo_data, split = TRUE,
                               measures = c("Observed","Shannon",
                                            "Simpson", "InvSimpson"))
alpha_div$ID <- rownames(alpha_div) %>%
  as.factor()
meta_data0 <- sample_data(phylo_data) %>%
  unclass() %>%
  data.frame() %>%
  left_join(alpha_div, by = "ID")

meta_data0$SubjectID3 <- as.factor(meta_data0$ID)
# now, extract the information from .biom/phyloseq for other analyses
otus <- psmelt(phylo_data) # uses the CLEAN PHYLO_DATA file to just get abundance data
otus <- otus[ otus$Sample %in% meta.data$ID, ]
otus <- reshape(otus, idvar = "OTU", timevar = "Sample", direction = "wide")
OTU <- otus$OTU
Kingdom <- otus[, paste0('Kingdom.', meta.data$ID[1])]
Phylum <- otus[, paste0('Phylum.', meta.data$ID[1])]
Class <- otus[, paste0('Class.', meta.data$ID[1])]
Order <- otus[, paste0('Order.', meta.data$ID[1])]
Family <- otus[, paste0('Family.', meta.data$ID[1])]
Genus <- otus[, paste0('Genus.', meta.data$ID[1])]
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
                        meta.dat = meta_data0)
names(microbiome.data) <- c("otu.tab", "otu.name","abund.list",
                            "meta.dat")

microbiome.data.tcga.WGS <- microbiome.data
phylo.data.tcga.WGS <- phylo.data


# remove unnecessary items
remove(abund.list, biom.file, meta, meta.data, observed.counts, otu.name, otus, tree.file, Class, Family, Genus, i, Kingdom, new.packages, Order, OTU,packages, Phylum, meta.data1, meta.data2, meta.data3, otu.tab, phylo.data0, microbiome.data, otus.RNAseq, otus.WGS, meta.data.RNAseq, meta.data.WGS, totals, prevdf, prevdf1, phylo.data, phylo_data, phylo_data0, phylo_data1,phylo_data2, filterPhyla, genusNames, keepTaxa, prevalenceThreshold, alpha_div, meta_data0)

