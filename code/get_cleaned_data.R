# ============================================= #
# script: get_data.R
# Project: Fiber Intervention Study
# Author(s): R.N. Padgett
# ============================================= #
# Data Created: 2020-01-30
# Date Modified: 2020-04-10
# By: R. Noah Padgett
# ============================================= #
# ============================================= #
# Purpose:
# This R script is for loading and formating
#   the data file for use in analyses.

# need to make sure the microbiome code is read in
source("code/microbiome_statistics_and_functions.R")

# read in analysis
# from data-for-analaysis tab of
analysis_data <- readxl::read_xlsx("data/analysis-data/Data_Fiber_2020_02_06.xlsx", sheet="DataforAnalysis",na = ".")
colnames(analysis_data) <- c("SubjectID", colnames(analysis_data)[2:72])

hei_data <- readxl::read_xlsx("data/analysis-data/FFQ_HEI_Fiber.xlsx")
asa24_data <- readxl::read_xlsx("data/analysis-data/ASA24_fiber_F_V_intake_by_week.xlsx")
asa24_data$Week <- fct_recode(factor(asa24_data$RecallNo), "1"="1", "4"="2", "8"="3", "12"="4")
#asa24_data$Week <- factor(recode(asa24_data$RecallNo, `1`=1, `2`=4, `3`=8, `4`=12))

# supplement IDS
# *highest intake = 7, supplement
# highest cookie intake = 3
# Gut health = 1 = no issues
# 2 = some discomfort
# 3 = sig discomfort/constipation
# 4 = sig discomfort/diarrhea


# get microbiome data
biom_file  <- import_biom("data/analysis-data/OTU_Table.biom")
tree_file <- read_tree("data/analysis-data/OTU_Table.tre")
meta_data <- read.table("data/analysis-data/Metadata file_nomiss_microbiomeIDs.txt", sep="\t", header=T)
  colnames(meta_data)[1] <- "ID"
  rownames(meta_data) <- meta_data$ID
  meta_data$SubjectID2 <- as.numeric(as.factor(meta_data$SubjectID))
  meta_data$Week <- as.factor(meta_data$Week)
  ## merge meta_dat with analysis data, hei data and ASA24 data
  meta_data <- left_join(meta_data, analysis_data) #, by="SubjectID")
  meta_data <- left_join(meta_data, hei_data) #, by="SubjectID")
  meta_data <- left_join(meta_data, asa24_data) #, by=c("SubjectID","Week"))
  #colnames(meta_data) <-  str_remove(colnames(meta_data), ".x" )
  #colnames(meta_data) <-  str_remove(colnames(meta_data), ".y" )

meta <- sample_data(meta_data)
sample_names(meta) <- meta_data$ID
phylo_data0 <- merge_phyloseq(biom_file, tree_file, meta)
colnames(phylo_data0@tax_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
# phylo_data <- rarefy_even_depth(phylo_data, sample.size = 7004, rngseed=20191210)
# phylo object ready for phyloseq related analyses (alphe, beta, etc..)

## The following does everything necessary from the data processing file
# to get the data ready for analysis

# rarify
# rarified to an even depth of
phylo_data0 <- rarefy_even_depth(phylo_data0, replace = T, rngseed = 20200101)

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


filterPhyla <- totals$Phylum[totals$Total <= 5, drop=T] # drop allows some of the attributes to be removed
phylo_data1 <- subset_taxa(phylo_data0, !Phylum %in% filterPhyla)

# Remove phylum
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(phylo_data1, "Phylum"))

# prevalence threshold
prevalenceThreshold <- 0.01*nsamples(phylo_data1)

# execute the filtering to this level
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
phylo_data2 <- prune_taxa(keepTaxa, phylo_data1)

genusNames <- get_taxa_unique(phylo_data2, "Genus")

# Merged data
phylo_data = tax_glom(phylo_data2, "Genus", NArm = TRUE)

# compute alpha diversity metrics

alpha_div <- estimate_richness(phylo_data, split = TRUE,
                               measures = c("Observed", "Chao1",
                                            "ACE", "Shannon",
                                            "Simpson", "InvSimpson",
                                            "Fisher"))

alpha_div$ID <- rownames(alpha_div) %>%
  as.factor()
meta_data0 <- sample_data(phylo_data) %>%
  unclass() %>%
  data.frame() %>%
  left_join(alpha_div, by = "ID")

meta_data0$SubjectID3 <- as.factor(meta_data0$SubjectID)

# now, extract the information from .biom/phyloseq for other analyses
otus <- psmelt(phylo_data) # uses the cleaned data in phyloseq file
otus <- otus[ otus$Sample %in% meta_data$ID, ]
otus <- reshape(otus, idvar = "OTU", timevar = "Sample", direction = "wide")
OTU <- otus$OTU
Kingdom <- otus[, paste0('Kingdom.', meta_data$ID[1])]
Phylum <- otus[, paste0('Phylum.', meta_data$ID[1])]
Class <- otus[, paste0('Class.', meta_data$ID[1])]
Order <- otus[, paste0('Order.', meta_data$ID[1])]
Family <- otus[, paste0('Family.', meta_data$ID[1])]
Genus <- otus[, paste0('Genus.', meta_data$ID[1])]
# Next get just the observed counts of the OTUs
observed_counts <- otus[, grep("Abundance",colnames(otus)) ]
i <- 1
for(i in 1:length(colnames(observed_counts))){
  colnames(observed_counts)[i] <- substring(colnames(observed_counts)[i], 11)
}
# OTU Names object
otu.name <- cbind(Kingdom,Phylum,Class,Order,Family,Genus)
rownames(otu.name) <- OTU
# Abundance List Objects
abund.list <- list( cbind(Kingdom,observed_counts),
                    cbind(Phylum,observed_counts),
                    cbind(Class,observed_counts),
                    cbind(Order,observed_counts),
                    cbind(Family,observed_counts),
                    cbind(Genus,observed_counts),
                    cbind(OTU,observed_counts))
names(abund.list) <- c("Kingdom","Phylum","Class","Order","Family","Genus","OTU")
i <- 1
for(i in 1:7){
  abund.list[[i]] <- abundance_list_create(abund.list[[i]],abund.list[[i]][,1])
}


microbiome_data <- list(otu.tab = abund.list$OTU,
                        otu.name,abund.list,
                        meta.dat = meta_data0,
                        tree = tree_file)
names(microbiome_data) <- c("otu.tab", "otu.name","abund.list",
                            "meta.dat","tree")

# remove unnecessary items
remove(abund.list, analysis_data, biom_file, meta, meta_data, observed_counts, otu.name, otus, tree_file, Class, Family, Genus, i, Kingdom, new.packages, Order, OTU,packages, Phylum, phylo_data0,phylo_data1, phylo_data2, prevdf, prevdf1, totals, filterPhyla, genusNames, keepTaxa,prevalenceThreshold, alpha_div,meta_data0, hei_data)
