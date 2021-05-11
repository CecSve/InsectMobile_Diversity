### load libraries ######################
library(tidyverse)
library(data.table)
library(seqinr)

# The four libraries that were sequenced from Germany in 2018 were demultiplexed, processed with DADA2. Then the output tables from the DADA2 output were merged prior to being processed with LULU. 

# first we read in the lulufied RDS file
lulified_nochim <-
  readRDS("data/Germany/lulified_nochim.RDS") # read in the lulufied RDS file for the German libraries 2018
otus <- lulified_nochim[["curated_table"]] # extract the otutable
#otus <- rownames_to_column(otus, var = "otuid")
asvtable <- otus[, grepl("X",colnames(otus))] # keep only the German samples - the ones with X*

# extract the otu IDs that are lulufied
lulufied_otus <- lulified_nochim[["curated_otus"]]

# subset the fasta file to only contain lulufied otus - get the fastafile from the DADA2 output
fasta <- read.fasta(file = "data/Germany/DADA2fw_nochim.otus", seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
fasta_sub <- fasta[names(fasta) %in% lulufied_otus]
subs <- fasta_sub[which(getLength(fasta_sub)>200)] # only select sequences with a length of more than _ bases - depends on primer, fwh = 200, art = 156, ins = 189
write.fasta(sequences = subs, names = names(subs), nbchar = 80, file.out = "data/Germany/lulufied_otus_lengthcorrected.fasta") 

# subset data to assign taxonomy
length(subs)
length(subs)
fasta_1 <- subs[c(1:5000)]
fasta_2 <- subs[c(5001:10328)]

# save the two separate fasta files to assign taxonomy with the GBIF tool (can only accept 6000 sequences at a time)
write.fasta(sequences = fasta_1, names = names(fasta_1), nbchar = 80, file.out = "data/Germany/lulufied_otus_lengthcorrected_1.fasta") 
write.fasta(sequences = fasta_2, names = names(fasta_2), nbchar = 80, file.out = "data/Germany/lulufied_otus_lengthcorrected_2.fasta") 

# the fasta file has been subset to only contain sequences with a length over 200 bp, so the asv table should only have asvs with read length above 200 bp
keep <- names(subs)
test <- asvtable
asvtable <- test[(rownames(test) %in% keep), ] 

### UPDATE FROM HERE ####################

# get metadata
data <- read.delim()

data <- data %>% filter(PCRID %in% colnames(asvtable)) # only retain metadata for the samples that are in the sequence table
keep <- data$PCRID
asvs <- asvtable[, (names(asvtable) %in% keep)] # subset metadata to match sequence data

# is there taxa in all samples?
colSums(asvs)
min(colSums(asvs))

# is some taxa not present in some samples?
rowSums(asvs)
min(rowSums(asvs))

asvs <- asvs[rowSums(asvs) > 0, ]

colSums(asvs)
min(colSums(asvs))

# taxonomy data 
taxonomy_1 <- read.delim("data/Germany/blastresult_19.csv", sep = ",")
taxonomy_2 <- read.delim("data/Germany/blastresult_20.csv", sep = ",")

# merge taxonomy data
taxonomy <- rbind(taxonomy_1, taxonomy_2)

# match asv ids to taxonomy
keep <- rownames(asvs)
taxonomy <- taxonomy[taxonomy$occurrenceId %in% keep, ]

taxa <- taxonomy
taxonomy <- taxonomy %>% separate(classification, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), "_") #split the taxonomy string into ranks using the dplyr and tidyr package

# choose only OTUs that have 99 % or higher ID match and choose only the OTUs that are assigned to class Insecta. To make sure the otuids are not deleted, we nedd to make the rownames into a column and then revert back to rownames
taxonomy_class <-
  taxonomy %>% filter(class == "Insecta") %>% column_to_rownames('occurrenceId') 

taxonomy_99 <-
  taxonomy %>% filter(class == "Insecta") %>% filter(identity >= 99) %>% column_to_rownames('occurrenceId')

# only keep insect asvs
keep <- rownames(taxonomy_class)
asvtable <- asvs[(rownames(asvs) %in% keep), ]

# save outputs
write.table(asvtable,file="cleaned-data/DE_asvtable_2018_data.txt",sep="\t")
write.table(data,file="cleaned-data/DE_metadata_2018_sequenced.txt",sep="\t")
write.table(taxonomy_class,file="cleaned-data/DE_taxonomy_Insecta.txt",sep="\t")
write.table(taxonomy_99,file="cleaned-data/DE_taxonomy_Insecta_99.txt",sep="\t")
