# Description

# RNAseq analysis pipeline.

# Author: MVAP
# Version: 2023-03-10


# Functions ----

# . ipak ----
# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  
  if (length(new.pkg)){ 
    install.packages(new.pkg, dependencies = TRUE)
  }
  sapply(pkg, require, character.only = TRUE)
}

# Packages ----

# Remember to insert %>% is ctrl + Shift +M

# . CRAN packages ----

packages <-c("tidyverse", "here", "janitor", "skimr", "devtools", "reshape2", "corrplot", "ggrepel")

ipak(packages)

# . Bioconductor packages ----

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}

BiocManager::install("DESeq2", force = TRUE)
BiocManager::install("tximport", force = TRUE)

library('DESeq2')
library('tximport')

# Load files ----

#  . Read Transcript-gene annotation user file ----

tx2gene <- read_csv(here("input", "Transcrip_gene_annotationfile.csv"), na = "NA")

tx2gene_df <- tx2gene %>% 
  dplyr::select(transcript_id, gene_id) # subset transcript and gene id for tximport package

# . Load Salmon files ----

salmon.dirs<- list.files(pattern = "^B138", ) # select directories with pattern

# List of files with the associated sample name
fileList <- tibble(file = salmon.dirs, sample = gsub(".*-1A_", "", salmon.dirs)) %>% # split file name to obtain sample name
  separate(sample, into = c("condition", "rep"), sep = -1, remove = FALSE) %>%  # separate the last character in column sample to obtain sample number, repeat number 
  mutate(batch = ifelse(str_detect(condition, "R2|R6"), 2, 1)) %>% # split by sample type liquid or plant
  mutate(across(c(condition, batch), factor)) # reorder by condition and sample type
  

files <- file.path(fileList$file, "quant.sf") # construct path to abundance .sf file 
names(files) <- fileList$sample # assigned to each sample file its name
all(file.exists(files)) # check all files exists

# . Load transcript-level abundance estimation files
txi <- tximport(files, type="salmon", tx2gene = tx2gene_df) # gene level conversion
gxi <- tximport(files, type="salmon", txOut=TRUE) # transcript level

# Explore Raw Data ----

# . Supl. Figure XX Pearson correlation raw Counts ----

x <- data.matrix(txi$counts) # corr only accepts vector or matrix variable types
corpears <- cor(x, method = "pearson")

# Save graph as png
png("plots/PearsonCorrelation_Counts.png", 
    width = 2500, 
    height = 2500, 
    res=300)

corrplot(corpears,
         order = "hclust", # order samples by cluster
         method = "color", 
         is.corr = FALSE, # use min and max cor as limit
         outline = T, 
         addgrid.col = "darkgray", 
         cl.pos = "b", # legend at the bottom
         tl.col ="black", # label color
         tl.cex = 0.7, # label font size
         addCoef.col = "white", # add coeficient in white
         number.cex = 0.6, # font size coeficient
         col = colorRampPalette(c("red", "darkgreen"))(200), # vector of colors
         )

dev.off()

# Wrangling data ----

# . Extract TPMs and Counts ----

# TPMS

tpms<- as.data.frame(txi$abundance) # tpms are in abundance column at variable txi. Save as df to modify later
colnames(tpms) <- paste(colnames(tpms),"TPMs",sep="_") # add suffix to column names

# Each condition was repeated 3 times.
# calculate Mean and standard deviation (sd) of each condition

for (b in unique(fileList$condition)){
  tpms <- tpms %>% 
    replace(is.na(.), 0) %>% 
    # To calculate the mean and sd for all conditions in a dinamically way using {} and :=
    mutate("Mean_{b}_TPMs" := rowMeans(dplyr::select(., starts_with(b)), na.rm = TRUE)) %>% 
    mutate("SD_{b}_TPMs" := rowSds(as.matrix(dplyr::select(., starts_with(b))), na.rm = TRUE)) # for rowSds to work data must be matrix or similar check info
}

# keep only mean and sd TPMs samples for later
tpms_sum <- tpms %>% 
  dplyr::select(matches("Mean|SD")) %>% # select specific columns
  as_tibble(rownames = "geneID") # df to tibble


# DE analysis ----

# . Load combinations ----

experiments <- read_csv(here("input", "combinaciones.txt.csv"))
  
exp_design <- formula(~condition) # formula with variable to analyze

# Prepare data for DESeq2
dds.raw <- DESeqDataSetFromTximport(txi, colData = fileList, design = exp_design)

# Remove genes with less than 1 count across all samples
dds.raw <- dds.raw[rowSums(counts(dds.raw)) > 1,]


# . Sup. Fig XX Principal component analysis (PCA) ----

#transform the raw count data
vst.counts <- vst(dds.raw) 

rv <- rowVars(assay(vst.counts)) # Calculate variance

# select top 500 genes with higher variance
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]

# PCA analysis top 500 genes with higher variance
vst.pca <- prcomp(t(assay(vst.counts)[select, ]), scale = FALSE) # data are normalized therefore scale=FALSE

# get PC1 and PC2 %variance
proportionvariances <- ((vst.pca$sdev^2) / (sum(vst.pca$sdev^2)))*100

# Plot PCA
ggplot(as.data.frame(vst.pca$x),
       aes(x = PC1, y = PC2,label = rownames(vst.pca$x), color = fileList$condition)) +
  geom_point(size = 3) +
  geom_text_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 24)) +
  # axis names
  xlab(paste("PC1, ", round(proportionvariances[1], 2), "%")) +
  ylab(paste("PC2, ", round(proportionvariances[2], 2), "%")) +
  theme(legend.position="none") # remove legend

ggsave(here("plots", "pcavst.png"), units = "in", dpi = 300, width = 6, height = 6)


# . DEG analysis with DESeq2----

dds <- DESeq(dds.raw)

# . DESeq2 results by comparison ----

# declare new tibble
standard_deg <- tibble::rownames_to_column(tpms, "geneID") %>% 
  dplyr::select(geneID)

# Extract Foldchange and adj. pvalue per comparison of interest
for (i in 1:nrow(experiments)){
  
  l <- experiments[i,]
  result <- results(dds, contrast = c("condition", as.character(l[,1]), as.character(l[,2]))) # get data for each combination of interest with the option contrast
  
  name <- paste(l[,1], l[,2], sep = "_VS_")
  cnames <- paste(name, colnames(result), sep = "_") # update result column names
  
  j <- cbind(rownames(result), result[,2], result[,6]) # extract and combine gene ID, Foldchange and adj. pvalue columns from DESeq2 data type
  
  colnames(j) <- c("geneID", cnames[2], cnames[6]) # assigned updated column names
  
  standard_deg <- standard_deg %>% left_join(as_tibble(j), by = "geneID") # join data
}

# Combine mean and sd with foldchange and adj. pvalue by gene id
stats_deg <- tpms_sum %>% 
  left_join(standard_deg, copy = TRUE, by = "geneID")

# . DEGs regulation ----

FC <- 2 # user input positive foldchange minimum threshold
pval <- 0.05 # user input adjusted pvalue maximum threshold
TPM <- 2 # user input positive TPMs minimum threshold

# declare new tibble
deg_status <- tibble::rownames_to_column(tpms, "geneID") %>% 
  dplyr::select(geneID)

# Determine if a gene is up, down or not regulated based on TPM, Foldchange and adj. pvalue
for (i in 1:nrow(experiments)){
  
  l <- experiments[i,]
  
  # define column names
  n1 <- paste("Mean", l[,1], "TPMs", sep = "_")
  n2 <- paste("Mean", l[,2], "TPMs", sep = "_")
  name <- paste(l[,1], l[,2], sep = "_VS_")
  status <- paste(name, "status", sep = "_")
  
  
  deg_s <- stats_deg[, grep (paste("geneID", n1, n2, name, sep = "|"), colnames(stats_deg))] %>% 
    # assigned status column name dinamicaly
    mutate("{status}" := case_when(
      (.[[2]] > TPM |.[[3]]>TPM) & .[[4]] > FC & .[[5]] < pval ~ "upGene", 
      (.[[2]]> TPM | .[[3]] > TPM) & -FC > .[[4]] & .[[5]] < pval ~ "downGene", 
      .default = "noReg"
      )
    )
  # join df by common column
  deg_status <- deg_status %>% 
    left_join(deg_s %>% dplyr::select(matches("geneID|status")), by = "geneID")
}

# Combined all files
stats_deg_status<- stats_deg %>%
  left_join(deg_status, by = "geneID") %>%  # combined Folchange and pvalue data with status
  left_join(tx2gene %>% 
              dplyr::select(gene_id, Chrm.ID, seqnames, strand, start_gene, end_gene, transcript_number, product) %>% # combine with gene information
              distinct(), 
            by = c("geneID"= "gene_id"))

# Save results to file
write.csv(stats_deg_status, 
          here("output", "liq-plantaR2R6_Mean_FC_adjpval_status.csv"),
          row.names = FALSE)

