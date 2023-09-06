# Description ----

# ChIP-Seq and RNA-Seq data per base coverage of selected regions

# Author: MVAP
# Version: 2023-02-14

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

packages <-c("tidyr","here", "janitor", "skimr", "dplyr", "stringr", "GenomicRanges", 
             "rtracklayer", "valr", "GenomicFeatures", "extrafont")

ipak(packages)

# . Bioconductor packages ----

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}

BiocManager::install("karyoploteR", force = TRUE)

library("karyoploteR")

# . Select user font ----

loadfonts(device="win")

windowsFonts(sans = windowsFont("Arial"))

# Read Files ----

# . Load selected genes----
# One line-one gene

regions <- read.csv(
  here("input", "selectedRegions.txt"), 
  header = FALSE, 
  sep = "\t")

# . Load GFF file ----
features <- import(
  here("input", "GCF_000149955.1_ASM14995v2_genomic.gff")
  )


# . Load file with chromosome info ----
# (name, length, type: CORE/LS, etc)

genome_data <- read.csv(
  here("input", "GCF_000149955.1_ASM14995v2_assembly_report_coreLS.csv")
  )
#head(genome_data)

# Wrangling Data ----

# . clean and modify chromosome data ----
# to draw chromosomes/contig 

genome_sum <- genome_data %>%
  dplyr::select(RefSeq.Accn, Sequence.Length, Chrm.ID) %>% # select columns 
  group_by(RefSeq.Accn, Chrm.ID) %>% # group by molecule
  summarise_all(sum) %>% # calculate length per molecule
  filter(RefSeq.Accn != "na") %>% # Remove other type of chromosome molecule, ex. Mitochondrial
  mutate(start=1) %>%  # all molecules start at 1
  mutate(name= case_when(Chrm.ID != "unplaced" ~ Chrm.ID, 
                         TRUE ~ RefSeq.Accn)
         ) # add new variable columns name, if Chrm.ID is not unplaced then used Chrm.ID otherwise RefSeq.Accn

# GRanges vector with chromosome information
genome_chr<- GRanges(seqnames = Rle(genome_sum$RefSeq.Accn),
                     ranges = IRanges(start = genome_sum$start,
                                    end = genome_sum$Sequence.Length),
                     name = genome_sum$name)

# . Define cytobands ----

cyto<- genome_data %>%
  as_tibble() %>% # transformed into tibble format 
  dplyr::select(RefSeq.Accn, Sequence.Name, Chrm.ID,Sequence.Length, start, end, status) %>% # select columns
  mutate(Name.short = str_replace(Sequence.Name, "^supercont2.", "")) %>% # new column Name.short and if column Chr.ID starts with sufix replace sufix with nothing
  mutate(Name.short= ifelse(str_detect(Chrm.ID, "^Chromosome"), Chrm.ID , Name.short)) %>% # if value in column Chrm.ID start with "Chromosome then use it in name.Short 
  # to not combine unplaced contings in one line adjust start and end postions accordingly
  mutate(n_end = case_when(Chrm.ID != "unplaced" & status == "CORE" ~ end,
                           Chrm.ID != "unplaced" & status == "LS" ~ end,
                           Chrm.ID == "unplaced" & status == "LS" ~ Sequence.Length)
         ) %>%  
  mutate(s_start= case_when(Chrm.ID != "unplaced" & status == "CORE" ~ as.numeric(start),
                            Chrm.ID != "unplaced" & status == "LS" ~ as.numeric(start),
                            Chrm.ID == "unplaced" & status == "LS" ~ as.numeric(1))
         ) %>% 
  # add colors to areas in chrosomosomes
  mutate(gieStain=case_when(Sequence.Name == "gap" ~ "gneg", #gaps --> black
                            status == "CORE" ~ "gpos25", # CORE regions --> light grey
                            status == "LS" ~ "gpos75") #LS regions --> dark grey
         )%>%
  filter(RefSeq.Accn != "na") %>% # remove those witout refseqid aka mitochondrial region. 
  mutate(name= case_when(Chrm.ID != "unplaced" ~ Chrm.ID, 
                         TRUE ~ RefSeq.Accn)
         )

# GRanges vector with cytoband information

cytobands<- GRanges(seqnames = Rle(cyto$RefSeq.Accn),
                    ranges = IRanges(start = cyto$s_start, 
                                     end = cyto$n_end),
                    type = cyto$status,
                    gieStain = cyto$gieStain)


# . Obtain data to draw transcripts ----

# From gff to database txdb
txdb<- makeTxDbFromGFF(file = "input/GCF_000149955.1_ASM14995v2_genomic.gff", 
                       organism = "Fusarium oxysporum f. sp. lycopersici")


# get all genes in gff file 
feature.genes <- features[features$type == "gene"]

for (n in 1:nrow(regions)){
  
  message("GOT a new region")
  
  # get info for out gene of interest
  
  f_genes <- feature.genes[feature.genes$Name==regions[n, 1]] # get gene info
  
  f_genes_df <- as.data.frame(f_genes) # transform to obtain data
  
  # notice that regardless chain start < end therefore no need for if loop
  f_genes_df.start <- f_genes_df$start -regions[n, 2] # left position + 1000bp window
  f_genes_df.end <- f_genes_df$end + regions[n, 3] # right position + 1000 bp window
  
  # Granges object of our area of interest in the genome
  f_genes_range <- GRanges(seqnames = Rle(f_genes_df$seqnames),
                           ranges = IRanges(start = f_genes_df.start,
                                          end = f_genes_df.end))
  
  message("Got a new focus area")
  
  # get plot focus  
  kp <- plotKaryotype(genome = genome_chr, zoom = f_genes_range)
  
  # get all genes in TxDb object
  all.genes <- genes(txdb)
  
  # find genes in selected region for plotting
  res<- subsetByOverlaps(all.genes, kp$plot.region)
  
  # gene's start position
  pos_start <- f_genes_df.start
  # gene's end position
  pos_end <- f_genes_df.end
  
  # check if selected genes fir in the plotted region
  if (f_genes_df.start > min(res@ranges@start)){
    pos_start <- min(res@ranges@start)
  }
  
  if (f_genes_df.end < max(as.data.frame(res@ranges)$end)){
    pos_end <- max(as.data.frame(res@ranges)$end)
  }
  
  # Reset GRanges vector with updated positions
  f_genes_range <- GRanges(seqnames = Rle(f_genes_df$seqnames),
                           ranges = IRanges(start = pos_start,
                                          end = pos_end))
  
  # df plotting region
  kp_df <- as.data.frame(kp$plot.region)
  
  # update GRanges plotting region positive strand
  kp$plot.region <- GRanges(seqnames = Rle(kp_df$seqnames),
                           ranges = IRanges(start = kp_df$start,
                                          end = kp_df$end),
                           strand = Rle("+"))
  
  
  
  positive<-list()
  
  # Get genes in focus area and their info (transcript/exon/intron) positive strand
  positive <- makeGenesDataFromTxDb(txdb,
                                    karyoplot = kp, 
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
  
  
  # update GRanges vector plotting region negative strand
  kp$plot.region<- GRanges(seqnames = Rle(kp_df$seqnames),
                           ranges=IRanges(start= kp_df$start,
                                          end=kp_df$end),
                           strand=Rle("-"))
  
  negative <- list()
  
  # Get genes in focus area and their info (transcript/exon/intron) negative strand
  negative <- makeGenesDataFromTxDb(txdb, 
                                    karyoplot = kp, 
                                    plot.transcripts = TRUE,
                                    plot.transcripts.structure = TRUE)
  
  
  message("GOT a new database")
  
  # merged transcripts per gene per strand
  positive$genes$name <- positive$genes$gene_id
  positive<- mergeTranscripts(positive)
  
  negative$genes$name <- negative$genes$gene_id
  negative<- mergeTranscripts(negative)
  
  message("merged transcripts")  
  
  # Make figure ----
  
  # save to tiff otherwise don't run command
  tiff(here("plots",paste0(regions[n,1], "_methylation_noname_Cluster.tiff")), 
       width = 6000,
       height = 3000, 
       res = 300, 
       compression = "none")
  
  message("New file")      
  
  # adjust plot area parameters
  pp <- getDefaultPlotParams(plot.type=2)
  pp$leftmargin <- 0.15
  pp$topmargin <- 15
  pp$bottommargin <- 15
  pp$ideogramheight <- 5
  pp$data1inmargin <- 10
  pp$data2inmargin <- 10
  pp$data2height<- 100
  
  # Plot chromosome with cytobands, zoom region 
  kp <- plotKaryotype(genome = genome_chr, 
                      zoom = f_genes_range, 
                      plot.type = 2, 
                      cex = 3, 
                      plot.params = pp,
                      cytobands = cytobands, 
                      labels.plotter = NULL)
  
  # Add title (gene name)
  kpAddMainTitle(kp, 
                 regions[n, 1], 
                 cex = 3)
  
  #Add Chromosome name 
  kpAddChromosomeNames(kp, 
                       chr.names =  genome_chr[genome_chr@seqnames == f_genes_range@seqnames]$name, 
                       cex = 2)
  
  # Add ticks with position in chromosome
  kpAddBaseNumbers(kp, 
                   tick.dist = 2000, 
                   digits = 4,
                   add.units = TRUE, 
                   cex = 1, 
                   tick.len = 3)
  
  # Plot genes with introns and non-coding regions
  
  # Add positive strand (upper region)
  suppressWarnings(kpPlotGenes(kp, 
                               data = positive, 
                               r0 = 0, 
                               r1 = 0.15, 
                               add.gene.names = FALSE,
                               gene.name.position = "top",
                               gene.name.col = "black", 
                               gene.name.cex=1,
                               introns.col = "grey", 
                               coding.exons.col = "red3", 
                               non.coding.exons.height = 0.2)
                   )
  
  # Add negative strand (lower region)
  suppressWarnings(kpPlotGenes(kp,
                               data = negative, 
                               data.panel = 2, 
                               r0 = 0.1, 
                               r1 = 0.4,
                               add.gene.names = FALSE,
                               gene.name.position = "top",
                               gene.name.col = "black", 
                               gene.name.cex=1,
                               introns.col = "grey",
                               coding.exons.col = "red3", 
                               non.coding.exons.height = 0.2)
                   )
  
  
  # RNA-Seq data (bam file extension)
  
  # add left side label
  kpAddLabels(kp,
              labels = "RNA-Seq", 
              label.margin = 0.12, 
              srt = 90, 
              pos = 3, 
              cex = 3, 
              r0 = 0.6, 
              r1 = 1)
  
  # plot bam file 1
  kp <- kpPlotBAMCoverage(kp, 
                          data = "E:/workingDirectory/RafaPalos/genome_mapping/B138_FOX4287_EKRN210002787-1A_wt_0Cu2_fol4287_bwaRefseq_sort.bam", 
                          col = "darkblue", 
                          border = NA, 
                          r0 = 0.8, 
                          r1 = 1, 
                          data.panel = 1)
  
  # get max coverage value
  computed.ymax <- ceiling(kp$latest.plot$computed.values$max.coverage)
  
  # plot axis
  kpAxis(kp, 
         r0 = 0.8, 
         r1 = 1, 
         ymax = computed.ymax, 
         data.panel = 1, 
         tick.pos = computed.ymax, 
         cex = 1)
  
  # add sample label
  kpAddLabels(kp, 
              "wt -Cu", 
              r0 = 0.8, 
              r1 = 1, 
              label.margin = 0.01, 
              data.panel = 1, 
              cex = 2)
  
  # plot bam file 2
  kp <- kpPlotBAMCoverage(kp, 
                          data = "E:/workingDirectory/RafaPalos/genome_mapping/B138_FOX4287_EKRN210002793-1A_DmA_0Cu2_fol4287_bwaRefseq_sort.bam", 
                          col = "red", 
                          border = NA, 
                          r0 = 0.8, 
                          r1 = 0.6, 
                          ymax = computed.ymax, 
                          data.panel = 1)
  
  
  kpAxis(kp, 
         r0 = 0.8, 
         r1 = 0.6, 
         ymax = computed.ymax, 
         data.panel = 1, 
         tick.pos = computed.ymax, 
         cex = 1)
  
  kpAddLabels(kp,
              expression(italic('mac1')*Delta*' -Cu'), 
              r0 = 0.8,
              r1 = 0.6,  
              label.margin = 0.01,
              data.panel = 1, 
              cex=2)
  
  
  # ChIP-Seq data
  
  kpAddLabels(kp, 
              labels = "ChIP-seq", 
              label.margin = 0.12, 
              srt = 90, 
              pos = 3, 
              cex = 3, 
              r0 = 0.1, 
              r1 = 0.5)
  
  kp <- kpPlotBigWig(kp, 
                     data = "input/bigwigFiles/macAst_noCu_1.wig.bw", 
                     ymax = "visible.region", 
                     r0 = 0.3, 
                     r1 = 0.5, 
                     col = "grey", 
                     data.panel = 1)
  
  computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
  
  kpAxis(kp, 
         ymin = 0, 
         ymax = computed.ymax, 
         tick.pos = computed.ymax,
         r0 = 0.3, 
         r1 = 0.5, 
         cex = 1, 
         data.panel = 1)
  
  kpAddLabels(kp, 
              labels = expression(italic('mac1')^st*' -Cu'), 
              r0 = 0.3, 
              r1 = 0.4,
              cex = 2, 
              label.margin = 0.01, 
              data.panel = 1)
  
  
  kp <- kpPlotBigWig(kp, 
                     data = "input/bigwigFiles/macAst_plusCu_1.wig.bw", 
                     ymax = computed.ymax, 
                     r0 = 0.3, 
                     r1 = 0.1, 
                     col = "steelblue1", 
                     data.panel = 1)
  
  computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
  
  kpAxis(kp, 
         ymin = 0, 
         ymax = computed.ymax,
         tick.pos = computed.ymax,
         r0 = 0.3, 
         r1 = 0.1, 
         cex = 1, 
         data.panel = 1 )
  
  kpAddLabels(kp, 
              labels = expression(italic('mac1')^st*' +Cu'),
              r0 = 0.3, 
              r1 = 0.2,
              cex = 2, 
              label.margin = 0.01, 
              data.panel = 1)
  
  dev.off()
  message("finished figure")
}
  

