# Heatmaps for manuscript

# Author: MVAP
# Version: 2023-07-04

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

packages <-c("tidyverse","here", "janitor", "skimr", "dplyr", "circlize", "extrafont")

ipak(packages)

# . Bioconductor packages ----

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}

BiocManager::install("ComplexHeatmap", force = TRUE)

library('ComplexHeatmap')

# . Select user font ----

loadfonts(device="win")

windowsFonts(sans = windowsFont("Arial"))

# Load files ----

# DEG output

stats_deg_status <- read.csv(here("output", "liq-plantaR2R6_Mean_FC_adjpval_status.csv"))

# New gene names 

geneNames <- read.delim(here("input", "Genenames_Heatmap.txt"), sep = "\t")

# Liquide macA/wt -Cu/+Cu ----

# . Wrangling data ----

# Select data based on column name 
fc_ups <- stats_deg_status[!grepl("FOXG_09162|FOXG_12319", stats_deg_status$geneID), # remove selected genes by row
                           grep("^geneID|wt_0Cu_VS_wt_1Cu|DmA_0Cu_VS_DmA_1Cu", colnames(stats_deg_status))] # select columns

# Keep only rows with regulated genes in WT -Cu / WT +Cu 
fc_ups.f <- fc_ups[fc_ups$wt_0Cu_VS_wt_1Cu_status != "noReg", ]

# Delete adjusted p-value column 
fc_ups.df <- fc_ups.f[, !grepl("padj", colnames(fc_ups.f))]

# Reorder columns 
fc_ups.rdf <- fc_ups.df[,c(1, 3, 2, 5, 4)]


# Replace the gene id with the final name and add the column as rowname

fc_ups.rdfn <- fc_ups.rdf %>%
  left_join(geneNames, by = "geneID") %>% # add tag column based on same geneID
  dplyr::select(-geneID) %>% # remove geneID column
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "Tag") # convert column to rownames 

# Reorder rows desc based on WT -Cu / Wt + Cu Log(2)Foldchange 
fc_ups.rdfn<- fc_ups.rdfn %>% arrange(desc(wt_0Cu_VS_wt_1Cu.log2FoldChange))


# . Heatmap ----

# Define colors. Max and Min are max and min from viridis palette
col_fun = colorRamp2(c(-15,0,5), c("#440154FF", "white","#FDE725FF"))

# Define column labels
my_titles <- c("wt",
               expression(italic('mac1')*plain(Delta)))


# save to file 
png("plots/liq_macwt_01Cu.png", 
    width = 2000, 
    height =2500, 
    pointsize=11, 
    family = "sans", 
    res=300)

ht = Heatmap(fc_ups.rdfn[,1:2], # select fold change columns 
             
             # do not cluster or reorder rows and columns, Keep defined order
             cluster_rows = FALSE, 
             row_dend_reorder = FALSE, 
             cluster_columns = FALSE, 
             column_dend_reorder = FALSE, 
             
             # Label font size and family 
             column_names_gp = grid::gpar(fontsize = 11, 
                                          fontfamily = "sans"),
             row_names_gp = grid::gpar(fontsize = 11, 
                                       fontfamily = "sans"),
             
             # Color range
             col = col_fun,
             
             # Cell size
             width = unit(6, "cm"), 
             height= unit(15,"cm"),
             
             # Sample names
             column_labels = my_titles,
             column_names_rot = 360,
             
             # Add p-value =< 0.05
             cell_fun = function(j, i, x, y, width, height, fill) {
               
               small_mat <- fc_ups.rdfn[,3:4]
               
               if(small_mat[i, j] != "noReg")
                 
                 grid.text("*", 
                           x, 
                           y, 
                           gp = gpar(fontsize = 15, 
                                     col= "red")
                           )
               },
             
             # Add legend
             heatmap_legend_param = list(direction = "horizontal",
                                         title = bquote(bold(Log[2](FC-Cu/+Cu))), 
                                         legend_gp = gpar(fontsize = 11, 
                                                          fontfamily = "sans"),
                                         title_position = "lefttop",
                                         use_raster = TRUE,
                                         raster_resize_mat = TRUE, 
                                         raster_device = "png"),
             
             # Draw rectangles around cell 
             rect_gp = gpar(col = "black",
                lwd = 1,
                lty = "dotted")
             )

# draw figure with legend on top
draw(ht, heatmap_legend_side = "top" )

dev.off()


# macA and Wt plant ----

# . Wrangling data ----

# Select data
fc_ups <- stats_deg_status[!grepl("FOXG_09162|FOXG_12319", stats_deg_status$geneID), # remove selected genes by row
                           grep("^geneID|wt_0Cu_VS_wt_1Cu|R2_wt_VS_wt_1Cu|R6_wt_VS_wt_1Cu|R2_DmA_VS_R2_wt|R6_DmA_VS_R6_wt",
                                colnames(stats_deg_status)) # select columns
                           ]

# Keeg only up regulated genes
fc_ups.f <- fc_ups[fc_ups$wt_0Cu_VS_wt_1Cu_status == "upGene", ]

# Delete adjusted pvalue column
fc_ups.df <- fc_ups.f[,!grepl("padj",colnames(fc_ups.f))]

# Replace the gene id with the final name and add the column as rowname
fc_ups.dfn <- fc_ups.df %>%
  left_join(geneNames, by = "geneID") %>% # add tag column based on same geneID
  dplyr::select(-geneID) %>% # remove geneID column
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "Tag") # convert column to rownames 

# select user specific samples
fc_ups.rdfn <- fc_ups.dfn %>% arrange(desc(wt_0Cu_VS_wt_1Cu.log2FoldChange)) %>% 
  dplyr::select(-starts_with("wt_0Cu_VS_wt_1Cu"))

# . WT samples ----

# Select wt samples and reorder columns
fc_ups.rdfn_wt <- fc_ups.rdfn[,c(2, 4, 6, 8)]

# . . Heatmap ----

# Define colors. Max and Min are max and min from viridis palette
col_fun = colorRamp2(c(-10,0,5), c("#440154FF", "white","#FDE725FF"))

# Define column labels
my_titles <- c("wt 2dpi /\nwt +Cu", 
               "wt 6dpi /\nwt +Cu")

# save to
png("plots/plant_wt_1Cu.png", 
    width = 2000,
    height =2500, 
    pointsize = 11, 
    family = "sans",
    res = 300)

# Heatmap
htwt = Heatmap(fc_ups.rdfn_wt[,1:2],
             
             # do not cluster or reorder rows and columns, Keep defined order
             cluster_rows = FALSE, 
             row_dend_reorder = FALSE,
             cluster_columns = FALSE, 
             column_dend_reorder = FALSE, 
             
             # Label font size and family
             column_names_gp = grid::gpar(fontsize = 11, 
                                          fontfamily = "sans"),
             row_names_gp = grid::gpar(fontsize = 11,
                                       fontfamily = "sans"),
             
             # Color range
             col = col_fun,
             
             # Cell size
             width = unit(6, "cm"), 
             height= unit(12,"cm"),
             
             #Sample names
             column_labels = my_titles,
             column_names_rot = 90,
             
             # Add p-value =< 0.05
             cell_fun = function(j, i, x, y, width, height, fill) {
               
               small_mat <- fc_ups.rdfn_wt[,3:4]
               
               if(small_mat[i, j] != "noReg")
                 grid.text("*", 
                           x, 
                           y, 
                           gp = gpar(fontsize = 15, 
                                     col= "red")
                           )
               },
             
             # Add legend
             heatmap_legend_param = list(direction = "horizontal",
                                         title = bquote(bold(Log[2](FC))), 
                                         legend_gp = gpar(fontsize = 11, 
                                                          fontfamily= "sans"),
                                         title_position = "lefttop",
                                         use_raster = TRUE, 
                                         raster_resize_mat = TRUE, 
                                         raster_device= "png"),
             
             # Draw rectangle around cell
             rect_gp = gpar(col = "grey",
                            lwd = 1, 
                            lty = "dotted")
             )

# draw figure with legend on top
draw(htwt, heatmap_legend_side = "top")

dev.off()

# . macA samples ----

# Select macA samples and reorder columns
fc_ups.rdfn_mac <- fc_ups.rdfn[,c(1, 3, 5, 7)]

# . . Heatmap ----

# Define colors. Max and Min are max and min from viridis palette
col_fun = colorRamp2(c(-10,0,5), c("#440154FF", "white","#FDE725FF"))

# Define column labels
my_titles <- c(expression(atop(italic('mac1')*plain(Delta)~" 2dpi /","wt 2dpi")),
               expression(atop(italic('mac1')*plain(Delta)~" 6dpi /", "wt 6dpi")))

# save to
png("plots/plant_macA_wt.png", 
    width = 2000,
    height =2500, 
    pointsize = 11, 
    family = "sans",
    res = 300)

# Heatmap
htmac = Heatmap(fc_ups.rdfn_mac[,1:2],
             
             # do not cluster or reorder rows and columns, Keep defined order
             cluster_rows = FALSE, 
             row_dend_reorder = FALSE,
             cluster_columns = FALSE, 
             column_dend_reorder = FALSE, 
             
             # Label font size and family
             column_names_gp = grid::gpar(fontsize = 11, 
                                          fontfamily = "sans"),
             row_names_gp = grid::gpar(fontsize = 11,
                                       fontfamily = "sans"),
             # Color range
             col = col_fun,
             
             # Cell size
             width = unit(6, "cm"), 
             height= unit(12,"cm"),
             
             # Sample names
             column_labels = my_titles,
             column_names_rot = 90,
             
             # Add p-value =< 0.05
             cell_fun = function(j, i, x, y, width, height, fill) {
               
               small_mat <- fc_ups.rdfn_mac[,3:4]
               
               if(small_mat[i, j] != "noReg")
                 grid.text("*", 
                           x, 
                           y, 
                           gp = gpar(fontsize = 15, 
                                     col= "red")
                           )
               },
             
             # Add legend
             heatmap_legend_param = list(direction = "horizontal",
                                         title = bquote(bold(Log[2](FC))), 
                                         legend_gp = gpar(fontsize = 11, 
                                                          fontfamily= "sans"),
                                         title_position = "lefttop",
                                         use_raster = TRUE, 
                                         raster_resize_mat = TRUE, 
                                         raster_device= "png"),
             
             # Draw rectangle around cell
             rect_gp = gpar(col = "grey",
                            lwd = 1, 
                            lty = "dotted")
             )

# draw figure with legend on top
draw(htmac, heatmap_legend_side = "top")

dev.off()

