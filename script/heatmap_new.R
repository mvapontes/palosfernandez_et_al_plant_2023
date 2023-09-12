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
stats_deg_status <- read.csv(here("output", "2023-0705_Rafa_liq-plantaR2R6_Mean_FC_adjpval_status.csv"))

# New gene names 

geneNames <- read.delim(here("input", "Genenames_Heatmap.txt"), sep = "\t")



# Liq macA/wt----

# . Wrangling data ----

# Select data based on colum name 
fc_ups <-  stats_deg_status[, grep("^geneID|wt_0Cu_VS_wt_1Cu|DmA_0Cu_VS_wt_0Cu|DmA_1Cu_VS_wt_1Cu", 
                                   colnames(stats_deg_status))
]

# Remove selected genes
fc_ups.f <- fc_ups[!grepl("FOXG_09162|FOXG_12319",
                          fc_ups$geneID),]

# Keep only regulated genes in WT -Cu / WT +Cu 
fc_ups.f <- fc_ups[fc_ups$wt_0Cu_VS_wt_1Cu_status!="noReg",]

# Delete adjusted p-value column 
fc_ups.f <- fc_ups.f[,!grepl("padj",colnames(fc_ups.f))]

# Reorder columns 
fc_ups.f <- fc_ups.f[,c(1,4,3,2,7,6,5)]


# Replace the gene id with the final name and add the column as rowname

fc_ups.df <- fc_ups.f %>%
  left_join(geneNames, by = "geneID") %>% # add tag column based on same geneID
  dplyr::select(-geneID) %>% # remove geneID column
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "Tag") # convert column to rownames 

# Reorder rows desc based on WT -Cu / Wt + Cu Log(2)Foldchange 
fc_ups.df<- fc_ups.df %>% arrange(desc(wt_0Cu_VS_wt_1Cu.log2FoldChange))


# . Heatmap ----

# Define colors. Max and Min are max and min from viridis palette
col_fun = colorRamp2(c(-15,0,5), c("#440154FF", "white","#FDE725FF"))

# Define column labels
my_titles <- c("wt -Cu /\nwt +Cu",
               expression(atop(italic('mac1')*plain(Delta)~"-Cu /","wt -Cu")),
               expression(atop(italic('mac1')*plain(Delta)~"+Cu /","wt +Cu")))


# save to file 
png("plots/2023-08-29_heatmapSelectedGenesWTcu_onlyliq_2.png", 
    width = 2000, 
    height =2500, 
    pointsize=11, 
    family = "sans", 
    res=300)

ht = Heatmap(fc_ups.df[,1:3], # select foldchange columns 
             # do not cluster or reorder rows and columns
             cluster_rows = FALSE, 
             row_dend_reorder = FALSE, 
             cluster_columns = FALSE, 
             column_dend_reorder = FALSE, 
             column_names_gp = grid::gpar(fontsize = 11, 
                                          fontfamily = "sans"),
             row_names_gp = grid::gpar(fontsize = 11, 
                                       fontfamily = "sans"),
             col = col_fun,
             width = unit(8, "cm"), 
             height= unit(15,"cm"),
             column_labels = my_titles,
             column_names_rot = 90,
             cell_fun = function(j, i, x, y, width, height, fill) {
               small_mat <- fc_ups.df[,4:6]
               if(small_mat[i, j] != "noReg")
                 grid.text("*", 
                           x, 
                           y, 
                           gp = gpar(fontsize = 15, col= "red")
                 )
             },
             heatmap_legend_param = list(direction = "horizontal",
                                         title = bquote(bold(Log[2](FC))), 
                                         legend_gp = gpar(fontsize = 11, 
                                                          fontfamily = "sans"),
                                         title_position = "lefttop",
                                         use_raster = TRUE,
                                         raster_resize_mat = TRUE, 
                                         raster_device = "png"))#,
# rect_gp = gpar(col = "black",
#                lwd = 1,
#                lty = "dotted"))

draw(ht, heatmap_legend_side = "top" )
dev.off()

# macA vs Wt ----

# . Wrangling data ----

# Select data
fc_ups <- stats_deg_status[!grepl("FOXG_09162|FOXG_12319", stats_deg_status$geneID),# remove selected genes by row
                           grep("^geneID|wt_0Cu_VS_wt_1Cu|R2_wt_VS_wt_1Cu|R6_wt_VS_wt_1Cu|DmA_0Cu_VS_wt_0Cu|DmA_1Cu_VS_wt_1Cu|R2_DmA_VS_R2_wt|R6_DmA_VS_R6_wt", colnames(stats_deg_status)) # select comparisons by column name
                           ]
# Keeg only up regulated genes
fc_ups.f <- fc_ups[fc_ups$wt_0Cu_VS_wt_1Cu_status == "upGene", ]

# Delete adjusted pvalue column
fc_ups.f <- fc_ups.f[,!grepl("padj",colnames(fc_ups.f))]

# Replace the gene id with the final name and add the column as rowname
fc_ups.df <- fc_ups.f %>%
  left_join(geneNames, by = "geneID") %>% # add tag column based on same geneID
  dplyr::select(-geneID) %>% # remove geneID column
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "Tag") # convert column to rownames 

# select user specific samples
fc_ups.df <- fc_ups.df %>% arrange(desc(wt_0Cu_VS_wt_1Cu.log2FoldChange)) %>% 
  dplyr::select(starts_with(c("R2_wt_VS_wt_1Cu", "R6_wt_VS_wt_1Cu","R2_DmA_VS_R2_wt","R6_DmA_VS_R6_wt")))

# reorder columns for graph
fc_ups.df <- fc_ups.df[,c(1, 3, 5, 7, 2, 4, 6, 8)]

# . Heatmap ----

# Define colors. Max and Min are max and min from viridis palette
col_fun = colorRamp2(c(-10,0,5), c("#440154FF", "white","#FDE725FF"))

# Define column labels
my_titles <- c("wt 2dpi /\nwt +Cu", 
               "wt 6dpi /\nwt +Cu",
               expression(atop(italic(mac1)*plain(Delta)~" 2dpi /","wt 2dpi")),
               expression(atop(italic(mac1)*plain(Delta)~" 6dpi /", "wt 6dpi")))

# save to
png("plots/2023-08-30_heatmapSelectedGenesWTDEGS_allcondition.png", 
    width = 3000,
    height =2500, 
    pointsize = 11, 
    family = "sans",
    res = 300)

# Heatmap
ht = Heatmap(fc_ups.df[,1:4],
             cluster_rows = FALSE,
             row_dend_reorder = FALSE,
#             cluster_rows = TRUE, 
#            clustering_distance_rows = "euclidean", 
#             clustering_method_rows = "complete",
#             row_dend_reorder = TRUE,
             cluster_columns = FALSE, 
             column_dend_reorder = FALSE, 
             column_names_gp = grid::gpar(fontsize = 11, 
                                          fontfamily = "sans"),
             row_names_gp = grid::gpar(fontsize = 11,
                                       fontfamily = "sans"),
             col = col_fun,
             width = unit(18, "cm"), 
             height= unit(12,"cm"),
             column_labels = my_titles,
             column_names_rot = 90,
             cell_fun = function(j, i, x, y, width, height, fill) {
               small_mat <- fc_ups.df[,5:8]
               if(small_mat[i, j] != "noReg")
                 grid.text("*", 
                           x, 
                           y, 
                           gp = gpar(fontsize = 15, 
                                     col= "red"))},
             
             heatmap_legend_param = list(direction = "horizontal",
                                         title = bquote(bold(Log[2](FC))), 
                                         legend_gp = gpar(fontsize = 11, 
                                                          fontfamily= "sans"),
                                         title_position = "lefttop",
                                         use_raster = TRUE, 
                                         raster_resize_mat = TRUE, 
                                         raster_device= "png"))#,
#             rect_gp = gpar(col = "grey", lwd = 1, lty = "dotted"))


draw(ht, heatmap_legend_side = "top" )
dev.off()
# 
# # Up regulated WT +/-Cu (all conditions) ----
# # Study expression of genes up regulated in comparison wt -Cu VS wt +Cu (wt_0Cu_VS_wt_1Cu)
# 
# # . Wrangling data ----
# 
# # Select data based on column name
# 
# fc_ups <-  stats_deg_status[, grep("^geneID|wt_0Cu_VS_wt_1Cu|R2_wt_VS_wt_1Cu|R6_wt_VS_wt_1Cu|DmA_0Cu_VS_wt_0Cu|DmA_1Cu_VS_wt_1Cu|R2_DmA_VS_R2_wt|R6_DmA_VS_R6_wt", colnames(stats_deg_status))
#                             ]
# # Select user gene
# 
# fc_ups.f <- fc_ups[grep("FOXG_11474|FOXG_17215|FOXG_03101|FOXG_18310|FOXG_02393",fc_ups$geneID),]
# 
# # Delete adjusted p-value column (padj), Keep foldchange and status
# 
# fc_ups.f <- fc_ups.f[,!grepl("padj",colnames(fc_ups.f))]
# 
# # Reorder columns, first foldchange and then adjpvalue per sample. 
# # it is necessary to plot the asterisc if the comparison has adjpvalue < 0.05
# 
# fc_ups.f <- fc_ups.f[,c(1, 8, 3, 7, 5, 4, 2, 6, 15, 10, 14, 12, 11, 9, 13)]
# 
# # Replace the gene id with the final name and add the column as rowname
# 
# fc_ups.df <- fc_ups.f %>%
#   left_join(geneNames, by = "geneID") %>% # add tag column based on same geneID
#   dplyr::select(-geneID) %>% # remove geneID column
#   `row.names<-`(., NULL) %>% 
#   column_to_rownames(var = "Tag") # convert column to rownames 
# 
# # . Heatmap ----
# 
# # Define colors. Max and Min are max and min from viridis palette
# col_fun = colorRamp2(c(-10,0,5), c("#440154FF", "white","#FDE725FF"))
# 
# # Define column labels
# my_titles <- c("wt -Cu /\nwt +Cu", 
#                "wt 2dpi /\nwt +Cu", 
#                "wt 6dpi /\nwt +Cu",
#                expression(atop(italic('mac1')*plain(Delta)~" -Cu /","wt -Cu")),
#                expression(atop(italic('mac1')*plain(Delta)~" +Cu /","wt +Cu")), 
#                expression(atop(italic('mac1')*plain(Delta)~" 2dpi /","wt 2dpi")),
#                expression(atop(italic('mac1')*plain(Delta)~" 6dpi /", "wt 6dpi")))
# 
# # save to 
# png("plots/2023-08-29_heatmapSelectedGenesUP_FC.png", width = 2000, height =1350, pointsize=11, family = "sans", res=300)
# 
# # Heatmap
# ht = Heatmap(fc_ups.df[,1:7], # select only Foldchange columns
#              cluster_rows = TRUE, # cluster rows complete euclidean
#              clustering_distance_rows = "euclidean",  
#              clustering_method_rows = "complete", 
#              row_dend_reorder = TRUE, # order based on the dendrogram
#              cluster_columns = FALSE, # do not cluster columns 
#              column_dend_reorder = FALSE, # keep user order
#              # row and columns name size and font style
#              column_names_gp = grid::gpar(fontsize = 11, 
#                                           fontfamily= "sans"),
#              row_names_gp = grid::gpar(fontsize = 11, 
#                                        fontfamily= "sans"),
#              col =col_fun, # grid colors
#              # width and height of the heatmap
#              width = unit(10, "cm"), 
#              height= unit(5,"cm"),
#              column_labels = my_titles, # use user labels for column names
#              # draw asterisk if adjusted p-value < 0.05
#              cell_fun = function(j, i, x, y, width, height, fill) {
#                small_mat <- fc_ups.df[,8:14] # select columns with adjusted p-value
#                if(small_mat[i, j] != "noReg") # keep only upGene or downGene
#                  grid.text("*", # draw asterisk 
#                            x, # in the middle of the grid
#                            y, 
#                            gp = gpar(fontsize = 11, # font size and color
#                                      col= "red"))},
#              heatmap_legend_param = list(direction = "horizontal", # legend horizontal
#                                          title = bquote(bold(Log[2](FC))), 
#                                          title_position = "lefttop", # title on the left of the legend
#                                          legend_gp = gpar(fontsize =11, # font and size
#                                                           fontfamily= "sans"),
#                                          use_raster = TRUE, # render as raster to reduce size
#                                          raster_resize_mat = TRUE, 
#                                          raster_device= "png")
#              )
# 
# # draw legend on top of heatmap
# draw(ht, heatmap_legend_side = "top")
# dev.off()
# 


# Liq macA/wt----

# . Wrangling data ----

# Select data based on colum name 
fc_ups <-  stats_deg_status[, grep("^geneID|wt_0Cu_VS_wt_1Cu|DmA_0Cu_VS_DmA_1Cu", 
                                   colnames(stats_deg_status))
]

# Remove selected genes
fc_ups.f <- fc_ups[!grepl("FOXG_09162|FOXG_12319",
                          fc_ups$geneID),]

# Keep only regulated genes in WT -Cu / WT +Cu 
fc_ups.f <- fc_ups[fc_ups$wt_0Cu_VS_wt_1Cu_status!="noReg",]

# Delete adjusted p-value column 
fc_ups.f <- fc_ups.f[,!grepl("padj",colnames(fc_ups.f))]

# Reorder columns 
fc_ups.f <- fc_ups.f[,c(1, 3, 2, 5, 4)]


# Replace the gene id with the final name and add the column as rowname

fc_ups.df <- fc_ups.f %>%
  left_join(geneNames, by = "geneID") %>% # add tag column based on same geneID
  dplyr::select(-geneID) %>% # remove geneID column
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "Tag") # convert column to rownames 

# Reorder rows desc based on WT -Cu / Wt + Cu Log(2)Foldchange 
fc_ups.df<- fc_ups.df %>% arrange(desc(wt_0Cu_VS_wt_1Cu.log2FoldChange))


# . Heatmap ----

# Define colors. Max and Min are max and min from viridis palette
col_fun = colorRamp2(c(-15,0,5), c("#440154FF", "white","#FDE725FF"))

# Define column labels
my_titles <- c("wt",
               expression(italic('mac1')*plain(Delta)))


# save to file 
png("plots/2023-08-29_heatmapSelectedGenesWTcu_onlyliq_2columns_2.png", 
    width = 2000, 
    height =2500, 
    pointsize=11, 
    family = "sans", 
    res=300)

ht = Heatmap(fc_ups.df[,1:2], # select foldchange columns 
             # do not cluster or reorder rows and columns
             cluster_rows = FALSE, 
             row_dend_reorder = FALSE, 
             cluster_columns = FALSE, 
             column_dend_reorder = FALSE, 
             column_names_gp = grid::gpar(fontsize = 11, 
                                          fontfamily = "sans"),
             row_names_gp = grid::gpar(fontsize = 11, 
                                       fontfamily = "sans"),
             col = col_fun,
             width = unit(6, "cm"), 
             height= unit(15,"cm"),
             column_labels = my_titles,
             column_names_rot = 360,
             cell_fun = function(j, i, x, y, width, height, fill) {
               small_mat <- fc_ups.df[,3:4]
               if(small_mat[i, j] != "noReg")
                 grid.text("*", 
                           x, 
                           y, 
                           gp = gpar(fontsize = 15, col= "red")
                 )
             },
             heatmap_legend_param = list(direction = "horizontal",
                                         title = bquote(bold(Log[2](FC-Cu/+Cu))), 
                                         legend_gp = gpar(fontsize = 11, 
                                                          fontfamily = "sans"),
                                         title_position = "lefttop",
                                         use_raster = TRUE,
                                         raster_resize_mat = TRUE, 
                                         raster_device = "png"),
 rect_gp = gpar(col = "black",
                lwd = 1,
                lty = "dotted"))

draw(ht, heatmap_legend_side = "top" )
dev.off()


# macA vs Wt ----

# . Wrangling data ----

# Select data
fc_ups <- stats_deg_status[!grepl("FOXG_09162|FOXG_12319", stats_deg_status$geneID),# remove selected genes by row
                           grep("^geneID|wt_0Cu_VS_wt_1Cu|R2_wt_VS_wt_1Cu|R6_wt_VS_wt_1Cu|R2_DmA_VS_R2_wt|R6_DmA_VS_R6_wt", colnames(stats_deg_status)) # select comparisons by column name
]
# Keeg only up regulated genes
fc_ups.f <- fc_ups[fc_ups$wt_0Cu_VS_wt_1Cu_status == "upGene", ]

# Delete adjusted pvalue column
fc_ups.f <- fc_ups.f[,!grepl("padj",colnames(fc_ups.f))]

# Replace the gene id with the final name and add the column as rowname
fc_ups.df <- fc_ups.f %>%
  left_join(geneNames, by = "geneID") %>% # add tag column based on same geneID
  dplyr::select(-geneID) %>% # remove geneID column
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "Tag") # convert column to rownames 

# select user specific samples
fc_ups.df <- fc_ups.df %>% arrange(desc(wt_0Cu_VS_wt_1Cu.log2FoldChange)) %>% 
  dplyr::select(-starts_with("wt_0Cu_VS_wt_1Cu"))

# reorder columns for graph
fc_ups.df <- fc_ups.df[,c(2, 4, 6, 8)]

# . Heatmap ----

# Define colors. Max and Min are max and min from viridis palette
col_fun = colorRamp2(c(-10,0,5), c("#440154FF", "white","#FDE725FF"))

# Define column labels
my_titles <- c("wt 2dpi /\nwt +Cu", 
               "wt 6dpi /\nwt +Cu")

# save to
png("plots/2023-09-12_heatmapSelectedGenesWTDEGS_allcondition_nuevafigura_WTCu.png", 
    width = 2000,
    height =2500, 
    pointsize = 11, 
    family = "sans",
    res = 300)

# Heatmap
ht = Heatmap(fc_ups.df[,1:2],
             cluster_rows = FALSE, 
             row_dend_reorder = FALSE,
             cluster_columns = FALSE, 
             column_dend_reorder = FALSE, 
             column_names_gp = grid::gpar(fontsize = 11, 
                                          fontfamily = "sans"),
             row_names_gp = grid::gpar(fontsize = 11,
                                       fontfamily = "sans"),
             col = col_fun,
             width = unit(6, "cm"), 
             height= unit(12,"cm"),
             column_labels = my_titles,
             column_names_rot = 90,
             cell_fun = function(j, i, x, y, width, height, fill) {
               small_mat <- fc_ups.df[,3:4]
               if(small_mat[i, j] != "noReg")
                 grid.text("*", 
                           x, 
                           y, 
                           gp = gpar(fontsize = 15, 
                                     col= "red"))},
             
             heatmap_legend_param = list(direction = "horizontal",
                                         title = bquote(bold(Log[2](FC))), 
                                         legend_gp = gpar(fontsize = 11, 
                                                          fontfamily= "sans"),
                                         title_position = "lefttop",
                                         use_raster = TRUE, 
                                         raster_resize_mat = TRUE, 
                                         raster_device= "png"),
             rect_gp = gpar(col = "grey", lwd = 1, lty = "dotted"))


draw(ht, heatmap_legend_side = "top" )
dev.off()

####  split graph ----

# Replace the gene id with the final name and add the column as rowname
fc_ups.df <- fc_ups.f %>%
  left_join(geneNames, by = "geneID") %>% # add tag column based on same geneID
  dplyr::select(-geneID) %>% # remove geneID column
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "Tag") # convert column to rownames 

# select user specific samples
fc_ups.df <- fc_ups.df %>% arrange(desc(wt_0Cu_VS_wt_1Cu.log2FoldChange)) %>% 
  dplyr::select(-starts_with("wt_0Cu_VS_wt_1Cu"))


# reorder columns for graph
fc_ups.df <- fc_ups.df[,c(1, 3, 5, 7)]

# . Heatmap ----

# Define colors. Max and Min are max and min from viridis palette
col_fun = colorRamp2(c(-10,0,5), c("#440154FF", "white","#FDE725FF"))

# Define column labels
my_titles <- c(expression(atop(italic('mac1')*plain(Delta)~" 2dpi /","wt 2dpi")),
               expression(atop(italic('mac1')*plain(Delta)~" 6dpi /", "wt 6dpi")))

# save to
png("plots/2023-09-12_heatmapSelectedGenesWTDEGS_allcondition_nuevafigura_Planta.png", 
    width = 2000,
    height =2500, 
    pointsize = 11, 
    family = "sans",
    res = 300)

# Heatmap
ht = Heatmap(fc_ups.df[,1:2],
             cluster_rows = FALSE, 
             row_dend_reorder = FALSE,
             cluster_columns = FALSE, 
             column_dend_reorder = FALSE, 
             column_names_gp = grid::gpar(fontsize = 11, 
                                          fontfamily = "sans"),
             row_names_gp = grid::gpar(fontsize = 11,
                                       fontfamily = "sans"),
             col = col_fun,
             width = unit(6, "cm"), 
             height= unit(12,"cm"),
             column_labels = my_titles,
             column_names_rot = 90,
             cell_fun = function(j, i, x, y, width, height, fill) {
               small_mat <- fc_ups.df[,3:4]
               if(small_mat[i, j] != "noReg")
                 grid.text("*", 
                           x, 
                           y, 
                           gp = gpar(fontsize = 15, 
                                     col= "red"))},
             
             heatmap_legend_param = list(direction = "horizontal",
                                         title = bquote(bold(Log[2](FC))), 
                                         legend_gp = gpar(fontsize = 11, 
                                                          fontfamily= "sans"),
                                         title_position = "lefttop",
                                         use_raster = TRUE, 
                                         raster_resize_mat = TRUE, 
                                         raster_device= "png"),
             rect_gp = gpar(col = "grey", lwd = 1, lty = "dotted"))


draw(ht, heatmap_legend_side = "top" )
dev.off()

