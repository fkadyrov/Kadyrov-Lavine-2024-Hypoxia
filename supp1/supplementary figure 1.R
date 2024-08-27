library(dplyr)
library(stringr)
library(viridis)
library(Seurat)
library(patchwork)
library(sctransform)
library(ggplot2)
library(harmony)
library(ggsci)
library(ggpubr)
library(pheatmap)
library(Matrix)
library(RColorBrewer)
library(scales)
library(data.table)
library(stats)
library(tibble)
library(tidyr)


setwd("E:/hif1a arg1 paper figures 2023/2024 REVISION/hif1a vs hif2a")

myeloid <- readRDS("./MonoTrack_Jan2023_UMAPmodel.rds")
DefaultAssay(myeloid) <- "RNA"
andrew_hif <- FetchData(myeloid, vars = c("Hif1a", "Arnt", "Epas1", "Arnt2", "Hif3a"), layer = "counts")
write.csv(andrew_hif, file = "./andrew_hif.csv", quote = FALSE)

htx <- readRDS("./new annotations hif1a_htx_arg1zsgr_d5.rds")
DefaultAssay(htx) <- "RNA"
htx_hif <- FetchData(htx, vars = c("Hif1a", "Arnt", "Epas1", "Arnt2", "Hif3a"), layer = "counts")
write.csv(htx_hif, file = "./htx_hif.csv", quote = FALSE)

rizzo <- readRDS("./rizzoint.rds")
DefaultAssay(rizzo) <- "RNA"
rizzo_hif <- FetchData(rizzo, vars = c("Hif1a", "Arnt", "Epas1", "Arnt2", "Hif3a"), layer = "counts")
write.csv(rizzo_hif, file = "./rizzo_hif.csv", quote = FALSE)










