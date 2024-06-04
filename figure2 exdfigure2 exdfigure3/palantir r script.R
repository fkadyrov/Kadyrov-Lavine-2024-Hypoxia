#cell cycle regression on mapped data

library(Seurat)
library(ggplot2)

myeloid <- readRDS("./WT_HIF1AKO_mapped.rds")

Idents(myeloid) <- "predicted.celltype"
ccmyeloid <- myeloid

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

library(gprofiler2)
mm.s.genes = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mm.g2m.genes = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

#assign cell cycle scores
ccmyeloid <- CellCycleScoring(object = ccmyeloid,
                              s.features = mm.s.genes, g2m.features = mm.g2m.genes, set.ident = TRUE)

head(ccmyeloid[[]])

ccmyeloid <- RunPCA(ccmyeloid, features = c(mm.s.genes, mm.g2m.genes))
DimPlot(ccmyeloid, reduction = "pca")


png(filename="./pre cell cycle regression pca.png", width=40, height=20, units="cm", res=600)
DimPlot(ccmyeloid, reduction = "pca", pt.size = 5)+
  scale_color_manual(values = c( '#0C7BDC', '#AA4499', 'orange2'))+
  coord_fixed(ratio = 1.3)
dev.off()

cc.myeloid <- ScaleData(ccmyeloid, vars.to.regress = c("S.Score", "G2M.Score"), 
                        features = rownames(ccmyeloid))

cc.myeloid <- RunPCA(cc.myeloid, features = VariableFeatures(cc.myeloid), nfeatures.print = 10)

cc.myeloid <- RunPCA(cc.myeloid, features = c(mm.s.genes, mm.g2m.genes))
DimPlot(cc.myeloid, reduction = "pca")

png(filename="./post cell cycle regression pca.png", width=20, height=20, units="cm", res=300)
DimPlot(cc.myeloid, reduction = "pca", pt.size = 5)+
  scale_color_manual(values = c( '#0C7BDC', '#AA4499', 'orange2'))+
  coord_fixed(ratio = 1.7)
dev.off()

# Save the scaled matrix
write.csv(as.matrix(cc.myeloid[["SCT"]]@scale.data), 
          file = "./myeloid_SCT_normalized.txt", quote = FALSE)

# Save the meta data
write.csv(cc.myeloid@meta.data, file = "./myeloid_meta.csv", quote = TRUE)

#input to palantir


#selecting monocyte based on monocyte z score

DefaultAssay(cc.myeloid) <- "SCT"
expdata <- GetAssayData(cc.myeloid)
zMono <- c("Plac8", "Ly6c2")
pops<-list(zMono)
#Z-Scores
z_scores<-NULL

for (i in 1:length(pops)) {
  genes <- pops[[i]]
  zz <- which(tolower(rownames(expdata)) %in% tolower(genes))
  av <- numeric(ncol(expdata))
  
  geneExp <- as.matrix(expdata[zz, ])
  geneExp <- t(scale(t(geneExp)))
  geneExp[is.nan(geneExp)] <- 0
  z_scores <- rbind(z_scores,(av + colSums(geneExp) / length(zz)))
}
cc.myeloid@meta.data$Mono_z<-z_scores[1,]

Idents(cc.myeloid) <- "condition"
cc.myeloid.wt <- subset(cc.myeloid, idents = "WT")

Idents(cc.myeloid.wt) <- "predicted.celltype"
myeloidwtmono <- subset(cc.myeloid.wt, idents = "Mono")

mono_z<- FetchData(myeloidwtmono, vars = "Mono_z") 
#CCACTACGTACATCCA-1_1


#start cell
png(filename="./palantir starting cell.png", width=20, height=20, units="cm", res=600)
DimPlot(myeloid, reduction = "ref.umap", pt.size = 3, label = FALSE,
        cells.highlight = "CCACTACGTACATCCA-1_1", cols.highlight = "red", cols = "gray", sizes.highlight = 3, order = TRUE)+
  coord_fixed()
dev.off()

#palantir results to umap

myeloid_meta <- read.csv2('./palantir_meta_data.csv', header=TRUE, sep=',', row.names = 1)
Myeloid_2 <- AddMetaData(myeloid, myeloid_meta)

Myeloid_2@meta.data$pseudotime <- as.numeric(as.character(Myeloid_2@meta.data$pseudotime))
Myeloid_2@meta.data$entropy <- as.numeric(as.character(Myeloid_2@meta.data$entropy))


myeloid_meta <- data.frame(myeloid_meta)
myeloid_meta

myeloid_meta$pseudotime <- as.numeric(myeloid_meta$pseudotime)
myeloid_meta$entropy <- as.numeric(myeloid_meta$entropy)
myeloid_meta$ClusterName <- as.character(myeloid_meta$ClusterName)



#plotting split of pseudotime entropy
#Ccl8
highlight_Ccl8 <- subset(myeloid_meta, ClusterName == "Ccl8")

png(filename="./Ccl8 palantir.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_meta, aes(x=pseudotime, y=entropy)) +
  ylim(0,1)+
  geom_point(size=3, color = "azure2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_Ccl8, alpha = 1,
             aes(x=pseudotime, y=entropy, colour = ClusterName),
             color = '#b66dff',
             size=3)+
  coord_fixed() +theme(legend.position = "right")
dev.off()

#MHCII
highlight_MHCII <- subset(myeloid_meta, ClusterName == "MHCII")

png(filename="./MHCII palantir.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_meta, aes(x=pseudotime, y=entropy))+
  ylim(0,1) +
  geom_point(size=3, color = "azure2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_MHCII, alpha = 1,
             aes(x=pseudotime, y=entropy, colour = ClusterName),
             color = '#009292',
             size=3)+
  coord_fixed() +theme(legend.position = "right")
dev.off()

#Trem2
highlight_Trem2 <- subset(myeloid_meta, ClusterName == "Trem2")

png(filename="./Trem2 palantir.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_meta, aes(x=pseudotime, y=entropy)) +
  ylim(0,1)+
  geom_point(size=3, color = "azure2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_Trem2, alpha = 1,
             aes(x=pseudotime, y=entropy, colour = ClusterName),
             color = '#920000',
             size=3)+
  coord_fixed() +theme(legend.position = "right")
dev.off()

#Gdf15
highlight_Gdf15 <- subset(myeloid_meta, ClusterName == "Gdf15")

png(filename="./Gdf15 palantir.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_meta, aes(x=pseudotime, y=entropy)) +
  ylim(0,1)+
  geom_point(size=3, color = "azure2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_Gdf15, alpha = 1,
             aes(x=pseudotime, y=entropy, colour = ClusterName),
             color = '#004949',
             size=3)+
  coord_fixed() +theme(legend.position = "right")
dev.off()

#Res
highlight_Res <- subset(myeloid_meta, ClusterName == "Res")

png(filename="./Res palantir.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_meta, aes(x=pseudotime, y=entropy)) +
  ylim(0,1)+
  geom_point(size=3, color = "azure2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_Res, alpha = 1,
             aes(x=pseudotime, y=entropy, colour = ClusterName),
             color = '#6db6ff',
             size=3)+
  coord_fixed() +theme(legend.position = "right")
dev.off()

#Arg1
highlight_Arg1 <- subset(myeloid_meta, ClusterName == "Arg1")

png(filename="./Arg1 palantir.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_meta, aes(x=pseudotime, y=entropy)) +
  ylim(0,1)+
  geom_point(size=3, color = "azure2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_Arg1, alpha = 1,
             aes(x=pseudotime, y=entropy, colour = ClusterName),
             color = '#006ddb',
             size=3)+
  coord_fixed() +theme(legend.position = "right")
dev.off()

#Prolif
highlight_Prolif <- subset(myeloid_meta, ClusterName == "Prolif")

png(filename="./Prolif palantir.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_meta, aes(x=pseudotime, y=entropy)) +
  ylim(0,1)+
  geom_point(size=3, color = "azure2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_Prolif, alpha = 1,
             aes(x=pseudotime, y=entropy, colour = ClusterName),
             color = '#490092',
             size=3)+
  coord_fixed() +theme(legend.position = "right")
dev.off()

#Mono
highlight_Mono <- subset(myeloid_meta, ClusterName == "Mono")

png(filename="./Mono palantir.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_meta, aes(x=pseudotime, y=entropy)) +
  ylim(0,1)+
  geom_point(size=3, color = "azure2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_Mono, alpha = 1,
             aes(x=pseudotime, y=entropy, colour = ClusterName),
             color = '#db6d00',
             size=3)+
  coord_fixed() +theme(legend.position = "right")
dev.off()

#IFN
highlight_IFN <- subset(myeloid_meta, ClusterName == "IFN")

png(filename="./IFN palantir.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_meta, aes(x=pseudotime, y=entropy)) +
  ylim(0,1)+
  geom_point(size=3, color = "azure2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_IFN, alpha = 1,
             aes(x=pseudotime, y=entropy, colour = ClusterName),
             color = '#ff6db6',
             size=3)+
  coord_fixed() +theme(legend.position = "right")
dev.off()

#cDC2
highlight_cDC2 <- subset(myeloid_meta, ClusterName == "cDC2")

png(filename="./cDC2 palantir.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_meta, aes(x=pseudotime, y=entropy)) +
  ylim(0,1)+
  geom_point(size=3, color = "azure2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_cDC2, alpha = 1,
             aes(x=pseudotime, y=entropy, colour = ClusterName),
             color = '#24ff24',
             size=3)+
  coord_fixed() +theme(legend.position = "right")
dev.off()

#FDL plots
fdl <- read.csv2('./fdl.csv', header=TRUE, sep=',', row.names = 1)
myeloid_fdl <- transform(merge(fdl,myeloid_meta,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
myeloid_fdl$x <- as.numeric(myeloid_fdl$x)
myeloid_fdl$y <- as.numeric(myeloid_fdl$y)
myeloid_fdl$pseudotime <- as.numeric(myeloid_fdl$pseudotime)
myeloid_fdl$entropy <- as.numeric(myeloid_fdl$entropy)
myeloid_fdl$cDC2 <- as.numeric(myeloid_fdl$cDC2)
myeloid_fdl$Trem2 <- as.numeric(myeloid_fdl$Trem2)

png(filename="./fdl entropy ggplot2.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_fdl, aes(x = x, y = y, color = entropy)) +
  geom_point(size = 3, alpha = 1)+
  scale_color_gradientn(colors=c("#0D0887FF", "#7E03A8FF", "#CC4678FF", "#F89441FF", "#F0F921FF"))+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_fixed(ratio=.8)
dev.off()

png(filename="./fdl pseudotime ggplot2.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_fdl, aes(x = x, y = y, color = pseudotime)) +
  geom_point(size = 3, alpha = 1)+
  scale_color_gradientn(colors=c("#0D0887FF", "#7E03A8FF", "#CC4678FF", "#F89441FF", "#F0F921FF"))+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_fixed(ratio=.8)
dev.off()

####cells on FDL layout
highlight_Ccl8 <- subset(myeloid_fdl, ClusterName == "Ccl8")

png(filename="./fdl ccl8 ggplot2.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_fdl, aes(x = x, y = y, color = ClusterName)) +
  geom_point(size = 3, alpha = 1, color ="azure2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_Ccl8,
             aes(x=x, y=y, colour = ClusterName),
             color = '#b66dff',
             size=3)+
  coord_fixed(ratio = .8) +theme(legend.position = "right")
dev.off()

highlight_MHCII <- subset(myeloid_fdl, ClusterName == "MHCII")

png(filename="./fdl MHCII ggplot2.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_fdl, aes(x = x, y = y, color = ClusterName)) +
  geom_point(size = 3, alpha = 1, color ="azure2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_MHCII,
             aes(x=x, y=y, colour = ClusterName),
             color = '#009292',
             size=3)+
  coord_fixed(ratio = .8) +theme(legend.position = "right")
dev.off()

highlight_Trem2 <- subset(myeloid_fdl, ClusterName == "Trem2")

png(filename="./fdl Trem2 ggplot2.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_fdl, aes(x = x, y = y, color = ClusterName)) +
  geom_point(size = 3, alpha = 1, color ="azure2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_Trem2,
             aes(x=x, y=y, colour = ClusterName),
             color = '#920000',
             size=3)+
  coord_fixed(ratio = .8) +theme(legend.position = "right")
dev.off()


highlight_Gdf15 <- subset(myeloid_fdl, ClusterName == "Gdf15")

png(filename="./fdl Gdf15 ggplot2.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_fdl, aes(x = x, y = y, color = ClusterName)) +
  geom_point(size = 3, alpha = 1, color ="azure2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_Gdf15,
             aes(x=x, y=y, colour = ClusterName),
             color = '#004949',
             size=3)+
  coord_fixed(ratio = .8) +theme(legend.position = "right")
dev.off()

highlight_Res <- subset(myeloid_fdl, ClusterName == "Res")

png(filename="./fdl Res ggplot2.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_fdl, aes(x = x, y = y, color = ClusterName)) +
  geom_point(size = 3, alpha = 1, color ="azure2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_Res,
             aes(x=x, y=y, colour = ClusterName),
             color = '#6db6ff',
             size=3)+
  coord_fixed(ratio = .8) +theme(legend.position = "right")
dev.off()

highlight_Arg1 <- subset(myeloid_fdl, ClusterName == "Arg1")

png(filename="./fdl Arg1 ggplot2.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_fdl, aes(x = x, y = y, color = ClusterName)) +
  geom_point(size = 3, alpha = 1, color ="azure2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_Arg1,
             aes(x=x, y=y, colour = ClusterName),
             color = '#006ddb',
             size=3)+
  coord_fixed(ratio = .8) +theme(legend.position = "right")
dev.off()

highlight_Prolif <- subset(myeloid_fdl, ClusterName == "Prolif")

png(filename="./fdl Prolif ggplot2.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_fdl, aes(x = x, y = y, color = ClusterName)) +
  geom_point(size = 3, alpha = 1, color ="azure2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_Prolif,
             aes(x=x, y=y, colour = ClusterName),
             color = '#490092',
             size=3)+
  coord_fixed(ratio = .8) +theme(legend.position = "right")
dev.off()

highlight_Mono <- subset(myeloid_fdl, ClusterName == "Mono")

png(filename="./fdl Mono ggplot2.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_fdl, aes(x = x, y = y, color = ClusterName)) +
  geom_point(size = 3, alpha = 1, color ="azure2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_Mono,
             aes(x=x, y=y, colour = ClusterName),
             color = '#db6d00',
             size=3)+
  coord_fixed(ratio = .8) +theme(legend.position = "right")
dev.off()

highlight_IFN <- subset(myeloid_fdl, ClusterName == "IFN")

png(filename="./fdl IFN ggplot2.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_fdl, aes(x = x, y = y, color = ClusterName)) +
  geom_point(size = 3, alpha = 1, color ="azure2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_IFN,
             aes(x=x, y=y, colour = ClusterName),
             color = '#ff6db6',
             size=3)+
  coord_fixed(ratio = .8) +theme(legend.position = "right")
dev.off()

highlight_cDC2 <- subset(myeloid_fdl, ClusterName == "cDC2")

png(filename="./fdl cDC2 ggplot2.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_fdl, aes(x = x, y = y, color = ClusterName)) +
  geom_point(size = 3, alpha = 1, color ="azure2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_cDC2,
             aes(x=x, y=y, colour = ClusterName),
             color = '#24ff24',
             size=3)+
  coord_fixed(ratio = .8) +theme(legend.position = "right")
dev.off()

#combined

png(filename = "./fdl predicted cell type ggplot2.png", width = 20, height = 20, units = "cm", res = 300)
ggplot(myeloid_fdl, aes(x = x , y = y, color = ClusterName))+
  geom_point(size = 2, alpha = 1)+
  scale_color_manual(values = c("#006ddb","#b66dff","#24ff24","#004949","#ff6db6","#009292","#db6d00","#490092","#6db6ff","#920000"
  ))+
  coord_fixed(ratio = .8) +theme(legend.position = "right")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

#cDC2 terminal
png(filename="./fdl summary cDC2 ggplot2.png", width=20, height=20, units="cm", res=300)
ggplot(myeloid_fdl, aes(x = x, y = y, color = cDC2)) +
  geom_point(size = 2, alpha = 1)+
  scale_color_gradientn(colors=c("#0D0887FF", "#7E03A8FF", "#CC4678FF", "#F89441FF", "#F0F921FF"))+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_fixed(ratio=.8)
dev.off()

#Trem2 terminal
png(filename="./fdl summary Trem2 ggplot2.png", width=20, height=20, units="cm", res=300)
ggplot(myeloid_fdl, aes(x = x, y = y, color = Trem2)) +
  geom_point(size = 2, alpha = 1)+
  scale_color_gradientn(colors=c("#0D0887FF", "#7E03A8FF", "#CC4678FF", "#F89441FF", "#F0F921FF"))+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_fixed(ratio=.8)
dev.off()

#monocytes, arg1, terminal populations: trem2 cdc2, on fdl
highlight_Mono <- subset(myeloid_fdl, ClusterName == "Mono")
highlight_Arg1 <- subset(myeloid_fdl, ClusterName == "Arg1")
highlight_Trem2 <- subset(myeloid_fdl, ClusterName == "Trem2")
highlight_cDC2 <- subset(myeloid_fdl, ClusterName == "cDC2")

png(filename="./fdl mono arg1 terminal.png", width=20, height=20, units="cm", res=600)
ggplot(myeloid_fdl, aes(x = x, y = y, color = ClusterName)) +
  geom_point(size = 4, alpha = 1, color ="azure2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = highlight_Mono,
             aes(x=x, y=y, colour = ClusterName),
             color = '#db6d00',
             size=5)+
  geom_point(data = highlight_Arg1,
             aes(x=x, y=y, colour = ClusterName),
             color = '#006ddb',
             size=5)+
  geom_point(data = highlight_Trem2,
             aes(x=x, y=y, colour = ClusterName),
             color = '#920000',
             size=5)+
  geom_point(data = highlight_cDC2,
             aes(x=x, y=y, colour = ClusterName),
             color = '#24ff24',
             size=5)+
  coord_fixed(ratio = .8) +theme(legend.position = "right")
dev.off()


