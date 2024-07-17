#Set working directory
setwd("/vol08/ngs/Seattle_Genomics/Collaborations/SG_Lifson/SG_Lifson_01/cornelius_analysis/de_work")

# Load libraries
library(NOISeq)
library(biomaRt)
source("./heatmap3LW_function.r")
library(edgeR)
library(tidyverse)
library(qusage)
library(corrplot)
library(factoextra)
library(ggfortify)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(svglite)

pca_fun <- function(exprs, labels, results_path,
                    base_file_name, target_columns,
                    figres = 100, size = 1, pca=FALSE, legend="right") {
    # Run PCA/SVD reduction
    if (isFALSE(pca)) {
        pca <- prcomp(t(exprs))
    }
    E <- get_eig(pca)
    cx <- sweep(t(exprs), 2, colMeans(t(exprs)), "-")
    sv <- svd(cx)


    vizualize_pca(
        file.path(results_path, paste0("svd_", base_file_name)),
        sv$u, labels[, target_columns[1]],
        labels[, target_columns[2]], figres, E, size, legend
    )
    vizualize_pca(
        file.path(results_path, paste0("pca_", base_file_name)),
        pca$x, labels[, target_columns[1]],
        labels[, target_columns[2]],
        figres, E, size, legend
    )
    vizualize_scree_plot(
        file.path(
            results_path,
            paste0("scree_", base_file_name)
        ), pca, figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    is_pc1_0 <- loadingscores$PC1 > 0
    is_pc2_0 <- loadingscores$PC2 > 0

    loadingscores <- loadingscores[is_pc1_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC1)), ]
    save_loading_scores(
        file.path(results_path, paste0("loadingscores_pc1", base_file_name, ".txt")),
        loadingscores["PC1"], figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[is_pc2_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC2)), ]
    save_loading_scores(
        file.path(results_path, paste0("loadingscores_pc2", base_file_name, ".txt")),
        loadingscores["PC2"], figres
    )
    return(pca)
}

vizualize_pca <- function(plot_file, PCA, class1, class2, figres, E, size, legend) {
    # Vizualize PCA  results
    library(Polychrome)
    minx <- min(PCA[, 1])
    maxx <- max(PCA[, 1])
    miny <- min(PCA[, 2])
    maxy <- max(PCA[, 2])
    # if (length(levels(factor(class2))) <= 3) {
    #     if (length(levels(factor(class1))) <= 6) {
    #         qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
    #             theme_Publication() +
    #             theme(legend.title = element_blank()) +
    #             scale_color_manual(values = c("Protected" = "red", "NonProtected" = "black", "Deceased" = "gray")) +
    #             scale_fill_manual(values =c("Protected" = "red", "NonProtected" = "black", "Deceased" = "gray")) +
    #             xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
    #             ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
    #             theme(legend.position = legend)
    #     } else {
    #         qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
    #             theme_Publication() +
    #             theme(legend.title = element_blank()) +
    #             scale_color_manual(values = c("IFN" = "pink", "CTL" = "black")) +
    #             scale_fill_manual(values = c("IFN" = "pink", "CTL" = "black")) +
    #             xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
    #             ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
    #             theme(legend.position = legend) +
    #             scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
    #     }
    # } else {
        P36 <- createPalette(length(levels(factor(class2))), c("#ff0000", "#00ff00", "#0000ff"))
        if (length(levels(factor(class1))) <= 6) {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = legend) +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36))
        } else if (length(levels(factor(class1))) > 6) {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = legend) +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    # }
    ggsave(plot_file, width = 6, height = 4, units = "in", dpi = 300)
}

theme_minimal_LW <- function(base_size = 14, base_family = "arial") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size, base_family = base_family)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.4, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(size = 10, face = "bold"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

theme_Publication <- function(base_size = 14, base_family = "arial") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size, base_family = base_family)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour = "#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.6, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face = "italic"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

scale_fill_Publication <- function(..) {
    library(scales)
    discrete_scale("fill", "Publication", manual_pal(values = c("#386cb0", "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33")), ..)
}

scale_colour_Publication <- function(..) {
    library(scales)
    discrete_scale("colour", "Publication", manual_pal(values = c("#386cb0", "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33")), ..)
}



# NOISeq analysis

# --Read in target files
message("STATUS: Load tables")
cm <- read.table("./count_matrix.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
target <- read.csv("./targetfile.csv", sep = ",", row.names = 1, as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)

# Rename sample IDs
newsampleIDs <- c()
for (i in colnames(cm)) {
    i <- str_remove(i, "_RNA\\d+_Lib\\d+\\S*$")
    i <- str_replace_all(i, "-", "_")
    newsampleIDs <- c(newsampleIDs, i)
}

colnames(cm) <- newsampleIDs

#Get necessary metadata
ensembl = useEnsembl(biomart="genes", dataset = "mmulatta_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
mmulata_BM <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "transcript_length", "percentage_gene_gc_content", "chromosome_name", "start_position", "end_position"), mart = ensembl)

# Remove duplicates and order the genes
all(row.names(cm) %in% mmulata_BM$ensembl_gene_id)
genomic_idx <- match(rownames(cm), mmulata_BM$ensembl_gene_id)
mmulata_BM_ordered <- mmulata_BM[genomic_idx,]
mmulata_BM_ordered <- na.omit(mmulata_BM_ordered)
cm <- cm[(rownames(cm) %in% mmulata_BM_ordered$ensembl_gene_id),]
all(row.names(cm) %in% mmulata_BM_ordered$ensembl_gene_id)

# Set up metadata files for downstream analysis
mmulata_BM_ordered <- data.frame(mmulata_BM_ordered, row.names = 1)
mylength <- mmulata_BM_ordered[,3, drop = FALSE]
mybiotypes <- mmulata_BM_ordered[,2, drop = FALSE]
mygc <- mmulata_BM_ordered[,4, drop = FALSE]
mychroms <- mmulata_BM_ordered[,5:7, drop = FALSE]

#Generate NOISeq object from raw counts

mydata <- readData(data = cm, length = mylength, gc = mygc, biotype = mybiotypes, chromosome = mychroms, factors = target[,c(1,5)])

# #QC of count data
mybiodetection <- dat(mydata, k = 0, type = "biodetection", factor = NULL)
# explo.plot(mybiodetection, samples = c(1,8), plottype = "persample")
# explo.plot(mybiodetection, samples = c(2,7), plottype = "persample")
# explo.plot(mybiodetection, samples = c(3,6), plottype = "persample")
# explo.plot(mybiodetection, samples = c(4,5), plottype = "persample")

mycountsbio = dat(mydata, factor = NULL, type= "countsbio")
# explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = 1, samples = 2, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = 1, samples = 3, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = 1, samples = 4, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = 1, samples = 5, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = 1, samples = 6, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = 1, samples = 7, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = 1, samples = 8, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "boxplot")
# explo.plot(mycountsbio,toplot=1,samples=NULL,plottype="barplot")

mysaturation = dat(mydata, k=0, ndepth = 7, type = "saturation")
# explo.plot(mysaturation, toplot = "protein_coding", samples = 1:8)

mylengthbias = dat(mydata, factor = "Time_Point", type="lengthbias")
# explo.plot(mylengthbias,samples=NULL,toplot="global")

mylengthbias = dat(mydata, factor = "Animal_ID", type="lengthbias")
# explo.plot(mylengthbias,samples=NULL,toplot="global")

myGCbias = dat(mydata, factor = "Time_Point", type = "GCbias")
# explo.plot(myGCbias, samples = NULL, toplot = "global")

myGCbias = dat(mydata, factor = "Animal_ID", type = "GCbias")
# explo.plot(myGCbias, samples = NULL, toplot = "global")

mycd = dat(mydata, type = "cd", norm =FALSE)
# explo.plot(mycd)

# PCA
myPCA = dat(mydata, type = "PCA")
par(mfrow = c(1,2))
explo.plot(myPCA, factor = "Animal_ID")
explo.plot(myPCA, factor = "Time_Point")

#Generate PDF report of the above
# QCreport(mydata, samples = NULL, factor = "Animal_ID", norm = FALSE)

# # Low-count filtering
myfilt = filtered.data(assayData(mydata)$exprs, norm = FALSE, depth = NULL, factor = target$Time_Point, method = 1, cv.cutoff = 100, cpm = 1, p.adj = "fdr")

# Make new NOISeq object from low count filtered expression data
mydata <- readData(data = myfilt, length = mylength, gc = mygc, biotype = mybiotypes, chromosome = mychroms, factors = target[,c(1,5)])

# Batch effect correction
mydatacorr = ARSyNseq(mydata, factor = "Animal_ID", batch = TRUE, logtransf = FALSE)
myPCA = dat(mydatacorr, type = "PCA")
par(mfrow = c(1,2))
explo.plot(myPCA, factor = "Animal_ID")
explo.plot(myPCA, factor = "Time_Point")

#Normalization
#setting long=1000 and lc=0 skips normalization for length
myTMM = tmm(assayData(mydatacorr)$exprs, long = 1000, lc = 0)

# Make new NOISeq object from normalized data
# Make new NOISeq object from low count filtered expression data
mydata <- readData(data = myTMM, length = mylength, gc = mygc, biotype = mybiotypes, chromosome = mychroms, factors = target[,c(1,5)])

#Re-run the above QC
# mybiodetection <- dat(mydata, k = 0, type = "biodetection", factor = NULL)
# explo.plot(mybiodetection, samples = c(1,8), plottype = "persample")
# explo.plot(mybiodetection, samples = c(2,7), plottype = "persample")
# explo.plot(mybiodetection, samples = c(3,6), plottype = "persample")
# explo.plot(mybiodetection, samples = c(4,5), plottype = "persample")

# mycountsbio = dat(mydata, factor = NULL, type= "countsbio")
# explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = 1, samples = 2, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = 1, samples = 3, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = 1, samples = 4, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = 1, samples = 5, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = 1, samples = 6, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = 1, samples = 7, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = 1, samples = 8, plottype = "boxplot")
# explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "boxplot")
# explo.plot(mycountsbio,toplot=1,samples=NULL,plottype="barplot")

# mysaturation = dat(mydata, k=0, ndepth = 7, type = "saturation")
# explo.plot(mysaturation, toplot = "protein_coding", samples = 1:8)

# mylengthbias = dat(mydata, factor = "Time_Point", type="lengthbias")
# explo.plot(mylengthbias,samples=NULL,toplot="global")

# mylengthbias = dat(mydata, factor = "Animal_ID", type="lengthbias")
# explo.plot(mylengthbias,samples=NULL,toplot="global")

# myGCbias = dat(mydata, factor = "Time_Point", type = "GCbias")
# explo.plot(myGCbias, samples = NULL, toplot = "global")

# myGCbias = dat(mydata, factor = "Animal_ID", type = "GCbias")
# explo.plot(myGCbias, samples = NULL, toplot = "global")

# mycd = dat(mydata, type = "cd", norm =FALSE)
# explo.plot(mycd)

# # PCA
# myPCA = dat(mydata, type = "PCA")
# par(mfrow = c(1,2))
# explo.plot(myPCA, factor = "Animal_ID")
# explo.plot(myPCA, factor = "Time_Point")

# #Generate PDF report of the above
# QCreport(mydata, samples = NULL, factor = "Animal_ID", norm = FALSE)


# Perform DE lets gooooooo

DE.plot(mynoiseqbio,q=0.95,graphic="expr",log.scale=TRUE)
DE.plot(mynoiseqbio,q=0.95,graphic="MD")

# Looked at simplifying pipeline. This sucked, don't try it.




mynoiseqbio7pp = noiseqbio(mydata, factor = "Time_Point", conditions = c("Day 7 post-prime", "Day 0 post-prime"), random.seed = 12345, norm = "n")
mynoiseqbio13pp = noiseqbio(mydata, factor = "Time_Point", conditions = c("Day 13 post-prime", "Day 0 post-prime"), random.seed = 12345, norm = "n")
mynoiseqbio11pb = noiseqbio(mydata, factor = "Time_Point", conditions = c("Day 11 post-boost", "Day 0 post-prime"), random.seed = 12345, norm = "n")


mynoiseq.deg7pp=degenes(mynoiseqbio7pp,q=0.95,M=NULL)
mynoiseq.deg13pp=degenes(mynoiseqbio13pp,q=0.95,M=NULL)
mynoiseq.deg11pb=degenes(mynoiseqbio11pb,q=0.95,M=NULL)

write.csv(mynoiseq.deg7pp, "./mynoiseq.deg7pp.csv")
write.csv(mynoiseq.deg13pp, "./mynoiseq.deg13pp.csv")
write.csv(mynoiseq.deg11pb, "./mynoiseq.deg11pb.csv")

# Generate more figuresssss
DE.plot(mynoiseqbio7pp,q=0.95,graphic="expr",log.scale=TRUE)
DE.plot(mynoiseqbio7pp,q=0.95,graphic="MD")

DE.plot(mynoiseqbio13pp,q=0.95,graphic="expr",log.scale=TRUE)
DE.plot(mynoiseqbio13pp,q=0.95,graphic="MD")

DE.plot(mynoiseqbio11pb,q=0.95,graphic="expr",log.scale=TRUE)
DE.plot(mynoiseqbio11pb,q=0.95,graphic="MD")

# Make a heatmap yo - not gonna work with the current data
asdf <- heatmap.L.4(mynoiseq.deg[,5, drop = FALSE],
    figmargins = c(20, 5),
    cutoff = 1, distmethod = "euclidean", cexcol = 2,
    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()

# Get gene lists for GSEA/ORA

D7pp_match <- match(rownames(mynoiseq.deg7pp), mmulata_BM_ordered$ensembl_gene_id)
D7pp_genes <- mmulata_BM_ordered[D7pp_match,]
D7pp_genes <- D7pp_genes[!(is.na(D7pp_genes$external_gene_name) | D7pp_genes$external_gene_name==""), ]

write.csv(D7pp_genes, "./D7pp_genes.csv")

D13pp_match <- match(rownames(mynoiseq.deg13pp), mmulata_BM_ordered$ensembl_gene_id)
D13pp_genes <- mmulata_BM_ordered[D13pp_match,]
D13pp_genes <- D13pp_genes[!(is.na(D13pp_genes$external_gene_name) | D13pp_genes$external_gene_name==""), ]

write.csv(D13pp_genes, "./D13pp_genes.csv")

D11pb_match <- match(rownames(mynoiseq.deg11pb), mmulata_BM_ordered$ensembl_gene_id)
D11pb_genes <- mmulata_BM_ordered[D11pb_match,]
D11pb_genes <- D11pb_genes[!(is.na(D11pb_genes$external_gene_name) | D11pb_genes$external_gene_name==""), ]

write.csv(D11pb_genes, "./D11pb_genes.csv")



####################################################
###### Trying this using EdgeR #####################
####################################################

# --Read in target files
message("STATUS: Load tables")
cm <- read.table("./count_matrix.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
target <- read.csv("./targetfile.csv", sep = ",", row.names = 1, as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)

# Rename sample IDs
newsampleIDs <- c()
for (i in colnames(cm)) {
    i <- str_replace_all(i, "-", "_")
    i <- str_replace_all(i, "_B_", "_")
    i <- str_replace_all(i, "3pb", "3pp")
    i <- str_remove(i, "^G..\\d+_Lifson01_AVP_081_")
    i <- str_remove(i, "_RNA\\d+_Lib\\d+\\S*$")
    
    newsampleIDs <- c(newsampleIDs, i)
}

colnames(cm) <- newsampleIDs

#Get necessary metadata
ensembl = useEnsembl(biomart="genes", dataset = "mmulatta_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
mmulata_BM <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "transcript_length", "percentage_gene_gc_content", "chromosome_name", "start_position", "end_position"), mart = ensembl)

# Remove duplicates and order the genes
all(row.names(cm) %in% mmulata_BM$ensembl_gene_id)
genomic_idx <- match(rownames(cm), mmulata_BM$ensembl_gene_id)
mmulata_BM_ordered <- mmulata_BM[genomic_idx,]
mmulata_BM_ordered <- na.omit(mmulata_BM_ordered)
cm <- cm[(rownames(cm) %in% mmulata_BM_ordered$ensembl_gene_id),]
all(row.names(cm) %in% mmulata_BM_ordered$ensembl_gene_id)

cm2 <- cm
cm2$external_gene_name <- mmulata_BM_ordered$external_gene_name


# asdf <- c("FA5T_D0pp", "FA5T_D7pp", "FA5T_D13pp", "FA5T_D11pb", "BH48_D11pb", 
#     "BH48_D13pp", "BH48_D7pp", "BH48_D0pp"
# )

asdf <- c("D0pp", "D7pp", "D13pp", "D11pb", "D11pb", 
    "D13pp", "D7pp", "D0pp"
)

# Account for experimental factors:
animal <- factor(substring(colnames(cm),1,4))
time <- factor(substring(colnames(cm), 6))



# Make DGE list object

d <- DGEList(counts = cm, group = time)

# Filter low read counts

dim(d)

keep <- filterByExpr(d)
d <- d[keep, , keep.lib.sizes=FALSE]

# Normalization (TMM)
d <- calcNormFactors(object = d)

plotMDS(d, col=as.numeric(d$samples$group), gene.selection="common")
legend("right", as.character(unique(d$samples$group)), col=1:4, pch=20)


# Make design matrix
design <- model.matrix(~animal + time, data=d$samples)
# design <- model.matrix(~group, data=d$samples)
design

# Reorder the matrix
d2 <- d
d2$samples$group <- relevel(d2$samples$group,ref="D11pb")

d2$samples$group <- relevel(d2$samples$group,ref="D13pp")

d2$samples$group <- relevel(d2$samples$group,ref="D7pp")

d2$samples$group <- relevel(d2$samples$group,ref="D0pp")

# Estimate common dispersion
d2 <- estimateDisp(d2, design)

plotMDS(d2, col=as.numeric(d2$samples$group))
plotBCV(d2)

#Perform DE test
fit <- glmQLFit(d2, design)
qlf<-glmQLFTest(fit, coef=c(5,4,3)) 
sig_stuff <- topTags(qlf, adjust.method = "none", p.value = 0.05, n = Inf)

sig_stuff2 <- sig_stuff$table

genomic_sig <- match(rownames(sig_stuff2), mmulata_BM$ensembl_gene_id)
mmulata_sig <- mmulata_BM[genomic_sig,]
sig_stuff2$external_gene_name <- mmulata_sig$external_gene_name

sig_stuff3 <- sig_stuff2
sig_stuff3$ensembl_id <- row.names(sig_stuff3)

write_csv(sig_stuff3, "sig_table.csv")

is.de <- decideTestsDGE(qlf, adjust.method = "none", p.value=0.05)
summary(is.de)

#get all DE genes
all_genes_DE <- topTags(qlf, n = Inf)
write.table(all_genes_DE, file = "DE_genes_all.txt")
table <- read.table("DE_genes_all.txt", sep = " ", header = TRUE)
write.csv(table, "./DE_genes_all.csv")

DE_all <- read.csv("./DE_genes_all.csv", sep = ",", row.names = 1, as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)
DE_all_match <- match(rownames(DE_all), mmulata_BM_ordered$ensembl_gene_id)
DE_all_names <- mmulata_BM_ordered[DE_all_match,]
#DE_all_names <- DE_all_names[!(is.na(DE_all_names$external_gene_name) | DE_all_names$external_gene_name==""),]
#DE_all_names <- as.matrix(DE_all_names)
DE_all <- cbind(DE_all, DE_all_names[,2])

colnames(DE_all)[colnames(DE_all) == "DE_all_names[, 1]"] = "external_gene_id"

# Volcano plots
# prep data
DE_all2 <- DE_all

DE_all2$logPValue <- lapply(DE_all2$PValue, function(x)-log10(x))
DE_all2$logPValue <- as.numeric(DE_all2$logPValue)

EnhancedVolcano(sig_stuff2,
                lab = sig_stuff2$external_gene_name,
                x = 'logFC.timeD7pp',
                y = 'PValue',
                title = 'Day 7 Post Prime',
                pCutoff = 0.05,
                legendLabels=c('','','p-value < 0.05',
                               'p-value < 0.05 & Log2 FC > 1'),
                legendPosition = 'top',
                drawConnectors = TRUE,
                max.overlaps = 17)

EnhancedVolcano(sig_stuff2,
                lab = sig_stuff2$external_gene_name,
                x = 'logFC.timeD13pp',
                y = 'PValue',
                title = 'Day 13 Post Prime',
                pCutoff = 0.05,
                drawConnectors = TRUE,
                max.overlaps = 39)

EnhancedVolcano(sig_stuff2,
                lab = sig_stuff2$external_gene_name,
                x = 'logFC.timeD11pb',
                y = 'PValue',
                title = 'Day 11 Post Boost',
                pCutoff = 0.05,
                drawConnectors = TRUE,
                max.overlaps = 20)


p1 <- ggplot(DE_all, aes(logFC.timeD7pp, -log(FDR, 10))) +
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))
p1


# Heatmap + stats

colcolorlistgrodwn <- c(rep("white"), rep("white"), rep("white"))
colcolormatrix <- as.matrix(colcolorlistgrodwn)

svglite("heatmap_final_transcriptome.svg", width = 10, height = 10)
asdf <- heatmap.L.4(sig_stuff2,
    figmargins = c(10, 5),
    cutoff = 1, distmethod = "spearman", cexcol = 2, colcolorlist = colcolormatrix,
    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()


# Based on the output from these, going with spearman clustering.
# It's also more robust to outliers which could be an issue with this dataset

# Get gene names for external ORA

for (cluster in unique(asdf$modulesrows)) {
    print(paste0("saving genes for ", cluster))
    genes <- which(asdf$modulesrows == cluster)
    gene_match <- match(names(genes), rownames(mmulata_BM_ordered))
    gene_names <- mmulata_BM_ordered[gene_match,]
    gene_names <- gene_names[!(is.na(gene_names$external_gene_name) | gene_names$external_gene_name==""),]
    write.csv(genes, paste0("cluster_",cluster,"_genes", ".csv"))
    write.csv(gene_names, paste0("cluster_",cluster,"_gene_names.csv"))
}

# Generate matrix of il15 DE genes that are significant

black <- read.csv("./gProfiler_il15_7-24-2023_black.csv", sep = ",", as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)
red <- read.csv("./gProfiler_il15_7-24-2023_red.csv", sep = ",", as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)
orange <- read.csv("./gProfiler_il15_7-24-2023_orange.csv", sep = ",", as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)
yellow <- read.csv("./gProfiler_il15_7-24-2023_yellow.csv", sep = ",", as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)

gene_list <- black[,10]
gene_list <- append(asdf, red[,10])
gene_list <- append(asdf, orange[,10])
gene_list <- append(asdf, yellow[,10])


gene_list2 <- str_split(asdf[1:6], ",")
gene_list2 <- unlist(qwerty)
gene_list2 <- unique(qwerty)

DE_all <- read.csv("./DE_genes_all.csv", sep = ",", row.names = 1, as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)
DE_all_match <- match(rownames(DE_all), rownames(mmulata_BM_ordered))
DE_all_names <- mmulata_BM_ordered[DE_all_match,]
DE_all_names <- DE_all_names[!(is.na(DE_all_names$external_gene_name) | DE_all_names$external_gene_name==""),]

qwerty2 <- DE_all_names[DE_all_names$external_gene_name %in% gene_list2,]
qwerty_match <- subset(DE_all, rownames(DE_all) %in% rownames(qwerty2) )

qwerty4 <- cbind(qwerty_match, qwerty2$external_gene_name)

colnames(qwerty4)[colnames(qwerty4) == "qwerty2$external_gene_name"] = "external_gene_name"

write.csv(qwerty4, "il15_de_genes.csv")

il15 <- read.csv("./il15_de_genes.csv", sep = ",", as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)

# This list was too big, need to match to the genes in the .gmt file that Leanne gave me

gmt <- read.gmt("./IL15_pathways.gmt")

head(gmt)
il15_down <- unlist(gmt$RM_IL15_DDE_OV_DWN_GENES)
il15_up <- unlist(gmt$RM_IL15_DDE_OV_UP_GENES)

il15_up2 <- il15[il15$external_gene_name %in% il15_up,]
il15_down2 <- il15[il15$external_gene_name %in% il15_down,]

write.csv(il15_up2, "il15_de_genes_up.csv")
write.csv(il15_down2, "il15_de_genes_down.csv")

il15_small <- rbind(il15_up2, il15_down2)
il15_small <- data.frame(il15_small, row.names = 1)
write.csv(il15_small, "il15_de_genes_small.csv")

# Generate heatmaps of the above
il15_small_hmap <- il15_small[,1:3]
il15_small_hmap <- data.frame(il15_small_hmap)
il15_small_hmap <- as.matrix(il15_small_hmap)

asdf <- heatmap.L.4(il15_small_hmap,
    figmargins = c(20, 5),
    cutoff = 1, distmethod = "pearson", cexcol = 2,
    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()

asdf <- heatmap.L.4(as.matrix(il15_small[,1:3]),
                    figmargins = c(20, 5),
                    cutoff = 1, distmethod = "pearson", cexcol = 2,
                    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()


asdf <- heatmap.L.4(il15_small_hmap,
    figmargins = c(20, 5),
    cutoff = 1, distmethod = "spearman", cexcol = 2,
    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()

asdf <- heatmap.L.4(il15_small_hmap,
    figmargins = c(20, 5),
    cutoff = 1, distmethod = "euclidean", cexcol = 2,
    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()

#Compare to il15 signature from Freddy's paper

DDE_IL15_overlap <- read.csv("/vol08/ngs/Seattle_Genomics/Collaborations/SG_Lifson/SG_Lifson_01/cornelius_analysis/de_work/DDE-IL15-overlap.csv", sep = ",", row.names = 1, as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)

DDE_IL15_overlap <- DDE_IL15_overlap[,c(5,9,13,17, 21, 25, 29, 33)]

il15_small_hmap_time_match <- il15_small[,c(1,3,8)]

DDE_IL15_overlap_small <- DDE_IL15_overlap[rownames(DDE_IL15_overlap) %in% rownames(il15_small_hmap_time_match),]

DDE_IL15_overlap_small <- as.matrix(DDE_IL15_overlap_small)

asdf <- heatmap.L.4(DDE_IL15_overlap_small,
    figmargins = c(20, 5),
    cutoff = 1, distmethod = "spearman", cexcol = 2,
    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()

il15_cccccombo <- cbind(il15_small_hmap_time_match, DDE_IL15_overlap_small)

il15_cccccombo <- data.frame(il15_cccccombo, row.names = il15_cccccombo[,3])
il15_cccccombo <- il15_cccccombo[,-3]

# Reorder columns for heatmap

il15_cccccombo <- il15_cccccombo[,c(1,2,7,8,3,4,9,10,5,6)]

colcolorlistgrodwn <- c(rep("white"), rep("white"), rep("white"), rep("white"), rep("white"), rep("white"), rep("white"), rep("white"), rep("white"), rep("white"))
colcolormatrix <- as.matrix(colcolorlistgrodwn)

pdf("il15combo_hmap.pdf", width = 15, height = 10)
asdf <- heatmap.L.4(as.matrix(il15_cccccombo),
    figmargins = c(20, 10),
    cutoff = 1, distmethod = "pearson", cexcol = 2, cexrow = 2, colcolorlist = colcolormatrix, labRow = TRUE,
    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()

for (cluster in unique(asdf$modulesrows)) {
  print(paste0("saving genes for ", cluster))
  genes <- which(asdf$modulesrows == cluster)
  gene_match <- match(names(genes), mmulata_BM_ordered$external_gene_name)
  gene_names <- mmulata_BM_ordered[gene_match,]
  gene_names <- gene_names[!(is.na(gene_names$external_gene_name) | gene_names$external_gene_name==""),]
  write.csv(genes, paste0("cluster_",cluster,"_genes_IL15.csv"))
  write.csv(gene_names, paste0("cluster_",cluster,"_gene_names_IL15.csv"))
}

il15_cor <- cor(il15_cccccombo)

corrplot(il15_cor, type = 'lower', col = colorRampPalette(c("blue", "white", "red" ))(100), diag = FALSE)

#Generate PCAs of heatmap data - This is for DDE data only

results <- prcomp(t(il15_cccccombo), scale = TRUE)


plot(results$x[,1], results$x[,2])
results.var <- results$sdev^2
results.var.per <- round(results.var/sum(results.var)*100, 1)
barplot(results.var.per, main= "Scree Plot", xlab = "Principal Component", ylab="Percent Variation")

results.data <- data.frame(Sample = rownames(results$x),
X=results$x[,1],
Y=results$x[,2])

ggplot(data=results.data, aes(x=X, y=Y, label = Sample)) +
geom_point(aes(color = Sample)) +
geom_text_repel()+
xlab(paste("PC1 - ", results.var.per[1],"%", sep = "")) +
ylab(paste("PC2 - ", results.var.per[2], "%", sep = "")) +
theme_bw() +
ggtitle("PCA of Lifson and Matched Barrenas Time Points, DDE Genes")


# Generate PCAs for full transcriptome
freddy_counts <- read.csv("/vol08/ngs/Seattle_Genomics/Collaborations/SG_Lifson/SG_Lifson_01/cornelius_analysis/GSE160562_O01_CountMatrix.csv", sep = ";", row.names = 1, as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)
freddy_counts_O_S <- freddy_counts[,-grep("^X", colnames(freddy_counts))]

target_freddy_o <- colnames(freddy_counts_O_S)

target_freddy <- read.table(text = target_freddy_o, sep = '_')

colnames(target_freddy) <- c("admin", "Animal_ID", "Time_Point", "cond")

rownames(target_freddy) <- target_freddy_o

target_freddy <- tidyr::separate(target_freddy, admin, c("admin", "num"), 1)

target_freddy_new <- target_freddy %>% rownames_to_column('row') %>% group_by(admin) %>% group_by(cond) %>% filter(Time_Point == "W0D7" | Time_Point == "W18D7") %>% column_to_rownames('row')



pca <- prcomp(t(freddy_counts_O_S))
autoplot(pca, data = target_freddy, color = "Time_Point")

freddy_counts_O_S_short <- freddy_counts_O_S[rownames(freddy_counts_O_S) %in% rownames(il15_cccccombo),]
pca <- prcomp(t(freddy_counts_O_S_short))
autoplot(pca, data = target_freddy, color = "Time_Point")
autoplot(pca, data = target_freddy, color = "admin")

pca_cm <- prcomp(t(cm))
autoplot(pca_cm, data = target, color = "Animal_ID") +
geom_text_repel(aes(label=Animal_ID))

head(cm)
cm_short <- cm[rownames(cm) %in% rownames(il15_cccccombo),]
head(cm_short)

pca_cm_short <- prcomp(t(cm_short))
autoplot(pca_cm_short, data = target, color = "Animal_ID") +
geom_text_repel(aes(label=Animal_ID))



# 3rd crack at generating these PCA plots
# Get unified and filtered list of genes from both count matrixes

all(row.names(cm) %in% row.names(freddy_counts_O_S))
all(row.names(freddy_counts_O_S) %in% row.names(cm))
cm_match <- match(rownames(cm), rownames(freddy_counts_O_S))
freddy_counts_O_S_ordered <- freddy_counts_O_S[cm_match,]
# mmulata_BM_ordered <- na.omit(mmulata_BM_ordered)
cm <- cm[(rownames(cm) %in% rownames(freddy_counts_O_S_ordered)),]
all(row.names(cm) %in% row.names(freddy_counts_O_S_ordered))
head(cm)
head(freddy_counts_O_S_ordered)


target_freddy_o <- colnames(freddy_counts_O_S_ordered)

target_freddy <- read.table(text = target_freddy_o, sep = '_')

colnames(target_freddy) <- c("admin", "Animal_ID", "Time_Point", "cond")

rownames(target_freddy) <- target_freddy_o

target_freddy <- tidyr::separate(target_freddy, admin, c("admin", "num"), 1)

target_freddy_new <- target_freddy %>% rownames_to_column('row') %>% group_by(admin) %>% group_by(cond) %>% filter(Time_Point == "W0D7" | Time_Point == "W18D7") %>% column_to_rownames('row')

freddy_match <- match(rownames(target_freddy_new), colnames(freddy_counts_O_S_ordered))

freddy_counts_O_S_ordered_match <- freddy_counts_O_S_ordered[,freddy_match]


combo_cm <- cbind(cm, freddy_counts_O_S_ordered_match)

pca_combo <- prcomp(t(combo_cm))

plot(pca_combo$x[,1], pca_combo$x[,2])
pca_combo.var <- pca_combo$sdev^2
pca_combo.var.per <- round(pca_combo.var/sum(pca_combo.var)*100, 1)
barplot(pca_combo.var.per, main= "Scree Plot", xlab = "Principal Component", ylab="Percent Variation")

pca_combo.data <- data.frame(Sample = rownames(pca_combo$x),
X=pca_combo$x[,1],
Y=pca_combo$x[,2])

ggplot(data=pca_combo.data, aes(x=X, y=Y, label = Sample)) +
geom_point() +
geom_text_repel()+
xlab(paste("PC1 - ", pca_combo.var.per[1],"%", sep = "")) +
ylab(paste("PC2 - ", pca_combo.var.per[2], "%", sep = "")) +
theme_bw() +
ggtitle("PCA of Lifson and Matched Barrenas Time Points, Full Transcriptome")


# Gonna normalize the combined count matrixes and generate a new pca

# Make DGE list object

time2 <- factor(substring(colnames(freddy_counts_O_S_ordered_match),11,15))

time <- as.character(time)
time2 <- as.character(time2)

time3 <- paste0(c(time,time2))

d <- DGEList(counts = combo_cm, group = time3)





# Make DGE list object

d <- DGEList(counts = combo_cm, group = time3)

# Filter low read counts

dim(d)

keep <- filterByExpr(d)
d <- d[keep, , keep.lib.sizes=FALSE]

# Normalization (TMM)
d <- calcNormFactors(object = d)

plotMDS(d, col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:4, pch=20)

# It worked yaaaaay

# Making more. Ugh.
d_pca <- prcomp(t(d$counts), scale = TRUE)


plot(d_pca$x[,1], d_pca$x[,2])
d_pca.var <- d_pca$sdev^2
d_pca.var.per <- round(d_pca.var/sum(d_pca.var)*100, 1)
barplot(d_pca.var.per, main= "Scree Plot", xlab = "Principal Component", ylab="Percent Variation")

d_pca.data <- data.frame(Sample = rownames(d_pca$x),
X=d_pca$x[,1],
Y=d_pca$x[,2])

ggplot(data=d_pca.data, aes(x=X, y=Y, label = Sample)) +
geom_point() +
geom_text_repel()+
xlab(paste("PC1 - ", d_pca.var.per[1],"%", sep = "")) +
ylab(paste("PC2 - ", d_pca.var.per[2], "%", sep = "")) +
theme_bw() +
ggtitle("PCA of Normalized Lifson and Barrenas Time Points, Full Transcriptome")


head(d$counts)
nrow(d$counts)
ncol(d$counts)

# Working on adding colors/shapes to PCA plots weeeeeee
admin <- c("L","L","L","L","L","L","L","L")
cond <- c("L","L","L","L","L","L","L","L")
target2 <- cbind(target, admin, cond)
target_combo <- rbind(target2[,43:44], target_freddy_new[,c(1,5)]) 
d_pca.data2 <- cbind(d_pca.data, target_combo)

ggplot(data=d_pca.data2, aes(x=X, y=Y, label = Sample)) +
  geom_point(aes(color = admin, shape = cond), size =3, alpha = 0.8) +
  #geom_text_repel()+
  xlab(paste("PC1 - ", d_pca.var.per[1],"%", sep = "")) +
  ylab(paste("PC2 - ", d_pca.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("PCA of Normalized Lifson and Barrenas Time Points, Full Transcriptome")

pca_combo.data2 <- cbind(pca_combo.data, target_combo)

ggplot(data=pca_combo.data2, aes(x=X, y=Y, label = Sample)) +
  geom_point(aes(color = admin, shape = cond), size =3, alpha = 0.8) +
  # geom_text_repel()+
  xlab(paste("PC1 - ", pca_combo.var.per[1],"%", sep = "")) +
  ylab(paste("PC2 - ", pca_combo.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("PCA of Lifson and Matched Barrenas Time Points, Full Transcriptome")

# Doing the same thing but for the time matched dde samples

admin <- c("L","L","S","S","S","S","O","O","O","O")
cond <- c("L","L","P","P","U","U","P","P","U","U")

results.data2 <- cbind(results.data, admin, cond)

ggplot(data=results.data, aes(x=X, y=Y, label = Sample)) +
  geom_point(aes(color = admin, shape = cond), size =3, alpha = 0.8) +
  # geom_text_repel()+
  xlab(paste("PC1 - ", results.var.per[1],"%", sep = "")) +
  ylab(paste("PC2 - ", results.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("PCA of Lifson and Matched Barrenas Time Points, DDE Genes")


# Remaking heat maps and correlation plots using all genes from Freddy list, not just those that are sig in our data
DE_all_time_match <- DE_all[,c(1,3)]

matching_DDE_inverse <- rownames(DDE_IL15_overlap) %in% rownames(DE_all_time_match)
matching_DDE_inverse_missing <- !(rownames(DDE_IL15_overlap) %in% rownames(DE_all_time_match))
matching_DDE <- rownames(DE_all_time_match) %in% rownames(DDE_IL15_overlap)
nrow(DE_all_time_match[matching_DDE,])
DDE_IL15_overlap[matching_DDE_inverse_missing,] # missing gene is: ENSMMUG00000050313
nrow(DDE_IL15_overlap[matching_DDE_inverse,])

IL15_combo_all_DDE <- cbind(DE_all_time_match[matching_DDE,], DDE_IL15_overlap[matching_DDE_inverse,])

IL15_combo_all_DDE <- as.matrix(IL15_combo_all_DDE)

colcolorlistgrodwn <- c(rep("white"), rep("white"), rep("white"), rep("white"), rep("white"), rep("white"), rep("white"), rep("white"), rep("white"), rep("white"))
colcolormatrix <- as.matrix(colcolorlistgrodwn)

IL15_combo_all_DDE <- IL15_combo_all_DDE[,c(1,2,7,8,3,4,9,10,5,6)]

par(mar= c(1,1,1,1))
pdf("il15all_DDE_hmap.pdf", width = 10, height = 10)
IL15_combo_all_DDE_hmap <- heatmap.L.4(IL15_combo_all_DDE,
                    figmargins = c(20, 5),
                    cutoff = 1, distmethod = "spearman", cexcol = 2, colcolorlist = colcolormatrix,
                    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()
IL15_combo_all_DDE <- as.data.frame(IL15_combo_all_DDE)

for (cluster in unique(IL15_combo_all_DDE_hmap$modulesrows)) {
  print(paste0("saving genes for ", cluster))
  genes <- which(IL15_combo_all_DDE_hmap$modulesrows == cluster)
  gene_match <- match(names(genes), mmulata_BM_ordered$ensembl_gene_id)
  gene_names <- mmulata_BM_ordered[gene_match,]
  gene_names <- gene_names[!(is.na(gene_names$external_gene_name) | gene_names$external_gene_name==""),]
  write.csv(genes, paste0("cluster_",cluster,"_genes_IL15_combo", ".csv"))
  write.csv(gene_names, paste0("cluster_",cluster,"_gene_names_IL15_combo.csv"))
}

il15_cor_all_DDE <- cor(IL15_combo_all_DDE)

corrplot(il15_cor_all_DDE, type = 'lower', col = colorRampPalette(c("blue", "white", "red" ))(100), diag = FALSE)

DE_all_w_names <- cbind(DE_all_names$external_gene_name, DE_all)


DE_all_w_names[matching_DDE,]
write.csv(DE_all_w_names[matching_DDE,], "./DE_all_matching_DDE.csv")

# Generate csv of all DDE genes
sig_stuff4 <- DE_all[matching_DDE,]
write.csv(sig_stuff4, "sig_stuff_4_185_dde.csv")

write.csv(DE_all, "full_DE_table.csv")
