setwd("/Users/jleibowitz/")
library(readxl)
library(RColorBrewer)
library(pheatmap)

#reading in sample metadata and antibody panel information
md<-read_excel("metaData.xlsx")
md$condition<-factor(md$condition,levels=c("naïve","SOD","KO"))

color_conditions <- c("#6A3D9A", "#FF7F00","#DC050C")
names(color_conditions) <- levels(md$condition)

panel<-read_excel("panel.xlsx")

setwd("/Users/jleibowitz/fcs_files")

#load in fcs_files
file_name<-dir()
library(flowCore)
fcs_raw<-read.flowSet(md$file_name,transformation=FALSE,truncate_max_range=FALSE)

# Replace problematic characters 
panel$Antigen <- gsub("-", "_", panel$Antigen)

#correct markernames
marknames<-panel$Antigen
names(marknames) <- panel$name
markernames(fcs_raw)<-marknames

panel_fcs <- pData(parameters(fcs_raw[[1]]))
head(panel_fcs)

# Replace problematic characters 
panel_fcs$desc <- gsub("-", "_", panel_fcs$desc)

# Lineage markers
(lineage_markers <- panel$Antigen[panel$Lineage == 1])

# Functional markers
(functional_markers <- panel$Antigen[panel$Functional == 1])

# Spot checks
all(lineage_markers %in% panel_fcs$desc)
all(functional_markers %in% panel_fcs$desc)


## arcsinh transformation and column subsetting
fcs <- fsApply(fcs_raw, function(x, cofactor = 500){
  colnames(x) <- panel_fcs$desc
  expr <- exprs(x)
  expr <- asinh(expr[, c(lineage_markers, functional_markers)] / cofactor)
  exprs(x) <- expr
  x
})
fcs

#different transforms and comparisons########
#biexp  <- biexponentialTransform("myTransform",a = 0.5, b = 1, c = 0.5, d = 1)
#log<-logicleTransform("myTransform")

#after.1 <- transform(fcs_raw, transformList(colnames(fcs_raw), biexp))

#biexp  <- biexponentialTransform("myTransform",w=10)
#after.2 <- transform(GvHD, transformList('FSC-H', biexp))

#opar = par(mfcol=c(3, 1))
#plot(density(exprs(fcs_raw[[1]])[, 19]), main="Original")
#plot(density(exprs(after.1[[1]])[, 19]), main="Standard Transform")
#plot(density(exprs(fcs[[1]])[, 10]), main="ARCsin Transform")


#fcs<-transform(fcs_raw, transformList(colnames(fcs_raw), biexp))
#colnames(fcs) <- as.character(panel_fcs$desc)
#fcs<-fcs[,10:29]


## Extract expression
expr <- fsApply(fcs, exprs)
dim(expr)

library(matrixStats)
rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1


## Generate sample IDs corresponding to each cell in the `expr` matrix
sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))

library(ggplot2)
library(reshape2)

ggdf <- data.frame(sample_id = sample_ids, expr)
ggdf <- melt(ggdf, id.var = "sample_id", 
             value.name = "expression", variable.name = "antigen")
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = expression, color = condition, 
                 group = sample_id)) +
  geom_density() +
  facet_wrap(~ antigen, nrow = 4, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
  guides(color = guide_legend(ncol = 1)) +
  scale_color_manual(values = color_conditions)


# Get the median marker expression per sample
library(dplyr)

expr_median_sample_tbl <- data.frame(sample_id = sample_ids, expr) %>%
  group_by(sample_id) %>% 
  summarize_all(funs(median))

expr_median_sample <- t(expr_median_sample_tbl[, -1])
colnames(expr_median_sample) <- expr_median_sample_tbl$sample_id

library(limma)
mds <- plotMDS(expr_median_sample, plot = FALSE)

library(ggrepel)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y, 
                   sample_id = colnames(expr_median_sample))
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = MDS1, y = MDS2, color = condition)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = sample_id)) +
  theme_bw() +
  scale_color_manual(values = color_conditions) 


############# CLUSTERING
library(FlowSOM)

fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)
set.seed(1234)
som <- BuildSOM(fsom, colsToUse = lineage_markers)

## Metaclustering into 20 clusters with ConsensusClusterPlus
library(ConsensusClusterPlus)

codes <- som$map$codes
plot_outdir <- "consensus_plots"
nmc <- 20

mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100, 
                           pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png", 
                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average", 
                           distance = "euclidean", seed = 1234)

## Get cluster ids for each cell
code_clustering1 <- mc[[nmc]]$consensusClass
cell_clustering1 <- code_clustering1[som$map$mapping[,1]]

color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", 
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", 
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", 
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999", 
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000", 
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

plot_clustering_heatmap_wrapper <- function(expr, expr01, 
                                            cell_clustering, color_clusters, cluster_merging = NULL){
  
  # Calculate the median expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% 
    summarize_all(funs(median))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% 
    summarize_all(funs(median))
  
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  
  # This clustering is based on the markers that were used for the main clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  
  labels_row <- paste0(rownames(expr_heat), " (", 
                       round(clustering_table / sum(clustering_table) * 100, 2), "%)")
  labels_col <- colnames(expr_heat)
  
  # Row annotation for the heatmap
  annotation_row <- data.frame(cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  
  color_clusters <- color_clusters[1:nlevels(annotation_row$cluster)]
  names(color_clusters) <- levels(annotation_row$cluster)
  annotation_colors <- list(cluster = color_clusters)
  annotation_legend <- FALSE
  
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$cluster_merging <- cluster_merging$new_cluster 
    color_clusters <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters) <- levels(cluster_merging$new_cluster)
    annotation_colors$cluster_merging <- color_clusters
    annotation_legend <- TRUE
  }
  
  # Colors for the heatmap
  color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  
  pheatmap(expr_heat, color = color, 
           cluster_cols = TRUE, cluster_rows = cluster_rows, 
           labels_col = labels_col, labels_row = labels_row, 
           display_numbers = TRUE, number_color = "black", 
           fontsize = 8, fontsize_number = 4,
           annotation_row = annotation_row, annotation_colors = annotation_colors, 
           annotation_legend = annotation_legend)
  
}

plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers], 
                                expr01 = expr01[, lineage_markers], 
                                cell_clustering = cell_clustering1, color_clusters = color_clusters)

screeplot(prcomp(expr_median))

## Find and skip duplicates
dups <- which(!duplicated(expr[, lineage_markers]))

## Data subsampling: create indices by sample
inds <- split(1:length(sample_ids), sample_ids) 

## How many cells to downsample per-sample
tsne_ncells <- pmin(table(sample_ids), 2000)  

## Get subsampled indices
set.seed(1234)
tsne_inds <- lapply(names(inds), function(i){
  s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
  intersect(s, dups)
})

tsne_inds <- unlist(tsne_inds)

tsne_expr <- expr[tsne_inds, lineage_markers]

## Run t-SNE
library(Rtsne)

set.seed(1234)
tsne_out <- Rtsne(tsne_expr, check_duplicates = FALSE, pca = FALSE)

## Plot t-SNE colored by CD4 expression
dr <- data.frame(tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2], 
                 expr[tsne_inds, lineage_markers])

ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = Ly6C)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_gradientn("Ly6C", 
                        colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))

dr$sample_id <- sample_ids[tsne_inds]
mm <- match(dr$sample_id, md$sample_id)
dr$condition <- md$condition[mm]
dr$cell_clustering1 <- factor(cell_clustering1[tsne_inds], levels = 1:nmc)

#
dr$patient_id <- md$patient_id[mm]
#

## Plot t-SNE colored by clusters
ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = cell_clustering1)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))


## individual################################################################

gg<-ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = cell_clustering1)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))

pdf("tSNE_bycondition.pdf",height=7,width=12)
gg + facet_wrap(~ condition) 
dev.off()

## Plot t-SNE colored by  expression for each antigen

for (i in 10:nrow(panel)){
  gg<-ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = get(panel$Antigen[i]))) +
    geom_point(size = 0.8) +
    theme_bw() +
    scale_color_gradientn(as.character(panel$Antigen[i]), 
                          colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))
  
  pdf(paste0("tSNE_",panel$Antigen[i],"_2000.pdf"),height=7,width=12)
  
  print(gg + facet_wrap(~ condition)) 
  dev.off()
}


gg<-ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = Clec7a)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_gradientn("Clec7a", 
                        colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))

pdf(paste0("tSNE_4C12.pdf"),height=7,width=12)
print(gg + facet_wrap(~ condition))
dev.off()




############################################################################################################

#HEATMAP FOR ALL MARKERS# SEE HEATMAP_AURORA.R FOR CLUSTER SPECIFIC HEATMAP

## Get subsampled indices
tsne_ncells <- pmin(table(sample_ids), 932) 

set.seed(1234)
tsne_inds <- lapply(names(inds), function(i){
  s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
  intersect(s, dups)
})
tsne_inds <- unlist(tsne_inds)


library(gplots)
ggdf <- data.frame(sample_id = sample_ids, expr)
centered_data <- data.frame((scale((data.matrix(ggdf[,2:20])), scale = T, center = T)))

clusterID<-data.frame()
for(i in 1:length(kept)){
  clusterID[i,1:2]<-keptclusters[which(keptclusters[,1]==kept[i]),]
}
colnames(clusterID)<-c("Index","cluster")




temp<-ggdf[,2:20]
for(p in 1:ncol(temp)){temp[,p]<-as.numeric(temp[,p])}
pca<-(prcomp(t(temp), scale.=TRUE))
screeplot(pca)
# perform k-means
set.seed(12345)
centers <- 6
final <- data.frame(centered_data)

# Assign genes to a cluster and sort clusters by 1-10
km <- kmeans(final, centers = centers,algorithm="Lloyd",iter.max=100)
m.kmeans<- cbind(final, cluster = km$cluster, sample_id=ggdf$sample_id)
order <- m.kmeans[order(m.kmeans$cluster),]


# arrange clusters based on size (highest --> low )
len <- c()
for(i in 1:centers){
  len <- c(len, length(which(order$cluster == i)))
}

len2 <- order(-len)

sort.final <- data.frame()
for ( x in len2){
  sort.final <- rbind(sort.final, order[which(order$cluster == x),])
}

# assign colors to clusters -- can change 
colours <-  colorRampPalette(brewer.pal(11,"Spectral"))(centers)[sort.final$cluster]

# generate heatmap with clusters 
sort.final<-data.matrix(sort.final)



#incl_gene_names<-vector(mode="character")
#incl_gene_names<-c("^P2ry12$","^Apoe$","^Spp1$","^Gpnmb$","^Bhlhe40$","^Entpd1$","^Fcrls$","^Tmem119$","^Olfml3$","^Hexb$","^Irf8$","^Tgfbr1$",
#                   "^Jun$","^Ccl2$","^C1qa$","^Bin1$","^Psen2$") #INPUT GENE NAMES TO INCLUDE ON HEAT MAP
#rowlabels<-vector(mode="character", length=nrow(sort.final))

#for (i in 1:length(incl_gene_names)){
#  rowlabels[grep(incl_gene_names[i],rownames(sort.final))]<-incl_gene_names[i]
#}

sort.final<-sort.final[order(sort.final[,(ncol(sort.final)-1)]),]

library(RColorBrewer)

annotate_row<-as.data.frame(as.factor(sort.final[1:nrow(sort.final),(ncol(sort.final)-1)]))
rownames(annotate_row)<-rownames(sort.final)

State<-as.factor(md$condition[match(sort.final[,"sample_id"],md$sample_id)])
annotate_row$State<-State
colnames(annotate_row)<-c("Cluster","Disease")


sort.final2<-sort.final[,order(colnames(sort.final)[1:19])]

#annotate_col<-as.data.frame(c(rep("Unstimulated",24),rep("TGFb",24),rep("IFNg",24),rep("TGFb_IFNg",24)))
#annotate_col<-as.data.frame(c(rep("Female",48),rep("Male",48)))
#annotate_col<-as.data.frame(c(rep("Healthy",15),rep("Mouse",13),rep("Progressive",15),rep("Relapsing",15)))
#colnames(annotate_col)<-"Samples"
#rownames(annotate_col)<-colnames(sort.final2)
#Patient<-as.factor(c(rep("H1",5),rep("H2",5),rep("H3",5),rep("H4",5),rep("M1",4),rep("M2",5),rep("M3",4),rep("M4",4),
#                     rep("P1",4),rep("P2",5),rep("P3",5),rep("P4",5),rep("R1",5),rep("R2",5),rep("R3",5),rep("R4",5)))
#annotate_col$Patient<-Patient

breaksList<-seq(-2, 2, by = 0.25)

pheatmap(sort.final2,cluster_rows = FALSE, cluster_cols = TRUE,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList, annotation_row=annotate_row, labels_row = rownames(sort.final2),fontsize_row = 10,
         #annotation_col = annotate_col, 
         show_colnames = TRUE,show_rownames = FALSE)



# save the clusters and genes z-scores 
write.xlsx(cbind(sort.final2,sort.final[,ncol(sort.final)]), '3_clusters_genes_withoutexp1.xlsx',
           row.names = TRUE)

pheatmap(t(sort.final2),cluster_rows = TRUE, cluster_cols = FALSE,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList, annotation_col=annotate_row, labels_row = rownames(t(sort.final2)),fontsize_row = 10,
         #annotation_col = annotate_col, 
         show_colnames = FALSE,show_rownames = TRUE)

length(annotate_row$Cluster[which(annotate_row$Cluster==2&annotate_row$Disease=="naïve")])/length(annotate_row$Cluster[which(annotate_row$Cluster==2)])
length(annotate_row$Cluster[which(annotate_row$Cluster==2&annotate_row$Disease=="SOD")])/length(annotate_row$Cluster[which(annotate_row$Cluster==2)])
length(annotate_row$Cluster[which(annotate_row$Cluster==2&annotate_row$Disease=="KO")])/length(annotate_row$Cluster[which(annotate_row$Cluster==2)])


# Simple Pie Chart
lbls <- c("naïve", "SOD", "KO")
clustersTrack<-data.frame()

for (i in 1:length(unique(annotate_row$Cluster))){
  clustersTrack[1:3,i]<-c(length(annotate_row$Cluster[which(annotate_row$Cluster==i&annotate_row$Disease=="naïve")])/length(annotate_row$Cluster[which(annotate_row$Cluster==i)]),
                          length(annotate_row$Cluster[which(annotate_row$Cluster==i&annotate_row$Disease=="SOD")])/length(annotate_row$Cluster[which(annotate_row$Cluster==i)]),
                          length(annotate_row$Cluster[which(annotate_row$Cluster==i&annotate_row$Disease=="KO")])/length(annotate_row$Cluster[which(annotate_row$Cluster==i)])
  )
}

par(mfrow=c(1,6) ) # 1 row and 3 columns for plots

pie(clustersTrack[,1], labels = paste0(lbls," ",round(clustersTrack[,1],digits=2)), xlab="Cluster 1")
pie(clustersTrack[,2], labels = paste0(lbls," ",round(clustersTrack[,2],digits=2)), xlab="Cluster 2")
pie(clustersTrack[,3], labels = paste0(lbls," ",round(clustersTrack[,3],digits=2)), xlab="Cluster 3")
pie(clustersTrack[,4], labels = paste0(lbls," ",round(clustersTrack[,4],digits=2)), xlab="Cluster 4")
pie(clustersTrack[,5], labels = paste0(lbls," ",round(clustersTrack[,5],digits=2)), xlab="Cluster 5")
pie(clustersTrack[,6], labels = paste0(lbls," ",round(clustersTrack[,6],digits=2)), xlab="Cluster 6")


