keptclusters<-(which(cell_clustering1==8|cell_clustering1==1|cell_clustering1==5|cell_clustering1==9|cell_clustering1==16|cell_clustering1==4|cell_clustering1==11))
keptclusters<-cbind(keptclusters,clusters=cell_clustering1[keptclusters])
tsne_ncells <- pmin(table(sample_ids), 932) 

set.seed(1234)
tsne_inds <- lapply(names(inds), function(i){
  s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
  intersect(s, dups)
})
tsne_inds <- unlist(tsne_inds)

kept<-intersect(tsne_inds,keptclusters[,1])

library(gplots)
ggdf <- data.frame(sample_id = sample_ids, expr)
ggdf<-ggdf[kept,]
centered_data <- data.frame((scale((data.matrix(ggdf[,2:20])), scale = T, center = T)))


clusterID<-data.frame()
for(i in 1:length(kept)){
  clusterID[i,1:2]<-keptclusters[which(keptclusters[,1]==kept[i]),]
}
colnames(clusterID)<-c("Index","cluster")



annotate_row<-as.data.frame(as.factor(clusterID$cluster))
rownames(annotate_row)<-clusterID$Index


#State<-as.factor(md$condition[match(ggdf$sample_id,md$sample_id)])
#annotate_row$State<-State
#colnames(annotate_row)<-c("Cluster","Disease")
colnames(annotate_row)<-c("Cluster")
annotate_row$Cluster<-factor(annotate_row$Cluster,levels=c(4,11,8,16,1,5,9))

sort.final<-centered_data[order(annotate_row$Cluster),]
#annotate_row<-annotate_row[order(annotate_row$Cluster),]

breaksList<-seq(-2, 2, by = 0.25)
annotation_colors<-c("#FDB462","#E6AB02","#33A02C","#7BAFDE","#E7298A","#882E72","#DC050C")
names(annotation_colors)<-unique(annotate_row)[,1]
my_color<-list(Cluster=annotation_colors)

pheatmap(t(sort.final),cluster_rows = TRUE, cluster_cols = FALSE,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList, annotation_col=annotate_row, labels_row = rownames(t(sort.final)),fontsize_row = 10,
         #annotation_col = annotate_col,
         annotation_colors = my_color,
         show_colnames = FALSE,show_rownames = TRUE)


#PIE CHART###############################################################################
lbls <- c("naïve", "SOD", "KO")
clustersTrack<-data.frame()
colbls<-paste("Cluster",unique(annotate_row$Cluster))

for (i in 1:length(unique(annotate_row$Cluster))){
  clustersTrack[1:3,i]<-c(length(annotate_row$Cluster[which(annotate_row$Cluster==unique(annotate_row$Cluster)[i]&annotate_row$Disease=="naïve")])/length(annotate_row$Cluster[which(annotate_row$Cluster==unique(annotate_row$Cluster)[i])]),
                          length(annotate_row$Cluster[which(annotate_row$Cluster==unique(annotate_row$Cluster)[i]&annotate_row$Disease=="SOD")])/length(annotate_row$Cluster[which(annotate_row$Cluster==unique(annotate_row$Cluster)[i])]),
                          length(annotate_row$Cluster[which(annotate_row$Cluster==unique(annotate_row$Cluster)[i]&annotate_row$Disease=="KO")])/length(annotate_row$Cluster[which(annotate_row$Cluster==unique(annotate_row$Cluster)[i])])
  )
}
rownames(clustersTrack)<-lbls
colnames(clustersTrack)<-colbls

par(mfrow=c(3,4) ) # 1 row and 3 columns for plots

pie(clustersTrack[,1], labels = paste0(lbls," ",round(clustersTrack[,1],digits=2)), xlab=colbls[1])
pie(clustersTrack[,2], labels = paste0(lbls," ",round(clustersTrack[,2],digits=2)), xlab=colbls[2])
pie(clustersTrack[,3], labels = paste0(lbls," ",round(clustersTrack[,3],digits=2)), xlab=colbls[3])
pie(clustersTrack[,4], labels = paste0(lbls," ",round(clustersTrack[,4],digits=2)), xlab=colbls[4])
pie(clustersTrack[,5], labels = paste0(lbls," ",round(clustersTrack[,5],digits=2)), xlab=colbls[5])
pie(clustersTrack[,6], labels = paste0(lbls," ",round(clustersTrack[,6],digits=2)), xlab=colbls[6])
pie(clustersTrack[,7], labels = paste0(lbls," ",round(clustersTrack[,7],digits=2)), xlab=colbls[7])
