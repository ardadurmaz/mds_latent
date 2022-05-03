library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(cluster)
library(randomForest)
library(ggplot2)

# Load Folds #
fold_cc <- sapply(1:5, simplify = FALSE, function(i){
  temp <- data.matrix(fread(sprintf('data/pca_cc_res_fold%d.tsv.gz', i), header = FALSE, sep='\t'))
  diag(temp) <- 1.0
  return(temp)
})

# Get Cluster Profiles #
sil_res <- sapply(fold_cc, simplify = FALSE, function(x){
  hc_res <- hclust(as.dist(1-x), method='ward.D2')
  local_sil_res <- sapply(2:20, function(k){
    temp <- summary(silhouette(cutree(hc_res, k), as.dist(1-x)))
    return(temp$avg.width)
  })
  return(local_sil_res)
})
plot_ft <- do.call('rbind', sapply(1:5, simplify = FALSE, function(i){
  return(data.table('NClust' = 2:20,
                    'Silhouette' = sil_res[[i]],
                    'Fold' = paste0('Fold-', i)))
}))

p <- ggplot(plot_ft, aes(x=NClust, y=Silhouette, color=Fold)) +
  xlab('# Clusters') +
  ylab('Average Silhouette Width') +
  geom_point() +
  geom_line() +
  theme_minimal() +
  theme(axis.title = element_text(face='bold', size=12))

ggsave(p, filename = 'plots/SilhouetteProfileFolds.pdf', width = 6, height = 4)


# RF-Model & Clusters #
n_cluster <- c(16, 14, 14, 14, 12)
rf_res <- sapply(1:5, simplify = FALSE, function(i){
  mut_mat <- data.matrix(fread(sprintf('data/mut_mat_v6_fold%d.tsv.gz', i), header=FALSE, sep='\t'))
  obs_ft <- fread(sprintf('data/mut_obs_v6_fold%d.tsv.gz', i), header=TRUE, sep='\t')
  feat_ft <- fread('data/mut_feat_v6.tsv.gz', header=TRUE, sep='\t')
  colnames(mut_mat) <- feat_ft$symbol  
  rownames(mut_mat) <- paste0('sample.', obs_ft$id)
  clust_res <- cutree(hclust(as.dist(1-fold_cc[[i]]), method = 'ward.D2'), k=n_cluster[i])
  obs_ft$cluster <- paste0('Cluster-', clust_res)
  local_rf <- randomForest(x=mut_mat, y=factor(obs_ft$cluster), mtry = 6, ntree = 500, nodesize=15)
  return(local_rf)
})

# Compare Clusters #
feat_ft <- fread('data/mut_feat_v6.tsv.gz', header=TRUE, sep='\t')
ari_mat <- matrix(0, ncol=5, nrow=5)
for(i in 1:5){
  mut_mat_i <- data.matrix(fread(sprintf('data/mut_mat_v6_fold%d.tsv.gz', i), header=FALSE, sep='\t'))
  colnames(mut_mat_i) <- feat_ft$symbol
  pred_i <- predict(rf_res[[i]], newdata=mut_mat_i)
  for(j in 1:5){
    pred_j <- predict(rf_res[[j]], newdata=mut_mat_i)
    ari_mat[i,j] <- mclust::adjustedRandIndex(pred_i, pred_j)
  }
}

# Plot #
cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.2f", ari_mat[i, j]), x, y, gp = gpar(fontsize = 10))
}
col_fun <- colorRamp2(breaks=seq(0,1,0.01), colors=viridis::viridis(length(seq(0,1,0.01))))
h <- Heatmap(ari_mat,
             col = col_fun,
             cell_fun = cell_fun,
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             show_column_names = TRUE,
             show_row_names = TRUE,
             row_labels = paste0('Fold-', 1:5),
             column_labels = paste0('Fold-', 1:5),
             row_names_gp = gpar(fontsize=12, fontface='bold'),
             column_names_gp = gpar(fontsize=12, fontface='bold'),
             na_col = 'white',
             heatmap_legend_param = list(title='Adjusted Rand Index', 
                                         title_position='leftcenter-rot', 
                                         legend_height=unit(4,'cm')),
             show_heatmap_legend = TRUE,
             use_raster=FALSE)
pdf('plots/ARI_Mat.pdf', width = 4, height = 3)
draw(h)
dev.off()