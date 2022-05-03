library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

setwd('~/Research/mds_latent/')

mut_mat <- data.matrix(fread('data/mut_mat_v6_train_imp.tsv.gz', header=FALSE, sep='\t'))
obs_ft <- fread('data/mut_obs_v6_train_annotated.tsv.gz', header=TRUE, sep='\t')
feat_ft <- fread('data/mut_feat_v6.tsv.gz', header=TRUE, sep='\t')
colnames(mut_mat) <- feat_ft$symbol
rownames(mut_mat) <- obs_ft$id

# Plot CC #
cc_mat <- data.matrix(fread('data/cc_res.tsv.gz', header=FALSE, sep='\t'))
hc.res <- hclust(as.dist(1-(cc_mat/100)), method = 'ward.D2')

# Heatmap 
man_colors <- obs_ft %>%
  distinct(cluster, color)

colAnnot <- HeatmapAnnotation(Cluster=obs_ft$cluster,
                              Diagnosis=obs_ft$diagnosis,
                              col = list(Cluster=setNames(man_colors$color, man_colors$cluster),
                                         Diagnosis=setNames(RColorBrewer::brewer.pal(3, 'Paired')[1:length(unique(obs_ft$diagnosis))], unique(obs_ft$diagnosis))),
                              which = 'row')

col_fun <- colorRamp2(breaks=seq(0,1,0.01), colors=viridis::viridis(length(seq(0,1,0.01))))
h <- Heatmap(as.matrix(cc_mat/100),
             col = col_fun,
             cluster_rows = hc.res,
             cluster_columns = hc.res,
             show_column_names = FALSE,
             show_row_names = FALSE,
             left_annotation = colAnnot,
             na_col = 'white',
             heatmap_legend_param = list(title='Consensus Matrix', 
                                         title_position='leftcenter-rot', 
                                         legend_height=unit(4,'cm')),
             show_heatmap_legend = TRUE,
             use_raster=FALSE)

png(sprintf('plots/cc_heatmap.png'), width = 720, height = 600)
draw(h)
dev.off()


colAnnot <- HeatmapAnnotation(Diagnosis=obs_ft$diagnosis,
                              col = list(Diagnosis=setNames(RColorBrewer::brewer.pal(3, 'Paired')[1:length(unique(obs_ft$diagnosis))], unique(obs_ft$diagnosis))),
                              which = 'row')
col_fun <- colorRamp2(breaks=seq(0,1,0.01), colors=viridis::viridis(length(seq(0,1,0.01))))
h <- Heatmap(mut_mat,
             col = col_fun,
             show_column_names = TRUE,
             show_row_names = FALSE,
             left_annotation = colAnnot,
             show_column_dend = FALSE,
             row_split = data.table('Clusters' = obs_ft$cluster),
             na_col = 'white',
             row_title_gp = gpar(fontsize=12),
             row_title_rot = 0,
             heatmap_legend_param = list(title='Imputed Matrix', 
                                         title_position='leftcenter-rot', 
                                         legend_height=unit(4,'cm')),
             show_heatmap_legend = TRUE,
             use_raster = FALSE)

png(sprintf('plots/feat_heatmap.png'), width = 580, height = 860)
draw(h)
dev.off()


# Frequency
freq_res <- sapply(split(1:nrow(mut_mat), obs_ft$cluster), simplify = FALSE, function(x){
  return(colSums(mut_mat[x,])/length(x))
})
freq_res <- do.call('rbind', freq_res)
freq_res <- freq_res[match(paste0('Cluster-', 0:13), table=rownames(freq_res)),]

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.3f", freq_res[i, j]), x, y, gp = gpar(fontsize = 8))
}

col_fun <- colorRamp2(breaks=seq(0,1,0.01), colors=colorRampPalette(c('white', 'firebrick'))(length(seq(0,1,0.01))))
h <- Heatmap(freq_res,
             col = col_fun,
             cell_fun = cell_fun,
             cluster_rows = FALSE,
             show_column_names = TRUE,
             show_row_names = TRUE,
             show_column_dend = FALSE,
             na_col = 'white',
             heatmap_legend_param = list(title='Mutation Frequency', 
                                         title_position='leftcenter-rot', 
                                         legend_height=unit(4,'cm')),
             show_heatmap_legend = TRUE,
             use_raster = FALSE)

pdf(sprintf('plots/freq_heatmap.pdf'), width = 16, height = 8)
draw(h)
dev.off()

# Survival
require(survival)
require(survminer)

man_colors <- obs_ft %>%
  distinct(color, cluster)

obs_ft$cluster <- factor(obs_ft$cluster, levels=paste0('Cluster-', 0:13))


fit <- survfit(Surv(os, event)~cluster, data=obs_ft)
p <- ggsurvplot(fit, conf.int = FALSE, pval = TRUE, risk.table = TRUE, palette = man_colors$color[match(paste0('Cluster-', 0:13), table=man_colors$cluster)])

comp_res <- matrix(0, ncol=14, nrow=14)
for(i in 1:13){
  for(j in (i+1):14){
    local_os <- subset(obs_ft, obs_ft$cluster %in% paste0('Cluster-', c(i-1, j-1)))
    fit <- survdiff(Surv(os, event) ~ cluster, data = local_os)
    comp_res[i,j] <- -1*pchisq(fit$chisq, 1, lower.tail = FALSE, log.p = TRUE)
  }
}

col_fun <- colorRamp2(breaks=seq(0,ceiling(max(comp_res)),0.001), colors=colorRampPalette(c('white', 'firebrick'))(length(seq(0,ceiling(max(comp_res)),0.001))))
h.ae <- Heatmap(comp_res,
                col = col_fun,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(j > i){
                    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "gray", fill = fill))
                    if(comp_res[i,j] > -log(0.001)){
                      grid.text(sprintf("%.1f***", comp_res[i, j]), x, y, gp = gpar(fontsize = 6))
                    }else if(comp_res[i,j] > -log(0.01)){
                      grid.text(sprintf("%.1f**", comp_res[i, j]), x, y, gp = gpar(fontsize = 6))
                    }else if(comp_res[i,j] > -log(0.05)){
                      grid.text(sprintf("%.1f*", comp_res[i, j]), x, y, gp = gpar(fontsize = 6))
                    }else(
                      grid.text(sprintf("%.1f", comp_res[i, j]), x, y, gp = gpar(fontsize = 6))
                    )
                  }else{
                    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "white", fill = 'white'))
                  }
                },
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_column_names = TRUE,
                show_row_names = TRUE,
                column_labels = paste0('Cluster-', 0:13),
                row_labels = paste0('Cluster-', 0:13),
                column_names_gp = gpar(fontface='bold', fontsize=8),
                row_names_gp = gpar(fontface='bold', fontsize=8),
                column_names_side = 'top',
                column_names_rot = 45,
                na_col = 'white',
                heatmap_legend_param = list(title='-log(P.Value)', 
                                            title_position='leftcenter-rot', 
                                            legend_height=unit(4,'cm')),
                show_heatmap_legend = TRUE,
                use_raster = FALSE)

pdf('plots/SurvivalHeatmap.pdf', width = 6, height = 4)
draw(h.ae)
dev.off()

# Cox-PH
obs_ft$sex[obs_ft$sex == ''] <- NA
obs_ft$bmFT <- log10(obs_ft$bm + 1.0)
cox_res <- coxph(Surv(os, event)~sex+age+cluster+cluster(cluster), data=obs_ft, x=TRUE, y=TRUE)
cox.zph(cox_res)
summary(cox_res)

local_ft <- data.table('sex'='M',
                       'age'=75,
                       'cluster'=paste0('Cluster-', 0:13))

plot(survfit(cox_res, newdata=local_ft),
     col = man_colors$color[match(levels(obs_ft$cluster), table=man_colors$cluster)],
     xlab = 'Time (Months)',
     ylab = 'S(t)',
     conf.int = TRUE,
     mark.time = FALSE,
     main = 'CoxPH-Baseline Survival (Age=75 & sex=M)',
     cumhaz=FALSE)

# RMS #
require(rms)

obs_ft$risk <- ifelse(obs_ft$cluster %in% paste0('Cluster-', c(3,9)), 'Low',
                      ifelse(obs_ft$cluster %in% paste0('Cluster-', c(1,5,7)), 'Int-Low',
                             ifelse(obs_ft$cluster %in% paste0('Cluster-', c(11,2,4,6,8,13)), 'Int-High',
                                    ifelse(obs_ft$cluster %in% paste0('Cluster-', c(0,10)), 'High', 'Poor'))))
dd <- datadist(obs_ft)
options(datadist='dd')

rms_cph_fit <- cph(Surv(os, event)~sex+age+risk, data=obs_ft, x=TRUE, y=TRUE, model=TRUE, se.fit=TRUE)
plot(summary(rms_cph_fit))

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

survplot(rms_cph_fit,
         age=75,
         sex = 'F',
         #cluster = paste0('Cluster-', c(3,10, 12)),
         risk = c('Low', 'Int-Low', 'Int-High', 'High', 'Poor'),
         col=RColorBrewer::brewer.pal(n=5, 'Dark2'),
         lty=2, 
         conf.int = TRUE,
         col.fill = add.alpha(RColorBrewer::brewer.pal(n=5, 'Dark2'), 0.3),
         xlab='Time (Months)', 
         n.risk = FALSE)


b <- bootcov(rms_cph_fit, B=1500, maxit=5000)
plot(summary(b))

# Clinical Features 
require(tableone)

tab1 <- CreateTableOne(vars = c('age', 'diagnosis', 'sex', 'wbc', 'hb', 'bm', 'anc', 'plt'), 
                       data = subset(obs_ft, !is.na(obs_ft$cluster)), 
                       strata = 'cluster', 
                       factorVars = c('sex', 'diagnosis'))

sink(sprintf('data/DescTable.tsv'))
print(tab1, nonnormal=c('wbc', 'hb', 'bm', 'anc', 'plt'))
sink()

# BM 
mean_levs <- obs_ft %>%
  group_by(cluster) %>%
  summarise(MedBM = median(bm, na.rm = TRUE))

obs_ft$cluster <- factor(obs_ft$cluster, levels=mean_levs$cluster[order(mean_levs$MedBM, decreasing = TRUE)])
p <- ggplot(obs_ft, aes(x=cluster, y=bm)) +
  xlab('Clusters') +
  ylab('% BM') +
  geom_boxplot(notch = TRUE) +
  theme_minimal() +
  theme(axis.title = element_text(face='bold', size=12))
ggsave(p, filename = 'plots/BM_Boxplot.pdf', width = 10, height = 6)

# Validation
require(randomForest)
rf_model <- readRDS('data/rf_model')

val_mat <- data.matrix(fread('data/mut_mat_v6_validation.tsv.gz', header=FALSE, sep='\t'))
val_obs_ft <- fread('data/mut_obs_v6_validation.tsv.gz', header=TRUE, sep='\t')
colnames(val_mat) <- feat_ft$symbol
rownames(val_mat) <- val_obs_ft$id
val_cluster <- predict(rf_model, newdata=val_mat)
val_obs_ft$cluster <- val_cluster

require(survival)
require(survminer)

fit <- survfit(Surv(os, event)~cluster, data=subset(val_obs_ft, val_obs_ft$cluster %in% paste0('Cluster-', c(0,1,12))))
p <- ggsurvplot(fit, conf.int = FALSE, pval = TRUE, risk.table = TRUE, palette = setNames(man_colors$color, paste0('cluster=', man_colors$cluster)))

# Circos Plot
man_colors <- obs_ft %>%
  distinct(cluster, color)
man_colors <- setNames(man_colors$color, man_colors$cluster)

for(c in 0:13){
  local_mut_mat <- mut_mat[obs_ft$cluster == paste0('Cluster-', c),]
  
  ## Extract Co-Occurrence ##
  freq <- (t(local_mut_mat) %*% local_mut_mat)/nrow(local_mut_mat)
  diag(freq) <- 0
  freq_ft <- reshape2::melt(as.matrix(freq))
  scale_ft <- function(x=NULL){
    return((x-min(x))/(max(x)-min(x)))
  }
  freq_ft$value <- scale_ft(exp(freq_ft$value))
  freq_ft <- freq_ft %>%
    filter(value > 0)
  freq_ft_cols <- sapply(freq_ft$value, function(alpha_val){
    temp <- col2rgb(man_colors[paste0('Cluster-', c)])[,1]
    return(rgb(temp[1]/255, temp[2]/255, temp[3]/255, alpha_val))
  })
  cyto_map <- setNames(c('normal',
                         'complex',
                         '-5/del(5q)', 
                         '-7/(del7q)', 
                         '-17/del(17p)', 
                         'del(20q)', 
                         '+8', 
                         '-Y', 
                         'other', 
                         expression(italic('ASXL1')),
                         expression(italic('CEBPA')),
                         expression(italic('CSF3R')),
                         expression(italic('CUX1')),
                         expression(italic('DDX41')),
                         expression(italic('DNMT3A')),
                         expression(italic('ETV6')),
                         expression(italic('EZH2')),
                         expression(italic('FLT3')),
                         expression(italic('GATA2')),
                         expression(italic('GNAS')),
                         expression(italic('JAK2')),
                         expression(italic('KIT')),
                         expression(italic('NOTCH1')),
                         expression(italic('NPM1')),
                         expression(italic('PHF6')),
                         expression(italic('RUNX1')),
                         expression(italic('SETBP1')),
                         expression(italic('SF3B1')),
                         expression(italic('SRSF2')),
                         expression(italic('TET2')),
                         expression(italic('TP53')),
                         expression(italic('U2AF1')),
                         expression(italic('WT1')),
                         expression(italic('ZRSR2')),
                         expression(italic('BCOR.BCORL1')),
                         expression(italic('IDH1.IDH2')),
                         expression(italic('NRAS.KRAS.PTPN11.NF1.CBL')),
                         expression(italic('RAD21.SMC1A.SMC3.STAG2')),
                         expression(italic('CALR.MPL'))),
                       colnames(mut_mat))
  
  pdf(sprintf('~/Research/mds_latent/plots/Cluster%d-Circos.pdf', c), width = 8, height = 8)
  circos.par("track.height" = 0.05, 
             "track.margin" = c(0.25, 0.25))
  circos.initialize(factors = colnames(mut_mat), xlim=c(0,1))
  circos.track(factors = colnames(mut_mat),
               ylim=c(0,1),
               bg.col = man_colors[paste0('Cluster-', c)],
               bg.border = 'black',
               panel.fun = function(x,y){
                 circos.text(CELL_META$xcenter, 
                             CELL_META$cell.ylim[2] + mm_y(1.75), ifelse(CELL_META$sector.index %in% names(cyto_map),
                                                                         cyto_map[which(names(cyto_map) == CELL_META$sector.index)],
                                                                         CELL_META$sector.index),
                             cex = 1,
                             facing = 'clockwise',
                             niceFacing = TRUE,
                             adj = c(0, 0.5))
               })
  for(i in 1:nrow(freq_ft)){
    circos.link(sector.index1 = freq_ft$Var1[i], c(0,1), sector.index2 = freq_ft$Var2[i], c(0,1), col = freq_ft_cols[i], rou = 0.675)  
  }
  dev.off()
  circos.clear()
}


