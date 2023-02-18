library(data.table)
library(dplyr)
library(randomForest)
library(caret)

setwd('/media/ardadurmaz/HelperWin/Research/mds_latent/')

mut_mat <- data.matrix(fread('data/mut_mat_v6_train_imp.tsv.gz', header=FALSE, sep='\t'))
obs_ft <- fread('data/mut_obs_v6_train_annotated.tsv.gz', header=TRUE, sep='\t')
feat_ft <- fread('data/mut_feat_v6.tsv.gz', header=TRUE, sep='\t')
colnames(mut_mat) <- feat_ft$symbol
rownames(mut_mat) <- obs_ft$id

risk_m <- readRDS('app_data/rf_model_risk.rds')
pred_res <- predict(risk_m, newdata = mut_mat)
obs_ft$risk <- pred_res
fwrite(obs_ft, file='data/mut_obs_v6_train_annotated_risk.tsv.gz', col.names = TRUE, row.names = FALSE, sep = '\t', append = FALSE, quote = FALSE)

# CV
split_idx <- createDataPartition(y=obs_ft$cluster, times=10)
cv_res <- sapply(2:6, function(m){
  temp <- sapply(split_idx, function(test_idx){
    train_idx <- setdiff(1:nrow(obs_ft), test_idx)
    rf <- randomForest(x=mut_mat[train_idx,], 
                       y=factor(obs_ft$cluster[train_idx]), 
                       ntree=1500, 
                       importance=FALSE, 
                       nodesize=15, 
                       mtry=m)
    pred_res <- predict(rf, newdata = mut_mat[test_idx,])
    crosstab <- table(pred_res, obs_ft$cluster[test_idx])
    return(sum(diag(crosstab))/sum(crosstab))
  })
  return(temp)
})

pdf('plots/CV_Res.pdf')
boxplot(cv_res, main='Hyperparameter Selection', xlab='', ylab='Accuracy', xaxt='n')
axis(side=1, at=1:5, labels=paste0('M-', 2:6))
dev.off()

imp_res <- sapply(1:10, simplify = FALSE, function(i){
  rf_model <- randomForest(x=mut_mat, 
                           y=factor(obs_ft$cluster), 
                           ntree=1500, 
                           importance=TRUE, 
                           localImp = TRUE,
                           nodesize=15, 
                           mtry=3)
  agg_imp <- do.call('cbind', sapply(split(1:nrow(obs_ft), obs_ft$cluster), simplify = FALSE, function(x){
    rowMeans(rf_model$localImportance[,x])
  }))
  temp <- melt(agg_imp)
  colnames(temp) <- c('feat', 'clust', 'imp')
  return(temp)
  #return(rf_model$importance[,ncol(rf_model$importance)-1])
})
imp_res_ft <- do.call('rbind', imp_res)
fwrite(imp_res_ft, file='data/LocalImportance.tsv.gz', col.names = TRUE, row.names = FALSE, sep='\t', append=FALSE, quote=FALSE)
imp_res_ft_agg <- imp_res_ft %>%
  group_by(feat, clust) %>%
  summarise(impM = mean(imp))

imp_res_ft_agg$clust <- factor(imp_res_ft_agg$clust, levels=paste0('Cluster-', 0:13))
man_colors <- obs_ft %>%
  distinct(cluster, color)

p <- ggplot(imp_res_ft_agg, aes(x=feat, y=impM, fill=clust)) +
  ylab('Mean Decrease Accuracy') +
  xlab('') +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~clust, ncol=1) +
  scale_fill_manual(values = setNames(man_colors$color, man_colors$cluster)) +
  guides(fill="none")
ggsave(p, filename = 'plots/LocalImportance.pdf', width = 8, height = 16)

plot_ft <- reshape2::melt(imp_res)
plot_ft$Var1 <- factor(plot_ft$Var1, levels = names(sort(rowMeans(imp_res), decreasing = TRUE)))
p <- ggplot(plot_ft, aes(x=Var1, y=value)) +
  xlab('') +
  ylab('Mean Decrease in Accuracy') +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1, angle = 45, face='bold', size=8))
ggsave(p, filename = 'plots/Importance.pdf')

fwrite(t(imp_res), file = 'data/Importance.tsv.gz', col.names = TRUE, row.names = FALSE, sep='\t', append=FALSE, quote=FALSE)

# Predict Test Data
rf_model <- randomForest(x=mut_mat, 
                         y=factor(obs_ft$cluster), 
                         ntree=1500, 
                         importance=TRUE, 
                         nodesize=15, 
                         mtry=3)
saveRDS(rf_model, file='data/rf_model')
mut_mat_test <- data.matrix(fread('data/mut_mat_v6_test_imp.tsv.gz', header=FALSE, sep='\t'))
obs_ft_test <- fread('data/mut_obs_v6_test.tsv.gz', header=TRUE, sep='\t')
colnames(mut_mat_test) <- feat_ft$symbol
rownames(mut_mat_test) <- obs_ft_test$id


pred_res <- predict(rf_model, newdata = mut_mat_test)
obs_ft_test$cluster <- pred_res
fwrite(obs_ft_test, file='data/mut_obs_v6_test_annotated.tsv.gz', sep='\t', append=FALSE, quote=FALSE)

# Risk Groups
obs_ft$cluster <- factor(obs_ft$cluster, levels=paste0('Cluster-', 0:13))
obs_ft$risk <- ifelse(obs_ft$cluster %in% paste0('Cluster-', c(3,9)), 'Low',
                      ifelse(obs_ft$cluster %in% paste0('Cluster-', c(1,5,7)), 'Int-Low',
                             ifelse(obs_ft$cluster %in% paste0('Cluster-', c(11,2,4,6,8,13)), 'Int-High',
                                    ifelse(obs_ft$cluster %in% paste0('Cluster-', c(0,10)), 'High', 'Poor'))))

split_idx <- createDataPartition(y=obs_ft$cluster, times=5)
cv_res_risk <- sapply(2:6, function(m){
  temp <- sapply(split_idx, function(test_idx){
    train_idx <- setdiff(1:nrow(obs_ft), test_idx)
    rf <- randomForest(x=mut_mat[train_idx,], 
                       y=factor(obs_ft$risk[train_idx]), 
                       ntree=1500, 
                       importance=FALSE, 
                       nodesize=15, 
                       mtry=m)
    pred_res <- predict(rf, newdata = mut_mat[test_idx,])
    crosstab <- table(pred_res, obs_ft$risk[test_idx])
    return(sum(diag(crosstab))/sum(crosstab))
  })
  return(temp)
})

pdf('plots/CV_Res_Risk.pdf')
boxplot(cv_res_risk, main='Hyperparameter Selection', xlab='', ylab='Accuracy', xaxt='n')
axis(side=1, at=1:5, labels=paste0('M-', 2:6))
dev.off()

risk_rf_model <- randomForest(x=mut_mat, 
                         y=factor(obs_ft$risk), 
                         ntree=1500, 
                         importance=TRUE, 
                         nodesize=15, 
                         mtry=3)

saveRDS(risk_rf_model, file='data/rf_model_risk')
mut_mat_test <- data.matrix(fread('data/mut_mat_v6_test_imp.tsv.gz', header=FALSE, sep='\t'))
obs_ft_test <- fread('data/mut_obs_v6_test.tsv.gz', header=TRUE, sep='\t')
colnames(mut_mat_test) <- feat_ft$symbol
rownames(mut_mat_test) <- obs_ft_test$id


# Survival
require(survival)
require(survminer)
require(tidyverse)

p <- ggsurvplot(survfit(Surv(os, event)~cluster, data=obs_ft_test), conf.int = FALSE)

obs_ft_noise <- obs_ft[!(obs_ft$id %in% obs_ft_ann$id),]
mut_mat_noise <- mut_mat[match(obs_ft_noise$id, table=obs_ft$id),]
colnames(mut_mat_noise) <- feat_ft$symbol
rownames(mut_mat_noise) <- obs_ft_noise$id

rf_model <- readRDS('HR_RF_Model.rds')
pred_res <- predict(rf_model, newdata = mut_mat)
obs_ft$cluster <- pred_res
obs_ft$type <- ifelse(obs_ft$cluster %in% c('Cluster-6', 'Cluster-5'), 'Low Risk', 'High Risk')
#obs_ft$type <- ifelse(obs_ft$cluster %in% c('Cluster-1'), 'Low Risk', 'High Risk')


do_test <- function(x=NULL){
  temp <- fisher.test(matrix(c(x$Found[x$type == 'High Risk'], x$NotFound[x$type == 'High Risk'], x$Found[x$type == 'Low Risk'], x$NotFound[x$type == 'Low Risk']), 
                             nrow=2, 
                             byrow = TRUE), 
                      alternative = 'two.sided')
  return(data.frame('p.val'=temp$p.value))
} 

# Histogram #
pval <- reshape2::melt(mut_mat) %>%
  left_join(obs_ft, by = c('Var1'='id')) %>%
  group_by(Var2, type) %>%
  summarise(Freq = sum(value)/n(),
            Found = sum(value == 1),
            NotFound = sum(value == 0)) %>%
  do(do_test(x=.))

pval$featFT <- sapply(1:nrow(pval), function(i){
  if(pval$p.val[i] < 1e-3){
    return(paste0(pval$Var2[i], ' ***'))
  }else if(pval$p.val[i] < 1e-2){
    return(paste0(pval$Var2[i], ' **'))
  }else if(pval$p.val[i] < 1e-1){
    return(paste0(pval$Var2[i], ' *'))
  }else{
    return(as.character(pval$Var2[i]))
  }
})

freq <- reshape2::melt(mut_mat) %>%
  left_join(obs_ft, by = c('Var1'='id')) %>%
  group_by(Var2, type) %>%
  summarise(Freq = sum(value)/n())

freq <- freq %>%
  left_join(pval, by='Var2')
freq$featFT <- factor(freq$featFT, levels = pval$featFT[order(pval$p.val, decreasing = FALSE)])

require(ggplot2)
p <- ggplot(freq, aes(x=featFT, y=Freq, fill=type)) +
  geom_col(position = position_dodge()) +
  xlab('') +
  ylab('Frequency') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face='bold', size=12),
        legend.title = element_blank())

ggsave(p, filename = 'plots/LR_HR_Histogram.pdf', width = 10, height = 6)

require(survival)
require(survminer)

obs_ft$event[obs_ft$os > 120] <- 0
obs_ft$os[obs_ft$os > 120] <- 120

fit <- survfit(Surv(os, event)~type, data=obs_ft)
p <- ggsurvplot(fit, 
                data = obs_ft,
                conf.int = FALSE, 
                pval = TRUE, 
                risk.table = TRUE)


pred_res <- predict(rf_model, newdata = mut_mat_noise)
obs_ft_noise$Cluster <- pred_res
obs_ft_merged <- rbind(obs_ft_ann, obs_ft_noise)
obs_ft_merged <- obs_ft_merged[match(obs_ft$id, table=obs_ft_merged$id),]

# Fix Sex
require(tableone)
obs_ft_merged$sex[obs_ft_merged$sex == 'm'] <- 'M'
obs_ft_merged$sex[obs_ft_merged$sex == 'f'] <- 'F'

tab1 <- CreateTableOne(vars = c('age', 'diagnosis', 'sex', 'wbc', 'hb', 'bm', 'anc', 'plt'), 
                       data = obs_ft_merged, 
                       strata = 'Cluster', 
                       factorVars = c('sex', 'diagnosis'))

sink('data/DescTable_HR_Merged.tsv')
print(tab1, nonnormal=c('wbc', 'hb', 'bm', 'anc', 'plt'))
sink()

obs_ft_merged$type <- ifelse(obs_ft_merged$id %in% obs_ft_ann$id, 'NonNoise', 'Noise')
fwrite(obs_ft_merged, file='data/HR_Annotated.tsv', sep = '\t', col.names = TRUE, row.names = FALSE, append=FALSE, quote=FALSE)

require(survival)
require(survminer)

obs_ft_merged$event[obs_ft_merged$os > 120] <- 0
obs_ft_merged$os[obs_ft_merged$os > 120] <- 120

matched_colors <- setNames(RColorBrewer::brewer.pal(length(unique(obs_ft_ann$Cluster)), 'Paired'), unique(obs_ft_ann$Cluster))

fit <- survfit(Surv(os, event)~Cluster, data=obs_ft_merged)
p <- ggsurvplot(fit, 
                data = obs_ft_merged, 
                palette = matched_colors, 
                color = 'Cluster', 
                conf.int = FALSE, 
                pval = TRUE, 
                risk.table = TRUE)

pdf('plots/Survival_HR_Updated.pdf', width = 10, height = 8)
print(p, newpage=FALSE)
dev.off()

# Survival Cross Comparison #
feat_ft <- fread('data/mut_feat_v4.tsv.gz', header=TRUE, sep = '\t')
mut_mat <- data.matrix(fread('data/mut_mat_v5_train_hr.tsv.gz', header=FALSE, sep='\t'))
obs_ft <- fread('data/mut_obs_v5_train_hr.tsv.gz', header=TRUE, sep='\t')
colnames(mut_mat) <- feat_ft$symbol
rownames(mut_mat) <- obs_ft$id
mut_mat_hr <- mut_mat
mut_mat <- data.matrix(fread('data/mut_mat_v5_train_lr.tsv.gz', header=FALSE, sep='\t'))
obs_ft <- fread('data/mut_obs_v5_train_lr.tsv.gz', header=TRUE, sep='\t')
colnames(mut_mat) <- feat_ft$symbol
rownames(mut_mat) <- obs_ft$id
mut_mat_lr <- mut_mat
hr_obs <- fread('data/HR_Annotated.tsv', header=TRUE, sep='\t')
lr_obs <- fread('data/LR_Annotated.tsv', header=TRUE, sep='\t')
hr_model <- readRDS('HR_RF_Model.rds')
lr_model <- readRDS('LR_RF_Model.rds')

hr_pred <- predict(hr_model, mut_mat_lr[match(lr_obs$id, table=rownames(mut_mat_lr)),])
lr_pred <- predict(lr_model, mut_mat_hr[match(hr_obs$id, table=rownames(mut_mat_hr)),])

# HR #
surv_ft <- rbind(data.table('os'=lr_obs$os,
                            'event'=lr_obs$event,
                            'cluster'=hr_pred,
                            'data'='LR_'),
                 data.table('os'=hr_obs$os,
                            'event'=hr_obs$event,
                            'cluster'=hr_obs$Cluster,
                            'data'='HR_')) %>%
  filter(cluster %in% paste0('Cluster-', c(1,5)))

surv_ft$event[surv_ft$os > 120] <- 0
surv_ft$os[surv_ft$os > 120] <- 120
surv_ft$type <- paste(surv_ft$data, surv_ft$cluster, sep='')

fit <- survfit(Surv(os, event)~type, data=surv_ft)
p <- ggsurvplot(fit, 
                pval = TRUE, 
                risk.table = TRUE)


pdf('plots/HR_Survival_LR_Updated.pdf', width = 12, height = 8)
print(p, newpage=FALSE)
dev.off()





local_obs <- fread('data/mut_obs_v5_train_lr_filtered_annotated.tsv.gz', header=TRUE, sep=',') %>%
  filter(Cluster %in% paste0('Cluster-', c(1,3,4)))





combined_obs <- rbind(data.table('os'=surv_ft$os,
                                 'event'=surv_ft$event,
                                 'type'=paste0('HR_', surv_ft$type)),
                      data.table('os'=local_obs$os,
                                 'event'=local_obs$event,
                                 'type'=paste0('LR_', ifelse(local_obs$Cluster %in% c('Cluster-1'), 'Good Prognosis',
                                                             ifelse(local_obs$Cluster %in% c('Cluster-3', 'Cluster-4'), 'Poor Prognosis', 'Filter'))))) %>%
  filter(!grepl(pattern='_Filter', type))




mut_mat <- mut_mat[match(obs_ft_ann$id, table=obs_ft$id),]

colnames(mut_mat) <- feat_ft$symbol
rownames(mut_mat) <- obs_ft_ann$id


fwrite(t(imp_res), file = 'data/HR_Importance.tsv.gz', col.names = TRUE, row.names = FALSE, sep='\t', append=FALSE, quote=FALSE)
data_ft <- as.data.table(mut_mat) %>%
  bind_cols(data.table('Cluster' = as.factor(obs_ft_ann$Cluster[match(obs_ft$id, table=obs_ft_ann$id)])))

rf_model <- rfsrc(Cluster~., data=data_ft, ntree=1500, nodesize=15, mtry=6, nodedepth = 3, importance='permute', splitrule = 'auc')

rf_model <- randomForest(x=mut_mat, 
                         y=factor(obs_ft_ann$Cluster), 
                         ntree=1500, 
                         importance=TRUE, 
                         nodesize=15, 
                         mtry=6)
saveRDS(rf_model, 'LR_RF_Model.rds')


# Validation Set #
val_mat <- data.matrix(fread('data/mut_mat_v5_validation.tsv.gz', header=FALSE, sep='\t'))
obs_ft <- fread('data/mut_obs_v5_validation.tsv.gz', header=TRUE, sep='\t')
feat_ft <- fread('data/mut_feat_v4.tsv.gz', header=TRUE)
colnames(val_mat) <- feat_ft$symbol
obs_ft_orig <- fread('data/mut_obs_v5_train_lr_filtered_annotated.tsv.gz', header=TRUE, sep=',')

rf_model <- readRDS('LR_RF_Model.rds')
idx <- obs_ft$diagnosis == 'LR-MDS'
pred <- predict(rf_model, val_mat[idx,])
local_obs <- obs_ft[idx,]
local_obs$cluster <- pred

require(survival)
require(survminer)
require(tidyverse)
matched_colors <- setNames(RColorBrewer::brewer.pal(length(unique(obs_ft_orig$Cluster)), 'Paired'), unique(obs_ft_orig$Cluster))

fit <- survfit(Surv(os, event)~cluster, data=local_obs)
p <- ggsurvplot(fit,  
                palette = matched_colors, 
                color = 'cluster', 
                conf.int = FALSE, 
                pval = TRUE, 
                risk.table = TRUE)

####
library(rms)
dd <- datadist(local_obs)
options(datadist = 'dd')

ph_fit <- rms::cph(Surv(os, event)~cluster, data=local_obs, x=TRUE, y=TRUE)
