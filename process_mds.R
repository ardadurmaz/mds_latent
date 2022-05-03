library(data.table)
library(skimr)
library(dplyr)
library(readxl)

mds <- fread('data/mds_final.csv', header=TRUE, sep=',', na.strings = c('', '-', 'NA', 'NaN', '--', 'N/A', 'none', 'NULL'), stringsAsFactors = FALSE, strip.white = TRUE)

obs_ft <- data.table('id'= paste0('sample.', as.integer(trimws(mds$`ID number`))),
                     'age' = as.numeric(trimws(mds$`Age @ sampling`)),
                     'diagnosis' = as.character(trimws(mds$`Diagnosis with Risk`)),
                     'sex' = as.character(trimws(mds$Sex)),
                     'os' = as.numeric(trimws(mds$`Overall survival time (months)`)),
                     'event' = as.integer(trimws(mds$`Status(death:1, alive:0) attenzione `)),
                     'bm' = as.numeric(trimws(mds$`BM blasts`)),
                     'wbc' = as.numeric(trimws(mds$`WBC(k/ul)`)),
                     'anc' = as.numeric(trimws(mds$ANC)),
                     'hb' = as.numeric(trimws(mds$`HB(g/dl)`)),
                     'plt' = as.numeric(trimws(gsub(pattern=',', replacement = '', mds$`PLT(k/ul)`))))

mds_ft <- data.table('normal' = as.integer(trimws(mds$Normal)),
                     'complex' = as.integer(trimws(mds$Complex)),
                     'del5' = as.integer(trimws(mds$`del(5)`)),
                     'del7' = as.integer(trimws(mds$`del(7)`)),
                     'del17' = as.integer(trimws(mds$`del(17)`)),
                     'del20' = as.integer(trimws(mds$`del(20)`)),
                     'trisomy8' = as.integer(trimws(mds$`(+8)`)),
                     'delY' = as.integer(trimws(mds$`del(Y)`)),
                     'other' = as.integer(trimws(mds$Other)),
                     'ASXL1' = as.integer(trimws(mds$ASXL1)),
                     'BCOR' = as.integer(trimws(mds$BCOR)),
                     'BCORL1' = as.integer(trimws(mds$BCORL1)),
                     'CALR' = as.integer(trimws(mds$CALR)),
                     'CBL' = as.integer(trimws(mds$CBL)),
                     'CEBPA' = as.integer(trimws(mds$CEBPA)),
                     'CSF3R' = as.integer(trimws(mds$CSF3R)),
                     'CUX1' = as.integer(trimws(mds$CUX1)),
                     'DDX41' = as.integer(trimws(mds$DDX41)),
                     'DNMT3A' = as.integer(trimws(mds$DNMT3A)),
                     'ETV6' = as.integer(trimws(mds$ETV6)),
                     'EZH2' = as.integer(trimws(mds$EZH2)),
                     'FLT3' = as.integer(trimws(mds$FLT3)),
                     'GATA2' = as.integer(trimws(mds$GATA2)),
                     'GNAS' = as.integer(trimws(mds$GNAS)),
                     'IDH1' = as.integer(trimws(mds$IDH1)),
                     'IDH2' = as.integer(trimws(mds$IDH2)),
                     'JAK2' = as.integer(trimws(mds$JAK2)),
                     'KIT' = as.integer(trimws(mds$KIT)),
                     'KRAS' = as.integer(trimws(mds$KRAS)),
                     'MPL' = as.integer(trimws(mds$MPL)),
                     'NF1' = as.integer(trimws(mds$NF1)),
                     'NOTCH1' = as.integer(trimws(mds$NOTCH1)),
                     'NPM1' = as.integer(trimws(mds$NPM1)),
                     'NRAS' = as.integer(trimws(mds$NRAS)),
                     'PHF6' = as.integer(trimws(mds$PHF6)),
                     'PTPN11' = as.integer(trimws(mds$PTPN11)),
                     'RAD21' = as.integer(trimws(mds$RAD21)),
                     'RUNX1' = as.integer(trimws(mds$RUNX1)),
                     'SETBP1' = as.integer(trimws(mds$SETBP1)),
                     'SF3B1' = as.integer(trimws(mds$SF3B1)),
                     'SMC1A' = as.integer(trimws(mds$SMC1A)),
                     'SMC3' = as.integer(trimws(mds$SMC3)),
                     'SRSF2' = as.integer(trimws(mds$SRSF2)),
                     'STAG2' = as.integer(trimws(mds$STAG2)),
                     'TET2' = as.integer(trimws(mds$TET2)),
                     'TP53' = as.integer(trimws(mds$TP53)),
                     'U2AF1' = as.integer(trimws(mds$U2AF1)),
                     'WT1' = as.integer(trimws(mds$WT1)),
                     'ZRSR2' = as.integer(trimws(mds$ZRSR2)))

mut_mat <- data.matrix(mds_ft)
mut_mat[mut_mat > 0] <- 1
rownames(mut_mat) <- obs_ft$id

# Merge Features #
mut_mat_init <- mut_mat[,!colnames(mut_mat) %in% c('BCOR', 'BCORL1', 'IDH1', 'IDH2', 'NRAS', 'KRAS', 'PTPN11', 'NF1', 'CBL', 'RAD21', 'SMC1A', 'SMC3', 'STAG2', 'CALR', 'MPL')]
mut_mat_second <- cbind(rowSums(mut_mat[,colnames(mut_mat) %in% c('BCOR', 'BCORL1')]),
                        rowSums(mut_mat[,colnames(mut_mat) %in% c('IDH1', 'IDH2')]),
                        rowSums(mut_mat[,colnames(mut_mat) %in% c('NRAS', 'KRAS', 'PTPN11', 'NF1', 'CBL')]),
                        rowSums(mut_mat[,colnames(mut_mat) %in% c('RAD21', 'SMC1A', 'SMC3', 'STAG2')]),
                        rowSums(mut_mat[,colnames(mut_mat) %in% c('CALR', 'MPL')]))
colnames(mut_mat_second) <- c('BCOR.BCORL1', 'IDH1.IDH2', 'NRAS.KRAS.PTPN11.NF1.CBL', 'RAD21.SMC1A.SMC3.STAG2', 'CALR.MPL')
mut_mat_combined <- cbind(mut_mat_init, mut_mat_second)
mut_mat_combined[mut_mat_combined > 1] <- 1

# Filter #
idx <- (rowSums(is.na(mut_mat_combined))/ncol(mut_mat_combined)) < 0.3
mut_mat_combined <- mut_mat_combined[idx,]
obs_ft <- obs_ft[idx,]

n_found <- TRUE
while(n_found){
  train_idx <- sample(1:nrow(mut_mat_combined), as.integer(nrow(mut_mat_combined)*0.8), replace = FALSE)
  if(all(colSums(mut_mat_combined, na.rm = TRUE) > 0)){
    n_found <- FALSE
  }
}
obs_ft$sex[obs_ft$sex == 'f'] <- 'F'
obs_ft$sex[obs_ft$sex == 'm'] <- 'M'
obs_ft$os[obs_ft$os > 150] <- NA
obs_ft$event[obs_ft$os > 120] <- 0
obs_ft$os[obs_ft$os > 120] <- 120

fwrite(mut_mat_combined[train_idx,], file='data/mut_mat_v6_train.tsv.gz', sep='\t', col.names = FALSE, row.names = FALSE, append=FALSE, quote=FALSE)
fwrite(obs_ft[train_idx,], file='data/mut_obs_v6_train.tsv.gz', sep='\t', col.names = TRUE, row.names = FALSE, append=FALSE, quote=FALSE)
fwrite(mut_mat_combined[-train_idx,], file='data/mut_mat_v6_test.tsv.gz', sep='\t', col.names = FALSE, row.names = FALSE, append=FALSE, quote=FALSE)
fwrite(obs_ft[-train_idx,], file='data/mut_obs_v6_test.tsv.gz', sep='\t', col.names = TRUE, row.names = FALSE, append=FALSE, quote=FALSE)
fwrite(data.table('symbol' = colnames(mut_mat_combined)), file='data/mut_feat_v6.tsv.gz', sep='\t', col.names = TRUE, row.names = FALSE, append=FALSE, quote=FALSE)


# ## Separate ##
# idx_hr <- obs_ft$diagnosis %in% c('HR-MDS', 's-AML')
# idx_lr <- obs_ft$diagnosis %in% c('LR-MDS')
# 
# freq_hr <- apply(mut_mat_combined[idx_hr,], 2, function(x){
#   sum(x, na.rm = TRUE)/(sum(!is.na(x)))
# })
# freq_lr <- apply(mut_mat_combined[idx_lr,], 2, function(x){
#   sum(x, na.rm = TRUE)/(sum(!is.na(x)))
# })
# 
# pdf('plots/FeatureHist.pdf')
# par(mfrow=c(2,1))
# hist(freq_hr, breaks=50)
# hist(freq_lr, breaks=50)
# dev.off()
# 
# idx <- freq_hr > 0.01 & freq_lr > 0.01
# mut_mat_combined <- mut_mat_combined[,idx]
# freq <- rowSums(is.na(mut_mat_combined))/ncol(mut_mat_combined)
# 
# pdf('plots/SampleHist.pdf')
# hist(freq, breaks=50)
# dev.off()
# 
# ## Combined ##
# mut_mat_combined <- cbind(mut_mat_combined, ifelse(obs_ft$diagnosis %in% c('HR-MDS', 's-AML'), 1, 0))
# colnames(mut_mat_combined)[ncol(mut_mat_combined)] <- 'HR'
# freq_res <- rowSums(is.na(mut_mat_combined))/ncol(mut_mat_combined)
# idx <- freq_res < 0.1
# mut_mat_combined <- mut_mat_combined[idx,]
# obs_ft <- obs_ft[idx,]
# 
# fwrite(mut_mat_combined, file='data/mut_mat_v6.tsv.gz', sep='\t', col.names = FALSE, row.names = FALSE, append=FALSE, quote=FALSE)
# fwrite(data.table('symbol' = colnames(mut_mat_combined)), file='data/mut_feat_v6.tsv.gz', sep='\t', col.names = TRUE, row.names = FALSE, append=FALSE, quote=FALSE)
# fwrite(obs_ft, file='data/mut_obs_v6.tsv.gz', sep='\t', col.names = TRUE, row.names = FALSE, append=FALSE, quote=FALSE)


# Validation Data
fix_val_bm <- function(s=NULL){
  s <- gsub(pattern='%|<', replacement = '', trimws(s))
  s_ft <- sapply(stringi::stri_split_fixed(s, pattern='-'), function(x){
    if(all(x != '') & all(!is.na(x))){
      if(length(x) == 1){
        return(as.numeric(x))
      }else if(length(x) == 2){
        return(mean(as.numeric(x), na.rm = TRUE))
      }else{
        return(NA)
      }
    }else{
      return(NA)
    }
  })
  return(s_ft)
}


mds <- fread('data/validation_data.csv', header=TRUE, sep=',', na.strings = c('', '-', 'NA', 'NaN', '--', 'N/A', 'none', 'NULL'), stringsAsFactors = FALSE, strip.white = TRUE)

val_obs_ft <- data.table('id'= paste0('val_sample.', as.integer(trimws(gsub(pattern='Validation', replacement = '', mds$`ID number`)))),
                         'age' = as.numeric(trimws(mds$`Age @ sampling`)),
                         'diagnosis' = as.character(trimws(mds$`Diagnosis with Risk`)),
                         'sex' = as.character(trimws(mds$Sex)),
                         'os' = as.numeric(trimws(mds$`Overall survival time (months)`)),
                         'event' = as.integer(trimws(mds$`Status(death:1, alive:0) attenzione `)),
                         'bm' = fix_val_bm(mds$`BM blasts`))

val_ft <- data.table('normal' = as.integer(trimws(mds$Normal)),
                     'complex' = as.integer(trimws(mds$Complex)),
                     'del5' = as.integer(trimws(mds$`del(5)`)),
                     'del7' = as.integer(trimws(mds$`del(7)`)),
                     'del17' = as.integer(trimws(mds$`del(17)`)),
                     'del20' = as.integer(trimws(mds$`del(20)`)),
                     'trisomy8' = as.integer(trimws(mds$`(+8)`)),
                     'delY' = as.integer(trimws(mds$`del(Y)`)),
                     'other' = as.integer(trimws(mds$Other)),
                     'ASXL1' = as.integer(trimws(mds$ASXL1)),
                     'BCOR' = as.integer(trimws(mds$BCOR)),
                     'BCORL1' = as.integer(trimws(mds$BCORL1)),
                     'CALR' = as.integer(trimws(mds$CALR)),
                     'CBL' = as.integer(trimws(mds$CBL)),
                     'CEBPA' = as.integer(trimws(mds$CEBPA)),
                     'CSF3R' = as.integer(trimws(mds$CSF3R)),
                     'CUX1' = as.integer(trimws(mds$CUX1)),
                     'DDX41' = as.integer(trimws(mds$DDX41)),
                     'DNMT3A' = as.integer(trimws(mds$DNMT3A)),
                     'ETV6' = as.integer(trimws(mds$ETV6)),
                     'EZH2' = as.integer(trimws(mds$EZH2)),
                     'FLT3' = as.integer(trimws(mds$FLT3)),
                     'GATA2' = as.integer(trimws(mds$GATA2)),
                     'GNAS' = as.integer(trimws(mds$GNAS)),
                     'IDH1' = as.integer(trimws(mds$IDH1)),
                     'IDH2' = as.integer(trimws(mds$IDH2)),
                     'JAK2' = as.integer(trimws(mds$JAK2)),
                     'KIT' = as.integer(trimws(mds$KIT)),
                     'KRAS' = as.integer(trimws(mds$KRAS)),
                     'MPL' = as.integer(trimws(mds$MPL)),
                     'NF1' = as.integer(trimws(mds$NF1)),
                     'NOTCH1' = as.integer(trimws(mds$NOTCH1)),
                     'NPM1' = as.integer(trimws(mds$NPM1)),
                     'NRAS' = as.integer(trimws(mds$NRAS)),
                     'PHF6' = as.integer(trimws(mds$PHF6)),
                     'PTPN11' = as.integer(trimws(mds$PTPN11)),
                     'RAD21' = as.integer(trimws(mds$RAD21)),
                     'RUNX1' = as.integer(trimws(mds$RUNX1)),
                     'SETBP1' = as.integer(trimws(mds$SETBP1)),
                     'SF3B1' = as.integer(trimws(mds$SF3B1)),
                     'SMC1A' = as.integer(trimws(mds$SMC1A)),
                     'SMC3' = as.integer(trimws(mds$SMC3)),
                     'SRSF2' = as.integer(trimws(mds$SRSF2)),
                     'STAG2' = as.integer(trimws(mds$STAG2)),
                     'TET2' = as.integer(trimws(mds$TET2)),
                     'TP53' = as.integer(trimws(mds$TP53)),
                     'U2AF1' = as.integer(trimws(mds$U2AF1)),
                     'WT1' = as.integer(trimws(mds$WT1)),
                     'ZRSR2' = as.integer(trimws(mds$ZRSR2)))

val_mat <- data.matrix(val_ft)

# Merge Features #
val_mat_init <- val_mat[,!colnames(val_mat) %in% c('BCOR', 'BCORL1', 'IDH1', 'IDH2', 'NRAS', 'KRAS', 'PTPN11', 'NF1', 'CBL', 'RAD21', 'SMC1A', 'SMC3', 'STAG2', 'CALR', 'MPL')]
val_mat_second <- cbind(rowSums(val_mat[,colnames(val_mat) %in% c('BCOR', 'BCORL1')]),
                        rowSums(val_mat[,colnames(val_mat) %in% c('IDH1', 'IDH2')]),
                        rowSums(val_mat[,colnames(val_mat) %in% c('NRAS', 'KRAS', 'PTPN11', 'NF1', 'CBL')]),
                        rowSums(val_mat[,colnames(val_mat) %in% c('RAD21', 'SMC1A', 'SMC3', 'STAG2')]),
                        rowSums(val_mat[,colnames(val_mat) %in% c('CALR', 'MPL')]))
colnames(val_mat_second) <- c('BCOR.BCORL1', 'IDH1.IDH2', 'NRAS.KRAS.PTPN11.NF1.CBL', 'RAD21.SMC1A.SMC3.STAG2', 'CALR.MPL')
val_mat_combined <- cbind(val_mat_init, val_mat_second)
val_mat_combined[val_mat_combined > 1] <- 1



feat_ft <- fread('data/mut_feat_v6.tsv.gz', header=TRUE)
val_mat_combined <- val_mat_combined[,match(feat_ft$symbol, table=colnames(val_mat_combined))]
val_mat_combined[is.na(val_mat_combined)] <- 0
fwrite(val_mat_combined, file='data/mut_mat_v6_validation.tsv.gz', sep='\t', col.names = FALSE, row.names = FALSE, append=FALSE, quote=FALSE)
fwrite(val_obs_ft, file='data/mut_obs_v6_validation.tsv.gz', sep='\t', col.names = TRUE, row.names = FALSE, append=FALSE, quote=FALSE)
