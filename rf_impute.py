#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 12:14:00 2021

@author: durmaz
"""

import numpy as np
from sklearn.experimental import enable_hist_gradient_boosting
from sklearn.ensemble import HistGradientBoostingClassifier


def local_impute(x=None, y=None):
    idx = np.argsort(np.sum(np.isnan(x), axis=0))
    for i in idx:
        if (np.sum(np.isnan(x[::,i])) > 0) | (np.sum(np.isnan(y[::,i])) > 0):
            train_idx = np.where(np.isnan(x[::,i]) == False)[0]
            pred_idx = np.where(np.isnan(x[::,i]) == True)[0]
            local_model = HistGradientBoostingClassifier(categorical_features=np.arange(x.shape[1]-1), 
                                                         max_bins=2, 
                                                         max_depth=3,
                                                         loss = 'binary_crossentropy',
                                                         max_iter = 500).fit(x[train_idx,][::,np.setdiff1d(np.arange(x.shape[1]), np.asarray([i]))],
                                                                             x[train_idx,][::,i])
            pred_score = local_model.score(x[train_idx,][::,np.setdiff1d(np.arange(x.shape[1]), np.asarray([i]))],
                                           x[train_idx,][::,i])
            print("Prediction Accuracy: {}".format(pred_score))
            
            if pred_idx.shape[0] > 0:
                pred = local_model.predict(x[pred_idx,][::,np.setdiff1d(np.arange(x.shape[1]), np.asarray([i]))])
                x[pred_idx, i] = pred
            
            y_pred = local_model.predict(y[::,np.setdiff1d(np.arange(x.shape[1]), np.asarray([i]))])
            y_pred_idx = np.where(np.isnan(y[::,i]) == True)[0]
            
            if y_pred_idx.shape[0] > 0:
                y[y_pred_idx,i]=y_pred[y_pred_idx]
        
    return x, y
    


if __name__ == '__main__':    
    
    # Load Data
    in_wd = "~/Research/mds_latent"
    
    local_dat = np.genfromtxt('{}/data/mut_mat_v6_train.tsv.gz'.format(in_wd), delimiter='\t', filling_values=np.nan, dtype=np.float32)
    local_dat_test = np.genfromtxt('{}/data/mut_mat_v6_test.tsv.gz'.format(in_wd), delimiter='\t', filling_values=np.nan, dtype=np.float32)

    local_dat_imp = local_impute(local_dat.copy(), local_dat_test.copy())
    np.savetxt('{}/data/mut_mat_v6_train_imp.tsv.gz'.format(in_wd), X=local_dat_imp[0], delimiter='\t')
    np.savetxt('{}/data/mut_mat_v6_test_imp.tsv.gz'.format(in_wd), X=local_dat_imp[1], delimiter='\t')
    
    
