# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 14:00:42 2021

@author: durmaz
"""


from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import numpy as np
from sklearn import mixture
import tensorflow as tf


class CustDenseTranspose(tf.keras.layers.Layer):
    def __init__(self, dense=None, l2=0.01, activation='linear', use_bias=True, **kwargs):
        super(CustDenseTranspose, self).__init__()
        self.dense = dense
        self.l2 = tf.Variable(name='l2_reg', initial_value=l2, trainable=False, dtype=tf.float32)
        self.act = activation
        self.use_bias = use_bias
            
    def build(self, input_shape):
        if self.use_bias:
            self.b = self.add_weight(shape=(tf.shape(self.dense.weights[0])[0],), initializer='zeros', trainable=True)
        else:
            self.b = self.add_weight(shape=(tf.shape(self.dense.weights[0])[0],), initializer='zeros', trainable=False)
    
    def call(self, inputs):
        l2_loss = self.l2 * tf.math.reduce_sum(tf.math.square(self.dense.weights[0]))
        self.add_loss(l2_loss)
        
        if self.act == 'linear':
            return tf.math.add(tf.matmul(inputs, self.dense.weights[0], transpose_b=True), self.b)
        elif self.act == 'tanh':
            return tf.math.tanh(tf.math.add(tf.matmul(inputs, self.dense.weights[0], transpose_b=True), self.b))
        elif self.act == 'sigmoid':
            return tf.math.sigmoid(tf.math.add(tf.matmul(inputs, self.dense.weights[0], transpose_b=True), self.b))


def cust_loss(y_true, y_pred):
    
    y_pred = tf.clip_by_value(y_pred, tf.keras.backend.epsilon(), 1.0-tf.keras.backend.epsilon())
    log_loss = tf.math.add(tf.math.multiply(y_true, tf.math.log(y_pred))*2.0, tf.math.multiply(1.0-y_true, tf.math.log(1.0-y_pred))*0.5)
    
    return -1.0*tf.math.reduce_sum(tf.math.reduce_mean(log_loss, axis=-1), axis=-1)



class AutoEncoder(tf.keras.Model):

    def __init__(self, enc, dec, name="autoencoder"):
        super(AutoEncoder, self).__init__(name=name)
        self.encoder = enc
        self.decoder = dec

    def call(self, inputs):
        # Reconstruction Loss
        z = self.encoder(inputs)
        reconstructed = self.decoder(z)
                
        return reconstructed


class EncoderModel(tf.keras.Model):

    def __init__(self, latent_dim, name="encoder"):
        super(EncoderModel, self).__init__(name=name)
        self.z = tf.keras.layers.Dense(latent_dim, 
                                       activation='linear', 
                                       name='mut_lat', 
                                       kernel_regularizer=tf.keras.regularizers.L1L2(l1=1e-3, l2=1e-3),
                                       use_bias=False)


    def call(self, inputs):
        lat_mut = self.z(inputs)
        
        return lat_mut


class DecoderModel(tf.keras.Model):

    def __init__(self, dense=None, name="decoder"):
        super(DecoderModel, self).__init__(name=name)
        self.out_mut = CustDenseTranspose(dense=dense[0], 
                                          l2=0.0, 
                                          activation='sigmoid', 
                                          use_bias=False, 
                                          name='mut_out')
        
                
    def call(self, inputs):
        d_mut = self.out_mut(inputs)
        
        return d_mut


def sp_ae(data=None, dim=None, doCV=False, n_epoch=None, keep=False):
    
    # Autoencoder
    BATCH_SIZE = 32
    if n_epoch is None:
        NUM_EPOCHS = 10000
    else:
        NUM_EPOCHS = n_epoch
        
    LAT_DIM = dim
    N_DIM_MUT=data.shape[1]
            
    # Model
    enc = EncoderModel(LAT_DIM)
    dec = DecoderModel(enc.layers)
    vae = AutoEncoder(enc, dec)
    
    # Callback
    lr_callback = tf.keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.1, patience=10, verbose=0, mode='auto', min_lr=1e-12)
    ear_callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=15, min_delta=1e-12, mode='auto', restore_best_weights=True)

    # Compile
    optimizer = tf.keras.optimizers.Adam(learning_rate=1e-4)
    vae.compile(optimizer=optimizer, loss=[cust_loss])
    rand_idx = np.random.permutation(np.arange(data.shape[0]))
    train_hist = vae.fit(data[rand_idx,], data[rand_idx,], epochs=NUM_EPOCHS, batch_size=BATCH_SIZE, shuffle=True, validation_split=0.1, callbacks=[lr_callback, ear_callback])
                    
    # Clear
    tf.keras.backend.clear_session()
    
    # Embedding
    z = enc.predict(data)

    return z, train_hist.history, enc, dec
    

def local_gmm(x=None):
    
    try:
        bic = []
        clust = []
        cv_types = ["spherical", "tied", "diag", "full"]
        for cv_type in cv_types:
            for n_components in range(1,20):
                # Fit a Gaussian mixture with EM
                gmm = mixture.GaussianMixture(
                    n_components=n_components, covariance_type=cv_type
                    )
                clust.append(gmm.fit_predict(x))
                bic.append(gmm.bic(x))
            
        bic = np.asarray(bic)
        clust = np.vstack(clust)
    
        clust_res = clust[np.argmin(bic),]
    
        return clust_res
    
    except ValueError:
        print("!Error in GMM")
        
        return None
    

if __name__ == '__main__':    
    
    in_wd = "~/Research/mds_latent"
    k_fold = sys.argv[1]
    
    mut_mat = np.loadtxt('{}/data/mut_mat_v6_fold{}.tsv.gz'.format(in_wd, k_fold), delimiter='\t', dtype=np.float32)
    
    #val_err = []
    #for d in np.asarray([16, 32, 64, 128]):
    #    embedding, hist, enc, dec = sp_ae(mut_mat, dim=d, n_epoch=1750, keep=False)
    #    val_err.append(np.asarray(hist['val_loss']))
    
    #fig, ax = plot.subplots(nrows=1, ncols=1)
    #ax.plot(np.vstack(val_err).T, labels=['Dim-{}'.format(i) for i in np.asarray([16,32,64,128])], legend='t')
    #fig.savefig('{}/plots/AE_CV_Res.pdf'.format(in_wd))
    #plot.show()
            
    # AE
    cc_mat = np.zeros((mut_mat.shape[0], mut_mat.shape[0]), dtype=np.float32)
    fr_mat = np.zeros((mut_mat.shape[0], mut_mat.shape[0]), dtype=np.float32)

    for i in range(100):
        local_rand = np.random.choice(np.arange(mut_mat.shape[0]), size=int(0.75*mut_mat.shape[0]), replace=False)
        embedding, hist, enc, dec = sp_ae(mut_mat[local_rand,], dim=32, n_epoch=7500, keep=False)
        
        # GMM
        clust_res = local_gmm(embedding)
        if clust_res is not None:
            clust_mat = np.zeros((mut_mat.shape[0], np.max(clust_res)+1))
            for j in range(local_rand.shape[0]):
                clust_mat[local_rand[j],] = np.eye(np.max(clust_res)+1)[clust_res[j],]
            cc_mat += np.dot(clust_mat, clust_mat.T)
            
            freq_mat = np.zeros((mut_mat.shape[0], 2))
            for j in range(local_rand.shape[0]):
                freq_mat[local_rand[j],] = np.eye(2)[0] 
            fr_mat += np.dot(freq_mat, freq_mat.T)
        
    freq_mat = cc_mat/(fr_mat+1)
    
    np.save('{}/data/pca_cc_res_fold{}.npy'.format(in_wd, k_fold), arr=freq_mat)
    np.savetxt('{}/data/pca_cc_res_fold{}.tsv.gz'.format(in_wd, k_fold), X=freq_mat, delimiter='\t')
