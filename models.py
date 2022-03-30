import tensorflow as tf
tf.compat.v1.disable_eager_execution()

import keras
from keras.layers import Input, Dense, Lambda
from keras.layers import BatchNormalization, Dropout, LeakyReLU
from keras.models import Model
from keras.losses import mse
from keras.optimizers import Adam
from keras import backend as K

from scipy import linalg as LA
import numpy as np
import pandas as pd
from sklearn.decomposition import KernelPCA

def PCA(data):
    
    data_raw = data
    #centering the data
    data -= np.mean(data, axis = 0)  

    cov = np.cov(data, rowvar = False)

    evals , evecs = LA.eigh(cov)

    # sort the eigenvalues (and eigenvectors accordingly) descending
    idx = np.argsort(evals)[::-1]
    evals = evals[idx]
    
    # sort the eigenvectors
    S = evecs[:,idx]
    
    A = np.dot(data, S) 
    
    np.testing.assert_array_almost_equal(np.matmul(A, LA.inv(S)), data_raw)
    
    # calculate loadings of each PCA component
    loadings = S * np.sqrt(evals)
    
    return([A, S, loadings])

    

class PCA_model:
    
    def __init__(self, reference_data, k):
        self.reference_data = reference_data
        self.k = k
        self.A, self.S, self.loadings = PCA(self.reference_data)
    
    def encode(self, target_data):

        # PCA encode target_data
        A2 = np.dot(target_data, self.S) 

        return(A2[:,:self.k])
    
    def reconstruct(self, target_data):

        # PCA encode target_data
        A2 = np.dot(target_data, self.S) 
        
        # reconstruct target_data based on learned S
        Xrec = np.matmul(A2[:, :self.k], LA.inv(self.S)[:self.k, ]) # Note: S.T == LA.inv(S) in this case, since eigen-vector matrix

        return(Xrec)
    
    def get_loadings(self):
        
        return(self.loadings[:,:self.k])



def KPCA(data, k, kernel, alpha, gamma, coef0, degree):

    data_raw = data
    
    #centering the data
    data -= np.mean(data, axis = 0)  

    kPCA = KernelPCA(n_components = k, kernel=kernel, alpha=alpha, gamma=gamma, coef0=coef0, degree=degree, fit_inverse_transform=True)

    kPCA.fit(data)
    
    return kPCA


class KPCA_model:
    
    def __init__(self, reference_data, k, kernel, alpha, gamma, coef0, degree):
        self.reference_data = reference_data
        self.k = k
        self.kernel= kernel
        self.alpha = alpha
        self.gamma = gamma
        self.coef0 = coef0 
        self.degree = degree
        self.model = KPCA(self.reference_data, self.k, self.kernel, self.alpha, self.gamma, self.coef0, self.degree)
    
    def encode(self, target_data):

        # KPCA transform target_data
        A2 = self.model.transform(target_data) 

        return(A2)
    
    def reconstruct(self, target_data):

        encoded = self.encode(target_data)

        return self.model.inverse_transform(encoded)


# Define sampling with reparameterization trick

class mtVAE():

    def __init__(self,
                input_shape,
                intermediate_dim,
                latent_dim,
                kl_beta=1e-2,
                learning_rate=1e-3):
        
        self.input_shape = input_shape
        self.intermediate_dim = intermediate_dim
        self.latent_dim = latent_dim
        self.kl_beta = kl_beta
        self.learning_rate = learning_rate
        
        # # =================
        # # Encoder
        # # =================

        # Definition
        self.input   = Input(shape=(self.input_shape,), name='encoder_input')
        self.x       = Dense(self.intermediate_dim)(self.input)
        self.x       = LeakyReLU()(self.x)


        self.mu      = Dense(self.latent_dim, name='latent_mu')(self.x)
        self.sigma   = Dense(self.latent_dim, name='latent_sigma')(self.x)

        # Define sampling with reparameterization trick
        def sample_z(args):
            mu, sigma = args
            batch     = K.shape(mu)[0]
            dim       = K.int_shape(mu)[1]
            eps       = K.random_normal(shape=(batch, dim))
            return mu + K.exp(sigma / 2) * eps


        # Use reparameterization trick 
        self.z       = Lambda(sample_z, output_shape=(self.latent_dim, ), name='z')([self.mu, self.sigma])

        # Instantiate encoder
        self.encoder = Model(self.input, [self.mu, self.sigma, self.z], name='encoder')

        # =================
        # Decoder
        # =================

        # Definition
        self.decoder_input   = Input(shape=(self.latent_dim, ), name='decoder_input')
        self.x     = Dense(self.intermediate_dim)(self.decoder_input)
        self.x     = LeakyReLU()(self.x)

        self.output  = Dense(self.input_shape)(self.x)

        # Instantiate decoder
        self.decoder = Model(self.decoder_input, self.output, name='decoder')

        # =================
        # VAE as a whole
        # =================

        # Instantiate VAE
        self.vae_outputs = self.decoder(self.encoder(self.input)[2])
        self.vae         = Model(self.input, self.vae_outputs, name='vae')


        # Define optimizer
        self.optimizer = Adam(learning_rate=self.learning_rate)
        
        # Define loss
        def kl_reconstruction_loss(true, pred):
            # Reconstruction loss
            reconstruction_loss = mse(true, pred)
            reconstruction_loss *= self.input_shape

            # KL divergence loss
            kl_loss = 1 + self.sigma - K.square(self.mu) - K.exp(self.sigma)
            kl_loss = K.sum(kl_loss, axis=-1)
            kl_loss *= -0.5

            # Total loss = 50% rec + 50% KL divergence loss
            return K.mean(reconstruction_loss + self.kl_beta*kl_loss)

        # Compile VAE
        self.vae.compile(optimizer=self.optimizer, loss=kl_reconstruction_loss, metrics = ['mse'])
    
    
    def train(self, train_data, val_data, n_epochs, batch_size, verbosity=2):
        
        self.vae.fit(train_data, train_data,
                     epochs = n_epochs, 
                     batch_size = batch_size, 
                     validation_data = (val_data, val_data),
                     verbose = verbosity)
    
    def encode(self, data):
        return self.encoder.predict(data)[2]
    
    def encode_mu(self, data):
        return self.encoder.predict(data)[0]
    
    def decode(self, data):
        return self.decoder.predict(data)
    
    def reconstruct(self, data):
        return self.decode(self.encode(data))
    
    def save_model(self, save_folder):
        
        vae_path     = save_folder + 'VAE.h5'
        encoder_path = save_folder + 'VAE_encoder.h5'
        decoder_path = save_folder + 'VAE_decoder.h5'
        
        self.vae.save(vae_path)
        self.encoder.save(encoder_path)
        self.decoder.save(decoder_path)
        
    
    def load_vae(self, save_path):

        # The two functions below have to be redefined for the loading
        # of the model. They cannot be methods of the mtVAE class for
        # some reason.
        # https://github.com/keras-team/keras/issues/13992
        
        # Define sampling with reparameterization trick
        def sample_z(args):
            mu, sigma = args
            batch     = K.shape(mu)[0]
            dim       = K.int_shape(mu)[1]
            eps       = K.random_normal(shape=(batch, dim))
            return mu + K.exp(sigma / 2) * eps
        
                # Define loss
        def kl_reconstruction_loss(true, pred):
            # Reconstruction loss
            reconstruction_loss = mse(true, pred)
            reconstruction_loss *= self.input_shape

            # KL divergence loss
            kl_loss = 1 + self.sigma - K.square(self.mu) - K.exp(self.sigma)
            kl_loss = K.sum(kl_loss, axis=-1)
            kl_loss *= -0.5

            # Total loss = 50% rec + 50% KL divergence loss
            return K.mean(reconstruction_loss + self.kl_beta*kl_loss)
        
        keras.losses.kl_reconstruction_loss = kl_reconstruction_loss
        self.vae = keras.models.load_model(save_path)
        self.vae.compile(optimizer=self.optimizer, 
                         custom_objects={'sample_z': sample_z}, 
                         loss=kl_reconstruction_loss, 
                         metrics = ['mse'])
        
    def load_encoder(self, save_path):
        self.encoder = keras.models.load_model(save_path)
        
    def load_decoder(self, save_path):
        self.decoder = keras.models.load_model(save_path)
        