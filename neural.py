from scipy.io import mmread
from scipy.sparse import save_npz
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE  # For 2D visualization
from keras.layers import Input, Dense, Layer
from keras.models import Model
import tensorflow.keras.backend as K
import os
import pickle

def assign_clusters(sample, k):
    path = 'data/sample'+str(sample)+'/matrix.mtx'
    sparse_matrix = mmread(path)
    X = sparse_matrix.toarray()
    X = X.transpose()

    # Preprocessing
    X = np.log1p(X)
    scaler = MinMaxScaler()
    X_scaled = scaler.fit_transform(X)

    input_dim = X_scaled.shape[1]
    encoding_dim = 32  # Latent space dimension
    n_clusters = k     # Adjust based on your data

    # Step 1: Pretrain Autoencoder
    input_layer = Input(shape=(input_dim,))
    encoder = Dense(encoding_dim, activation='relu')(input_layer)
    decoder = Dense(input_dim, activation='sigmoid')(encoder)
    autoencoder = Model(inputs=input_layer, outputs=decoder)
    autoencoder.compile(optimizer='adam', loss='mse')
    autoencoder.fit(X_scaled, X_scaled, epochs=100, batch_size=32, shuffle=True)

    # Extract encoded features
    encoder_model = Model(inputs=input_layer, outputs=encoder)
    encoded_data = encoder_model.predict(X_scaled)

    # Step 2: Initialize Clusters
    kmeans = KMeans(n_clusters=n_clusters, n_init=10)
    kmeans.fit(encoded_data)
    cluster_centers = kmeans.cluster_centers_

    # Step 3: Define Clustering Layer (same as before)
    class ClusteringLayer(Layer):
        def __init__(self, n_clusters, weights=None, alpha=1.0, **kwargs):
            super(ClusteringLayer, self).__init__(**kwargs)
            self.n_clusters = n_clusters
            self.alpha = alpha
            self.initial_weights = weights

        def build(self, input_shape):
            self.clusters = self.add_weight(
                shape=(self.n_clusters, input_shape[1]),
                initializer='glorot_uniform',
                name='clusters'
            )
            if self.initial_weights is not None:
                self.set_weights(self.initial_weights)
                del self.initial_weights
            self.built = True

        def call(self, inputs):
            q = 1.0 / (1.0 + (K.sum(K.square(K.expand_dims(inputs, axis=1) - self.clusters), axis=2) / self.alpha))
            q **= (self.alpha + 1.0) / 2.0
            q = K.transpose(K.transpose(q) / K.sum(q, axis=1))
            return q

    # Build and train clustering model (same as before)
    clustering_layer = ClusteringLayer(n_clusters, weights=[cluster_centers], name='clustering')(encoder)
    clustering_model = Model(inputs=input_layer, outputs=clustering_layer)
    clustering_model.compile(optimizer='adam', loss='kld')

    # Step 4: Define Target Distribution
    def target_distribution(q):
        weight = q ** 2 / q.sum(axis=0)
        return (weight.T / weight.sum(axis=1)).T

    # Train with target distribution
    max_iter = 20
    for epoch in range(max_iter):
        print('epoch '+str(epoch))
        q = clustering_model.predict(X_scaled, verbose=0)
        p = target_distribution(q)
        clustering_model.fit(X_scaled, p, epochs=1, batch_size=32, verbose=0)

    # Get final cluster labels
    q_final = clustering_model.predict(X_scaled, verbose=0)
    cluster_labels = np.argmax(q_final, axis=1)

    # Step 4: Plot Clusters in 2D using t-SNE
    tsne = TSNE(n_components=2, random_state=42)
    X_tsne = tsne.fit_transform(encoded_data)  # Use encoded features
    prefix_path = 'compressed/neural/sample'+str(sample)+'/k'+str(n_clusters)
    if not os.path.isdir(prefix_path):
        os.mkdir(prefix_path)

    plt.figure(figsize=(8, 6))
    plt.scatter(
        X_tsne[:, 0], X_tsne[:, 1],
        c=cluster_labels, cmap='viridis', alpha=0.8
    )
    plt.title('t-SNE Visualization of Cell Clusters (k='+str(n_clusters)+') for scRNA Matrix '+str(sample))
    plt.xlabel('t-SNE Dimension 1')
    plt.ylabel('t-SNE Dimension 2')
    plt.colorbar(label='Cluster')
    plt.savefig(prefix_path+'/clusters.png')
    plt.show()

    # Pickle the array and save to a file
    filename = prefix_path + '/cluster_labels.pkl'
    with open(filename, 'wb') as file:
        pickle.dump(cluster_labels, file)