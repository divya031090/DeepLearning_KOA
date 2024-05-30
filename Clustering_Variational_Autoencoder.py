
import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow.keras.layers import Input, Dense, Lambda, concatenate
from tensorflow.keras.models import Model
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.cluster import KMeans
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from tensorflow.keras import backend as K
from tensorflow.keras.layers import BatchNormalization

# Load data from CSV files into DataFrames

df_domain = pd.read_csv('metabolites.csv')
df_domain1 = pd.read_csv('RNA_synovial.csv')
df_domain2 = pd.read_csv('RNA_urine.csv')
df_domain3 = pd.read_csv('RNA_plasma.csv')

# Assuming both datasets have the same features (columns)

df_domain_without_id = df_domain.drop(columns=['ID'])
df_domain1_without_id = df_domain1.drop(columns=['ID'])
df_domain2_without_id = df_domain2.drop(columns=['ID'])
df_domain3_without_id = df_domain3.drop(columns=['ID'])


# Standardize data

X_domain = StandardScaler().fit_transform(df_domain_without_id.values)
X_domain1 = StandardScaler().fit_transform(df_domain1_without_id.values)
X_domain2 = StandardScaler().fit_transform(df_domain2_without_id.values)
X_domain3 = StandardScaler().fit_transform(df_domain3_without_id.values)

# Concatenate the two datasets
X_combined = np.concatenate((X_domain,X_domain1,X_domain2,X_domain3), axis=1)

# Split the combined data into training and testing sets
X_combined_train, X_combined_test = train_test_split(X_combined, test_size=0.3, random_state=42)

# Set random seeds for reproducibility
tf.random.set_seed(42)
np.random.seed(42)

# Define dimensions for RNA and domain domains

input_dim_domain = X_domain.shape[1]
input_dim_domain1 = X_domain1.shape[1]
input_dim_domain2 = X_domain2.shape[1]
input_dim_domain3 = X_domain3.shape[1]

# Define a shared encoder for both domains
input_dim_combined = X_combined_train.shape[1]
latent_dim_shared = 3



shared_encoder_input = Input(shape=(input_dim_combined,))
shared_encoder_hidden = Dense(256, activation='relu')(shared_encoder_input)
#shared_encoder_hidden = BatchNormalization()(shared_encoder_hidden)
shared_z_mean = Dense(latent_dim_shared)(shared_encoder_hidden)
shared_z_log_var = Dense(latent_dim_shared)(shared_encoder_hidden)

def sampling(args):
    shared_z_mean, shared_z_log_var = args
    batch = K.shape(shared_z_mean)[0]
    dim = K.int_shape(shared_z_mean)[1]
    epsilon = K.random_normal(shape=(batch, dim))
    return shared_z_mean + K.exp(0.5 * shared_z_log_var) * epsilon

shared_z = Lambda(sampling, output_shape=(latent_dim_shared,))([shared_z_mean, shared_z_log_var])

# Define separate decoders for RNA and metabolite domains
def build_decoder(input_dim, latent_dim):
    decoder_input = Input(shape=(latent_dim,))
    decoder_hidden = Dense(256, activation='relu')(decoder_input)
    decoder_output = Dense(input_dim, activation='sigmoid')(decoder_hidden)
    return Model(decoder_input, decoder_output)

# Build decoders for RNA and metabolite domains

decoder_domain = build_decoder(input_dim_domain, latent_dim_shared)
decoder_domain1 = build_decoder(input_dim_domain1, latent_dim_shared)
decoder_domain2 = build_decoder(input_dim_domain2, latent_dim_shared)
decoder_domain3 = build_decoder(input_dim_domain3, latent_dim_shared)

# Combine Encoders (shared for both domains)
combined_encoder = Model(shared_encoder_input, shared_z)

# Combine Decoders (separate for each domain)
combined_decoder_domain3 = decoder_domain3(combined_encoder.output)
combined_decoder_domain = decoder_domain(combined_encoder.output)
combined_decoder_domain1 = decoder_domain1(combined_encoder.output)
combined_decoder_domain2 = decoder_domain2(combined_encoder.output)

# Define Joint Autoencoder model
joint_autoencoder = Model(shared_encoder_input, [ combined_decoder_domain,combined_decoder_domain1,combined_decoder_domain2, combined_decoder_domain3])

# Compile Joint Autoencoder
joint_autoencoder.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.0001), loss='binary_crossentropy')

# Train the joint autoencoder
joint_autoencoder.fit(
    X_combined_train,
    [
        X_combined_train[:, :input_dim_domain3],
        X_combined_train[:, input_dim_domain3: input_dim_domain3 + input_dim_domain],
        X_combined_train[:, input_dim_domain3 + input_dim_domain: input_dim_domain3 + input_dim_domain + input_dim_domain1],
        X_combined_train[:, input_dim_domain3 + input_dim_domain + input_dim_domain1: input_dim_domain3 + input_dim_domain + input_dim_domain1 + input_dim_domain2],
    ],
    epochs=100,
    batch_size=16,
    validation_data=(
        X_combined_test,
        [
            X_combined_test[:, :input_dim_domain3],
            X_combined_test[:, input_dim_domain3: input_dim_domain3 + input_dim_domain],
            X_combined_test[:, input_dim_domain3 + input_dim_domain: input_dim_domain3 + input_dim_domain + input_dim_domain1],
            X_combined_test[:, input_dim_domain3 + input_dim_domain + input_dim_domain1: input_dim_domain3 + input_dim_domain + input_dim_domain1 + input_dim_domain2],
        ],
    ),
)
# Get the encoded representations for both domains
encoded_combined = combined_encoder.predict(X_combined)


###Evaluate optimal clusters

import sklearn
import numpy as np
import tensorflow as tf
from tensorflow.keras.layers import Input, Dense, Lambda
from tensorflow.keras.models import Model
from tensorflow.keras import backend as K
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import make_blobs
from sklearn.model_selection import train_test_split
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from plotly.subplots import make_subplots


from sklearn.metrics import silhouette_score
inertia = []
for k in range(1, 11):
    kmeans = KMeans(n_clusters=k, init='k-means++', random_state=42)
    kmeans.fit(encoded_combined)
    inertia.append(kmeans.inertia_)

# Plot the elbow curve
plt.plot(range(1, 11), inertia, marker='o')
plt.title('Elbow Method for Optimal Number of Clusters in VAE Latent Space')
plt.xlabel('Number of Clusters')
plt.ylabel('Inertia (Within-Cluster Sum of Squares)')
plt.show()

# Determine the optimal number of clusters using silhouette analysis
silhouette_scores = []
for k in range(2, 11):
    kmeans = KMeans(n_clusters=k, random_state=42)
    labels = kmeans.fit_predict(encoded_combined)
    silhouette_avg = silhouette_score(encoded_combined, labels)
    silhouette_scores.append(silhouette_avg)

# Plot silhouette scores
plt.plot(range(2, 11), silhouette_scores, marker='o')
plt.title('Silhouette Score for Optimal Number of Clusters in VAE Latent Space')
plt.xlabel('Number of Clusters')
plt.ylabel('Silhouette Score')
plt.show()



# Apply KMeans clustering to the encoded representations
nan_mask = np.isnan(encoded_combined)
encoded_combined[nan_mask] = np.nanmean(encoded_combined)
kmeans_combined = KMeans(n_clusters=3, init='k-means++', random_state=42)
clusters = kmeans_combined.fit_predict(encoded_combined)
# Create a DataFrame with ID and Cluster Labels
result_df = pd.DataFrame({'ID': df_domain3['ID'], 'Cluster_Labels': clusters})

#

# Display the cluster counts
cluster_counts = result_df['Cluster_Labels'].value_counts()
print(cluster_counts)


#Plot results in an interactive plot



from scipy.spatial.transform import Rotation
fig = make_subplots(rows=1, cols=1, specs=[[{'type': 'scatter3d'}]])

# Add a small offset to each cluster along each dimension
offset = 0.67


for cluster_label in np.unique(clusters):
    cluster_mask = (clusters == cluster_label)
    
    # Add an offset to each dimension for each cluster
    #offset_values = np.ones((np.sum(cluster_mask), 3)) * offset * cluster_label
    offset_values = (np.random.rand(np.sum(cluster_mask), 3) - 0.7) * offset * 0.85 + offset * cluster_label
    
    #if cluster_label > 0:
        # Move subsequent clusters towards the center
   #     offset_values = offset_values * (0.23 ** cluster_label)
   # else:
   #     offset_values = offset_values * 0.45

    if cluster_label == 0:
        # Move subsequent clusters towards the center
         offset_values = offset_values * 0.52
    elif cluster_label==1:
        offset_values = offset_values * 0.47
    elif cluster_label==2:
        offset_values = offset_values * 0.48
        
   

    trace = go.Scatter3d(
        x=encoded_combined[cluster_mask, 0] + offset_values[:, 0],
        y=encoded_combined[cluster_mask, 1] + offset_values[:, 1],
        z=encoded_combined[cluster_mask, 2] + offset_values[:, 2],
        mode='markers',
        name=f'Cluster {cluster_label}',
        marker=dict(size=8, opacity=0.8),
        text=[f'Cluster: {cluster}' for cluster in clusters[cluster_mask]],
    )
    fig.add_trace(trace)

# Update layout for better visualization
fig.update_layout(
    scene=dict(
        xaxis_title='Latent Dimension 1',
        yaxis_title='Latent Dimension 2',
        zaxis_title='Latent Dimension 3',
        aspectmode='cube'
    ),
    scene_camera=dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=0.08, y=2.2, z=0.08)  # Adjust the eye position for a better view
    ),
    height=700
)

# Add legend
fig.update_layout(legend=dict(x=0, y=1.1, orientation='h'))

# Show the interactive plot
fig.show('notebook')






