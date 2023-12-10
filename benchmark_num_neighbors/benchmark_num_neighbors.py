import numpy as np
import pandas as pd
from sklearn.preprocessing import  MinMaxScaler
import os
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans
from sklearn.metrics.cluster import adjusted_rand_score as ari, normalized_mutual_info_score as nmi

import warnings
warnings.filterwarnings("ignore")


base_path = "/u/parishar/nobackup/DATASETS/exp_data/raw_data/"
centroid_path = "/u/parishar/nobackup/DATASETS/exp_data/benchmark_1_centroids/"


file_list = ["pollen_raw.csv", "usoskin_raw.csv", "mouse_pan.csv", 
            "QSTrachea_raw.csv", "Q10XSpleen_raw.csv"]

data_list = ["Pollen", "Usoskin", "Mouse_pan", "QSTrachea", "Q10XSpleen"]

label_list = ["labels_pollen.csv", "labels_usoskin.csv", 
              "labels_mouse_pan.csv", "labels_QSTrachea.csv", 
              "labels_Q10XSpleen.csv"]


# file_list = ["pollen_raw.csv", "darmanis_raw.csv", "usoskin_raw.csv", "mouse_pan.csv"]
# data_list = ["Pollen", "Darmanis", "Usoskin", "Mouse_pan"]
# label_list = ["labels_pollen.csv", "labels_darmanis.csv", "labels_usoskin.csv", 
#               "labels_mouse_pan.csv"]


# file_list = ["pollen_raw.csv"]
# data_list = ["Pollen"]
# label_list = ["labels_pollen.csv"]


# file_list = ["QSDiaphragm_raw.csv"]
# data_list = ["QSDiaphragm"]
# label_list = ["labels_QSDiaphragm.csv"]
# centroids_file = ["QSDiaphragmCentroids"]



data_num_clusters = {"Usoskin": 4, "Pollen": 11, 
             "Mouse_pan": 13, "Darmanis": 8, "Muraro": 9,
             "QSLimb": 6, "QSLung": 11, "Q10XSpleen": 5, 
             "QSTrachea": 4, "QSDiaphragm": 5}



def load_data(file_name, label_name):
    
    file_path = os.path.join(base_path, file_name)
    data = pd.read_csv(file_path, sep=",", header=None)
    labels = pd.read_csv(os.path.join(base_path, label_name), header=None)

    return data, labels


def pre_process_data(data, labels):

    data = np.array(data)

    print("Data Shape: ", data.shape)

    # Remove genes expressed in less than 1% of cell
    gene_expr_sum = np.sum(data, axis=0)
    limit = np.ceil(0.01 * data.shape[0])
    wch_genes = np.where(gene_expr_sum < limit)[0]
    
    if len(wch_genes) > 0:
        print("Genes to be removed:", len(wch_genes))
        data = np.delete(data, wch_genes, 1)

    rsum = np.sum(data, 1)
    wch_cells = np.where(rsum < 100)[0]
    

    if len(wch_cells) > 0:
        print("Cell to be removed: ", len(wch_cells))
        data = np.delete(data, wch_cells, 0)
        labels = np.delete(labels, wch_cells, 0)
        

    # Library size normalization
    data = (data/np.sum(data, 1)[:, None]) * np.median(np.sum(data, 1))
    labels = np.array(labels).reshape(data.shape[0],)
    data = pd.DataFrame(data)

    print("Data Shape: ", data.shape, len(labels))
    
    return data, labels



def find_parameters(mydata, n_nei):
    
    neigh = NearestNeighbors(n_neighbors=n_nei)
    neigh.fit(mydata)
    dists, _ = neigh.kneighbors(mydata, n_neighbors=n_nei, return_distance=True)

    coh_avg = dists[:, 0:n_nei].mean()
    sep_avg = coh_avg + 1

    return coh_avg, sep_avg



def transform_data(data, coh_avg, sep_avg):

    new_data = np.array(data[:])

    new_data = np.power(new_data, sep_avg)
    new_data = np.floor(np.emath.logn(coh_avg, new_data))

    new_data[new_data == np.inf] = 0
    new_data[new_data == -np.inf] = 0

    new_data = np.power(new_data, 2)
    new_data = pd.DataFrame(new_data)

    return new_data


# Lists to store results
out_neighbors = []
out_data = []
out_clusters = []
out_ari = []
out_nmi = []


num_rep = 20


# Generate seeds 
np.random.seed(9)
seeds = np.random.choice(2000, 200, replace=False).reshape(10, 20)


# Num neighbors
num_neighbors = [2, 5, 10, 15, 20, 25]
# num_neighbors = [2, 5, 10]


# Output pandas data frame
all_results = pd.DataFrame(columns=["Data", "Clusters", "Neighbors", "ARI", "NMI"])


for i in range(len(file_list)):

    print("Processing: ", data_list[i])

    file_name = file_list[i]
    label_name = label_list[i]
    
    # Num clusters
    num_clusters = data_num_clusters[data_list[i]]

    # Load
    data, labels = load_data(file_name, label_name)

    # Preprocess
    data, labels = pre_process_data(data, labels)

    # Standardize the data (to find parameters)
    scaler = MinMaxScaler()
    X_std = scaler.fit_transform(data)
    
    
    for n_nei in num_neighbors:
    
        # Find parameters
        coh_avg, sep_avg = find_parameters(X_std, n_nei)

        # DiVaS Transform
        new_data = transform_data(data, coh_avg, sep_avg)
        
        temp_transformed_ari = []
        temp_transformed_nmi = []
        
        temp_neighbors = []
        temp_cluster = []
        temp_data = []
    
        
        # Clustering  
        for rep in range(num_rep):
            
            kmeans = KMeans(n_clusters=num_clusters, random_state=seeds[i, rep], init="k-means++", n_init=5).fit(new_data)
                    
            temp_transformed_ari.append(ari(kmeans.labels_, labels))
            temp_transformed_nmi.append(nmi(kmeans.labels_, labels))
            
            temp_neighbors.append(n_nei)
            temp_data.append(data_list[i])
            temp_cluster.append(num_clusters)
    

        temp = pd.DataFrame(list(zip(temp_data, temp_cluster, temp_neighbors,
                                 temp_transformed_ari, temp_transformed_nmi)), 
                        
                        columns=["Data", "Clusters", "Neighbors", "ARI", "NMI"])
    
        all_results = pd.concat([all_results, temp], ignore_index=True, sort=False)

        out_data.append(data_list[i])
        out_neighbors.append(n_nei)
        out_clusters.append(num_clusters)

        out_ari.append(np.mean(temp_transformed_ari))
        out_nmi.append(np.mean(temp_transformed_nmi))



avg_results = pd.DataFrame(list(zip(out_data, out_clusters, out_neighbors, out_ari, out_nmi)), 
                           columns=["Data", "Clusters", "Neighbors", "ARI", "NMI"])


avg_results.to_csv("avg_results.csv", index=False)
all_results.to_csv("all_results.csv", index=False)

