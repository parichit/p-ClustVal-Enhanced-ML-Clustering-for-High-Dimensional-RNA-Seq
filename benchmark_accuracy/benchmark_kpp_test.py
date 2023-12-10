import numpy as np
import pandas as pd
from sklearn.preprocessing import  MinMaxScaler, StandardScaler
import os
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from sklearn.mixture import GaussianMixture
from sklearn.metrics.cluster import adjusted_rand_score as ari, normalized_mutual_info_score as nmi

import warnings
warnings.filterwarnings("ignore")


base_path = "/u/parishar/nobackup/DATASETS/exp_data/raw_data/"
centroid_path = "/u/parishar/nobackup/DATASETS/exp_data/benchmark_1_centroids/"


file_list = ["pollen_raw.csv", "darmanis_raw.csv", "usoskin_raw.csv", "mouse_pan.csv", 
            "Muraro_raw.csv", "QSLimb_raw.csv", "QSTrachea_raw.csv", "QSLung_raw.csv", "QSDiaphragm_raw.csv", 
             "Q10XSpleen_raw.csv"]

data_list = ["Pollen", "Darmanis", "Usoskin", "Mouse_pan", "Muraro", "QSLimb", "QSTrachea", "QSLung", 
             "QSDiaphragm", "Q10XSpleen"]

label_list = ["labels_pollen.csv", "labels_darmanis.csv", "labels_usoskin.csv", 
              "labels_mouse_pan.csv", "labels_Muraro.csv", "labels_QSLimb.csv", "labels_QSTrachea.csv", 
              "labels_QSLung.csv", "labels_QSDiaphragm.csv",
              "labels_Q10XSpleen.csv"]


# file_list = ["darmanis_raw.csv", "pollen_raw.csv", "usoskin_raw.csv"]
# data_list = ["Darmanis", "Pollen", "Usoskin"]
# label_list = ["labels_darmanis.csv", "labels_pollen.csv", "labels_usoskin.csv"]
# centroids_file = ["DarmanisCentroids", "PollenCentroids", "UsoskinCentroids"]

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



def load_centroids(centroid_data, data_name, rep):
    
    ind = centroid_data[((centroid_data["Data"] == data_name) & (centroid_data["Run"] == rep))]["Indices"].values[0].split("+")
    ind = [int(i) for i in ind]
    return ind


def pre_process_data(data, labels, logged_version):

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
out_alg = []
out_data = []
out_clusters = []
out_ari = []
out_nmi = []
out_ari_std = []
out_nmi_std = []


num_rep = 20
n_nei = 5


# Generate seeds 
np.random.seed(9)
seeds = np.random.choice(2000, 200, replace=False).reshape(10, 20)

# Output pandas data frame
all_results = pd.DataFrame(columns=["Data", "Algorithm", "Clusters", "ARI", "NMI"])


# Read the centroids file
centroids_data = pd.read_csv("kpp_centroid_indices.csv", header=0)


for i in range(len(file_list)):

    print("Processing: ", data_list[i])

    file_name = file_list[i]
    label_name = label_list[i]

    # Load
    data, labels = load_data(file_name, label_name)

    # Preprocess
    data, labels = pre_process_data(data, labels, False)

    # Num clusters
    num_clusters = data_num_clusters[data_list[i]]

    # Standardize the data
    scaler = MinMaxScaler()
    X_std = scaler.fit_transform(data)
    
    # Find parameters
    coh_avg, sep_avg = find_parameters(X_std, n_nei)

    # Transform
    new_data =  transform_data(data, coh_avg, sep_avg)

    
    # Clustering
    temp_raw_ari_ac = []
    temp_raw_nmi_ac = []
    
    temp_transformed_ari_ac = []
    temp_transformed_nmi_ac = []
    
    temp_raw_ari_em = []
    temp_raw_nmi_em = []
    
    temp_transformed_ari_em = []
    temp_transformed_nmi_em = []
    
    temp_alg = []
    temp_cluster = []
    temp_data = []
    
    
    # Agglomarative clustering on raw data  
    for rep in range(num_rep):
                
        ac = AgglomerativeClustering(n_clusters=num_clusters).fit(data)
        
        temp_raw_ari_ac.append(ari(ac.labels_, labels))
        temp_raw_nmi_ac.append(nmi(ac.labels_, labels))
        temp_alg.append("WAC")
        
        
    # Agglomarative clustering on transformed data  
    for rep in range(num_rep):
        
        ac = AgglomerativeClustering(n_clusters=num_clusters).fit(new_data)
        
        temp_transformed_ari_ac.append(ari(ac.labels_, labels))
        temp_transformed_nmi_ac.append(nmi(ac.labels_, labels))
        temp_alg.append("T-WAC")
        
    
    # Expectation Maximization on raw data  
    for rep in range(num_rep):
        
        centroids_indices = load_centroids(centroids_data, data_list[i], rep)
        
        raw_labels = GaussianMixture(n_components=num_clusters, covariance_type="diag", 
                                     means_init=data.iloc[centroids_indices], n_init=5).fit_predict(data)
        
        temp_raw_ari_em.append(ari(raw_labels, labels))
        temp_raw_nmi_em.append(nmi(raw_labels, labels))
        temp_alg.append("EM")
        
    
    
    for rep in range(num_rep):
        
        centroids_indices = load_centroids(centroids_data, data_list[i], rep)
        
        raw_labels = GaussianMixture(n_components=num_clusters, covariance_type="diag", 
                                     means_init=data.iloc[centroids_indices], n_init=5).fit_predict(new_data)
        
        temp_transformed_ari_em.append(ari(raw_labels, labels))
        temp_transformed_nmi_em.append(nmi(raw_labels, labels))
        temp_alg.append("T-EM")

        
        temp_data.append(data_list[i])
        temp_data.append(data_list[i])
        temp_data.append(data_list[i])
        temp_data.append(data_list[i])
        temp_data.append(data_list[i])

        temp_cluster.append(num_clusters)
        temp_cluster.append(num_clusters)
        temp_cluster.append(num_clusters)
        temp_cluster.append(num_clusters)
        temp_cluster.append(num_clusters)
    
    # print(temp_raw_ari)
    # print(temp_tra_ari)

    temp = pd.DataFrame(list(zip(temp_data, temp_alg, temp_cluster, 
                        temp_raw_ari_ac+temp_transformed_ari_ac+
                        temp_raw_ari_em+temp_transformed_ari_em,
                                 
                        temp_raw_nmi_ac+temp_transformed_nmi_ac+
                        temp_transformed_nmi_em+temp_transformed_nmi_em)),
                        columns=["Data", "Algorithm", "Clusters", "ARI", "NMI"])
    
    all_results = pd.concat([all_results, temp], ignore_index=True, sort=False)

    out_alg.append("WAC")
    out_alg.append("T-WAC")
    out_alg.append("EM")
    out_alg.append("T-EM")

    out_data.append(data_list[i])
    out_data.append(data_list[i])
    out_data.append(data_list[i])
    out_data.append(data_list[i])

    out_clusters.append(num_clusters)
    out_clusters.append(num_clusters)
    out_clusters.append(num_clusters)
    out_clusters.append(num_clusters)


    out_ari.append(np.mean(temp_raw_ari_ac))
    out_ari.append(np.mean(temp_transformed_ari_ac))
    out_ari.append(np.mean(temp_raw_ari_em))
    out_ari.append(np.mean(temp_transformed_ari_em))
    
    
    out_ari_std.append(np.std(temp_raw_ari_ac))
    out_ari_std.append(np.std(temp_transformed_ari_ac))
    out_ari_std.append(np.std(temp_raw_ari_em))
    out_ari_std.append(np.std(temp_transformed_ari_em))
    
    
    out_nmi.append(np.mean(temp_raw_nmi_ac))
    out_nmi.append(np.mean(temp_transformed_nmi_ac))
    out_nmi.append(np.mean(temp_raw_nmi_em))
    out_nmi.append(np.mean(temp_transformed_nmi_em))
    
    
    out_nmi_std.append(np.std(temp_raw_nmi_ac))
    out_nmi_std.append(np.std(temp_transformed_nmi_ac))
    out_nmi_std.append(np.std(temp_raw_nmi_em))
    out_nmi_std.append(np.std(temp_transformed_nmi_em))



avg_results = pd.DataFrame(list(zip(out_data, out_alg, out_clusters, out_ari, out_ari_std, out_nmi, out_nmi_std)), 
                           columns=["Data", "Algorithm", "Clusters", "ARI", 
                                    "ARI_STD", "NMI", "NMI_STD"])



avg_results.to_csv("avg_results_kpp_accuracy_algorithms.csv", index=False)
all_results.to_csv("all_results_kpp_accuracy_algorithms.csv", index=False)

