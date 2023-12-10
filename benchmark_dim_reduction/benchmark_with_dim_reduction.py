import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap.umap_ import UMAP as um
import pacmap
from sklearn.preprocessing import MinMaxScaler
import os
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics.cluster import adjusted_rand_score as ari, normalized_mutual_info_score as nmi


# Put the full path to the folder containing the CSV files here
base_path = "raw_data/"

# file_list = ["pollen_raw.csv", "darmanis_raw.csv", "usoskin_raw.csv", "mouse_pan.csv", 
#              "Muraro_raw.csv", "QSDiaphragm_raw.csv"]

# data_list = ["Pollen", "Darmanis", "Usoskin", "Mouse_pan", "Muraro", "QSDiaphragm"]

# label_list = ["labels_pollen.csv", "labels_darmanis.csv", "labels_usoskin.csv", 
#               "labels_mouse_pan.csv", "labels_Muraro.csv", "labels_QSDiaphragm.csv"]

file_list = ["pollen_raw.csv", "darmanis_raw.csv", "usoskin_raw.csv", "mouse_pan.csv", 
            "Muraro_raw.csv", "QSLimb_raw.csv", "QSTrachea_raw.csv", "QSLung_raw.csv", "QSDiaphragm_raw.csv", 
             "Q10XSpleen_raw.csv"]

data_list = ["Pollen", "Darmanis", "Usoskin", "Mouse_pan", "Muraro", "QSLimb", "QSTrachea", "QSLung", 
             "QSDiaphragm", "Q10XSpleen"]

label_list = ["labels_pollen.csv", "labels_darmanis.csv", "labels_usoskin.csv", 
              "labels_mouse_pan.csv", "labels_Muraro.csv", "labels_QSLimb.csv", "labels_QSTrachea.csv", 
              "labels_QSLung.csv", "labels_QSDiaphragm.csv",
              "labels_Q10XSpleen.csv"]


# file_list = ["pollen_raw.csv", "darmanis_raw.csv"]
# data_list = ["Pollen", "Darmanis"]
# label_list = ["labels_pollen.csv", "labels_darmanis.csv"]

# file_list = ["Q10XSpleen_raw.csv"]
# data_list = ["Q10XSpleen"]
# label_list = ["labels_Q10XSpleen.csv"]


def load_data(file_name, label_name):
    
    file_path = os.path.join(base_path, file_name)
    data = pd.read_csv(file_path, sep=",", header=None)
    labels = pd.read_csv(os.path.join(base_path, label_name), header=None)

    return data, labels


data_num_clusters = {"Usoskin": 4, "Pollen": 11, 
             "Mouse_pan": 13, "Darmanis": 8, "Muraro": 9,
             "QSLimb": 6, "QSLung": 11, "Q10XSpleen": 5, 
             "QSTrachea": 4, "QSDiaphragm": 5}



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
    new_data = np.emath.logn(coh_avg, new_data)

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
    # coh_list, sep_list = find_parameters(X_std, labels, n_nei, kmeans)
    coh_avg, sep_avg = find_parameters(X_std, n_nei)

    # Transform
    new_data =  transform_data(data, coh_avg, sep_avg)
    
    # Standardize the data
    data = np.log2(data+1)
    scaler = MinMaxScaler()
    X_std = scaler.fit_transform(data)

    # Clustering
    temp_raw_pca_ari = []
    temp_raw_tsne_ari = []
    temp_raw_umap_ari = []
    temp_raw_pca_nmi = []
    temp_raw_tsne_nmi = []
    temp_raw_umap_nmi = []
    
    temp_transformed_pca_ari = []
    temp_transformed_tsne_ari = []
    temp_transformed_umap_ari = []
    temp_transformed_pca_nmi = []
    temp_transformed_tsne_nmi = []
    temp_transformed_umap_nmi = []
    
    temp_alg = []
    temp_cluster = []
    temp_data = []
    
    
    # Perform PCA (keep components that provide 90% variance)
    pca = PCA(n_components=0.90)
    data_pca = pca.fit_transform(X_std)
    data_pca = pd.DataFrame(data_pca)
    
    pca = PCA(n_components=0.90)
    new_data_pca = pca.fit_transform(new_data)
    new_data_pca = pd.DataFrame(new_data_pca)
    
    
    
    # Perform T-SNE reduction
    data_tsne = TSNE(n_components=2, learning_rate='auto', init='random', random_state=10, perplexity=20).fit_transform(X_std)
    new_data_tsne = TSNE(n_components=2, learning_rate='auto', init='random', random_state=10, perplexity=20).fit_transform(new_data)
    
    data_tsne = pd.DataFrame(data_tsne)
    new_data_tsne = pd.DataFrame(new_data_tsne)
    
    
    # Perform UMAP reduction
    data_umap = um(n_neighbors=20, min_dist=0.0, n_components=2, random_state=100).fit_transform(X_std)
    new_data_umap = um(n_neighbors=20, min_dist=0.0, n_components=2, random_state=100).fit_transform(new_data)

    data_umap = pd.DataFrame(data_umap)
    new_data_umap = pd.DataFrame(new_data_umap)
    
    print(data_umap.shape, new_data_umap.shape)
    

    # Perform pcamap reduction
    embedding = pacmap.PaCMAP(n_components=2, n_neighbors=None, MN_ratio=0.5, FP_ratio=2.0) 
    data_pacmap = embedding.fit_transform(X_std, init="pca")
    new_data_pacmap = embedding.fit_transform(new_data, init="pca")
    
    data_pacmap = pd.DataFrame(data_pacmap)
    new_data_pacmap = pd.DataFrame(new_data_pacmap)


    
    for rep in range(num_rep):
        
        
        # PCA
        kmeans_pca = KMeans(n_clusters=num_clusters, random_state=seeds[i, rep], init="k-means++", n_init=5).fit(data_pca)
        
        temp_raw_pca_ari.append(ari(kmeans_pca.labels_, labels))
        temp_raw_pca_nmi.append(nmi(kmeans_pca.labels_, labels))
        temp_alg.append("PCA-Raw")
        
        
        # TSNE
        kmeans_tsne = KMeans(n_clusters=num_clusters, random_state=seeds[i, rep], init="k-means++", n_init=5).fit(data_tsne)
        
        temp_raw_tsne_ari.append(ari(kmeans_tsne.labels_, labels))
        temp_raw_tsne_nmi.append(nmi(kmeans_tsne.labels_, labels))
        temp_alg.append("TSNE-Raw")
        
        
        # UMAP 
        kmeans_umap = KMeans(n_clusters=num_clusters, random_state=seeds[i, rep], init="k-means++", n_init=5).fit(data_umap)
        
        temp_raw_umap_ari.append(ari(kmeans_umap.labels_, labels))
        temp_raw_umap_nmi.append(nmi(kmeans_umap.labels_, labels))
        temp_alg.append("UMAP-Raw")
        
        
        # PACMAP 
        # kmeans_pacmap = KMeans(n_clusters=num_clusters, random_state=seeds[i, rep], init="k-means++", n_init=5).fit(data_pacmap)
        
        # temp_raw_pacmap_ari.append(ari(kmeans_pacmap.labels_, labels))
        # temp_raw_pacmap_nmi.append(nmi(kmeans_pacmap.labels_, labels))
        # temp_alg.append("PACMAP-Raw")
        
        
        
        
        #### Transformed data
        
        kmeans_transformed_pca = KMeans(n_clusters=num_clusters, random_state=seeds[i, rep], init="k-means++", n_init=5).fit(new_data_pca)
        
        temp_transformed_pca_ari.append(ari(kmeans_transformed_pca.labels_, labels))
        temp_transformed_pca_nmi.append(nmi(kmeans_transformed_pca.labels_, labels))
        temp_alg.append("PCA-Transformed")
        
        
        kmeans_transformed_tsne = KMeans(n_clusters=num_clusters, random_state=seeds[i, rep], init="k-means++", n_init=5).fit(new_data_tsne)
        
        temp_transformed_tsne_ari.append(ari(kmeans_transformed_tsne.labels_, labels))
        temp_transformed_tsne_nmi.append(nmi(kmeans_transformed_tsne.labels_, labels))
        temp_alg.append("TSNE-Transformed")
        
        
        kmeans_transformed_umap = KMeans(n_clusters=num_clusters, random_state=seeds[i, rep], init="k-means++", n_init=5).fit(new_data_umap)
        
        temp_transformed_umap_ari.append(ari(kmeans_transformed_umap.labels_, labels))
        temp_transformed_umap_nmi.append(nmi(kmeans_transformed_umap.labels_, labels))
        temp_alg.append("UMAP-Transformed")
        
        
        # kmeans_transformed_pacmap = KMeans(n_clusters=num_clusters, random_state=seeds[i, rep], init="k-means++", n_init=5).fit(new_data_pacmap)
        
        # temp_transformed_pacmap_ari.append(ari(kmeans_transformed_pacmap.labels_, labels))
        # temp_transformed_pacmap_nmi.append(nmi(kmeans_transformed_pacmap.labels_, labels))
        # temp_alg.append("PACMAP-Transformed")
        
        
        temp_data.append(data_list[i])
        temp_data.append(data_list[i])
        temp_data.append(data_list[i])
        # temp_data.append(data_list[i])
        temp_data.append(data_list[i])
        temp_data.append(data_list[i])
        temp_data.append(data_list[i])
        # temp_data.append(data_list[i])
        
        temp_cluster.append(num_clusters)
        temp_cluster.append(num_clusters)
        temp_cluster.append(num_clusters)
        # temp_cluster.append(num_clusters)
        temp_cluster.append(num_clusters)
        temp_cluster.append(num_clusters)
        temp_cluster.append(num_clusters)
        # temp_cluster.append(num_clusters)



    temp = pd.DataFrame(list(zip(temp_data, temp_alg, temp_cluster, 
                                 temp_raw_pca_ari+temp_raw_tsne_ari+
                                 temp_raw_umap_ari+
                                 temp_transformed_pca_ari+temp_transformed_tsne_ari+
                                 temp_transformed_umap_ari,
                                 temp_raw_pca_nmi+temp_raw_tsne_nmi+
                                 temp_raw_umap_nmi+
                                 temp_transformed_pca_nmi+temp_transformed_tsne_nmi+
                                 temp_transformed_umap_nmi
                                 )), 
                        columns=["Data", "Algorithm", "Clusters", "ARI", "NMI"])
    
    
    all_results = pd.concat([all_results, temp], ignore_index=True, sort=False)

    out_alg.append("PCA-Raw")
    out_alg.append("TSNE-Raw")
    out_alg.append("UMAP-Raw")
    out_alg.append("PCA-Transformed")
    out_alg.append("TSNE-Transformed")
    out_alg.append("UMAP-Transformed")

    out_data.append(data_list[i])
    out_data.append(data_list[i])
    out_data.append(data_list[i])
    # out_data.append(data_list[i])
    out_data.append(data_list[i])
    out_data.append(data_list[i])
    out_data.append(data_list[i])
    # out_data.append(data_list[i])

    out_clusters.append(num_clusters)
    out_clusters.append(num_clusters)
    out_clusters.append(num_clusters)
    # out_clusters.append(num_clusters)
    out_clusters.append(num_clusters)
    out_clusters.append(num_clusters)
    out_clusters.append(num_clusters)
    # out_clusters.append(num_clusters)

    out_ari.append(np.mean(temp_raw_pca_ari))
    out_ari.append(np.mean(temp_raw_tsne_ari))
    out_ari.append(np.mean(temp_raw_umap_ari))
    out_ari.append(np.mean(temp_transformed_pca_ari))
    out_ari.append(np.mean(temp_transformed_tsne_ari))
    out_ari.append(np.mean(temp_transformed_umap_ari))
    
    out_ari_std.append(np.std(temp_raw_pca_ari))
    out_ari_std.append(np.std(temp_raw_tsne_ari))
    out_ari_std.append(np.std(temp_raw_umap_ari))
    out_ari_std.append(np.std(temp_transformed_pca_ari))
    out_ari_std.append(np.std(temp_transformed_tsne_ari))
    out_ari_std.append(np.std(temp_transformed_umap_ari))
    
    
    out_nmi.append(np.mean(temp_raw_pca_nmi))
    out_nmi.append(np.mean(temp_raw_tsne_nmi))
    out_nmi.append(np.mean(temp_raw_umap_nmi))
    out_nmi.append(np.mean(temp_transformed_pca_nmi))
    out_nmi.append(np.mean(temp_transformed_tsne_nmi))
    out_nmi.append(np.mean(temp_transformed_umap_nmi))
    
    
    out_nmi_std.append(np.std(temp_raw_pca_nmi))
    out_nmi_std.append(np.std(temp_raw_tsne_nmi))
    out_nmi_std.append(np.std(temp_raw_umap_nmi))
    out_nmi_std.append(np.std(temp_transformed_pca_nmi))
    out_nmi_std.append(np.std(temp_transformed_tsne_nmi))
    out_nmi_std.append(np.std(temp_transformed_umap_nmi))



# res = {"Data": out_data, 
#        "Algorithm": out_alg,
#        "Clusters": out_clusters,
#        "ARI": out_ari,
#        "ARI_STD": out_ari_std,
#        "NMI": out_nmi,
#        "NMI_STD": out_nmi_std}


avg_results = pd.DataFrame(list(zip(out_data, out_alg, out_clusters, out_ari, 
                                    out_ari_std, out_nmi, out_nmi_std)), 
                           columns=["Data", "Algorithm", "Clusters", "ARI", 
                                    "ARI_STD", "NMI", "NMI_STD"])


avg_results.to_csv("avg_results_benchmark_dim_reduction.csv", index=False)
all_results.to_csv("all_results_benchmark_dim_reduction.csv", index=False)
