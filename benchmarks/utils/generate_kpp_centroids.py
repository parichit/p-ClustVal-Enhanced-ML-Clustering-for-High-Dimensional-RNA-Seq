from sklearn.cluster import kmeans_plusplus
import pandas as pd
import numpy as np
import os

# basePath = "/Users/schmuck/Library/CloudStorage/OneDrive-IndianaUniversity/PhD/DATASETS/real_data/Data_without_labels/data_for_padic/"
# out_path = "/Users/schmuck/Library/CloudStorage/OneDrive-IndianaUniversity/PhD/DATASETS/real_data/Data_without_labels/Acc_centroids/"
# transformed_base_path = "/Users/schmuck/Library/CloudStorage/OneDrive-IndianaUniversity/PhD/DATASETS/real_data/Data_without_labels/data_for_padic/transformed_data/"

base_path = "/u/parishar/nobackup/DATASETS/exp_data/processed_data/"
transformed_path = "/u/parishar/nobackup/DATASETS/exp_data/transformed_data/"
std_path = "/u/parishar/nobackup/DATASETS/exp_data/std_data/"
out_path = "/u/parishar/nobackup/DATASETS/exp_data/benchmark_1_centroids/"


file_list = ["Pollen.csv", "Darmanis.csv", "Usoskin.csv", "mouse_pan.csv", 
             "Muraro.csv", "QSDiaphragm.csv"]

data_list = ["Pollen", "Darmanis", "Usoskin", "Mouse_pan", "Muraro", "QSDiaphragm"]

out_list = ["PollenCentroids", "DarmanisCentroids", "UsoskinCentroids", 
            "Mouse_panCentroids", "MuraroCentroids", "QSDiaphragmCentroids"]


# file_list = ["QSDiaphragm.csv"]
# data_list = ["QSDiaphragm"]
# out_list = ["QSDiaphragmCentroids"]


np.random.seed(9)
seeds = np.random.choice(2000, 120, replace=False).reshape(6, 20)

num_rep = 20


for i in range(len(file_list)):

    file = file_list[i]  
    
    path = os.path.join(base_path, file)
    raw_data = np.array(pd.read_csv(path, low_memory=False), dtype="float32")
    
    path = os.path.join(std_path, file)
    std_data = np.array(pd.read_csv(path, low_memory=False), dtype="float32")

    path = os.path.join(transformed_path, file)
    transformed_data = np.array(pd.read_csv(path, low_memory=False), dtype="float32")
    

    if data_list[i] == "Pollen":
        clus = 11

    elif (data_list[i] == "Darmanis"):
        clus = 8

    elif (data_list[i] == "Usoskin"):
        clus = 4

    elif (data_list[i] == "Mouse_pan"):
        clus = 13
        
    elif (data_list[i] == "Muraro"):
        clus = 6

    elif (data_list[i] == "QSDiaphragm"):
        clus = 5

    print("Processing: ", file)
    
    for rep in range(num_rep):

        centers_init, _ = kmeans_plusplus(raw_data, n_clusters=clus, random_state=seeds[i, rep])
        centers_init = pd.DataFrame(centers_init)

        centers_init.to_csv(os.path.join(out_path, out_list[i] + "_raw_" + str(rep) + ".csv"), sep=",", index=False, 
                            header=False)


        centers_init.to_csv(os.path.join(out_path, out_list[i] + "_transformed_" + str(rep) + ".csv"), sep=",", index=False, 
                            header=False)
        
        
        

