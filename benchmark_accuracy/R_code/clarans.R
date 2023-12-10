library(mclust)
library(e1071)
library(aricode)
library(fastkmedoids)

BASE_PATH = "/u/parishar/scratch/DATASETS/exp_data/"


file_list <- c("BreastCancer.csv", "Gender.csv", "deng.csv", "pollen.csv")
data_list <- c("BreastCancer", "Gender", "deng", "pollen")
labels_list <- c("labels_BreastCancer.csv", "labels_Gender.csv", 
"labels_deng.csv",  "labels_pollen.csv")

# file_list <- c("Crop.csv", "Twitter.csv")
# data_list <- c("Crop", "Twitter")
# labels_list <- c("labels_Crop.csv", "labels_Twitter.csv")


rep <- 200
ari = c()

all_output <- data.frame()
avg_output <- data.frame()


set.seed(6)
seeds <- sample(200, rep, replace = FALSE)


for (i in 1:length(file_list)){

    filepath = file.path(BASE_PATH, file_list[i])
    data_name <- data_list[i]


    if (data_list[i] == "BreastCancer"){
        num_clusters = 3
    } else if (data_list[i] == "Gender"){
        num_clusters = 2
    } else if (data_list[i] == "deng"){
        num_clusters <- 8
    } else if (data_list[i] == "pollen"){
        num_clusters <- 10
    } else if (data_list[i] == "Crop"){
        num_clusters <- 5
    } else if (data_list[i] == "Twitter"){
        num_clusters <- 2
    }

    # Load the data
    data <- read.csv2(filepath, sep=",", stringsAsFactors=FALSE)
    labels = read.csv2(file.path(BASE_PATH, labels_list[i]))
  
    print(paste("DATA:", data_name))

    for (k in 1:rep){

        dd <- as.vector(dist(data))
        n <- nrow(data)

        clusters <- fastclarans(dd, n, k=num_clusters, numlocal=2, maxneighbor=0.025, seed=seeds[k])
        
        # Calculate the statistics
        ari <- adjustedRandIndex(labels[, 1], clusters@assignment)
        nmi <- NMI(labels[, 1], clusters@assignment)
        temp <- c("Clarans", data_name, num_clusters, 0, ari, nmi)

        all_output <- rbind(all_output, temp)

    }

    avg_ari = mean(as.numeric(all_output[(all_output[, 2] == data_name & all_output[, 3] == num_clusters), 5]))
    avg_nmi = mean(as.numeric(all_output[(all_output[, 2] == data_name & all_output[, 3] == num_clusters), 6]))
    
    max_ari = max(as.numeric(all_output[(all_output[, 2] == data_name & all_output[, 3] == num_clusters), 5]))
    max_nmi = max(as.numeric(all_output[(all_output[, 2] == data_name & all_output[, 3] == num_clusters), 6]))

    # Prepare average results
    temp = c("Clarans", data_name, num_clusters, 0, avg_ari, avg_nmi, max_ari, max_nmi) 
    avg_output = rbind(avg_output, temp)
}

colnames(all_output) <- c("Algorithm", "Data", "Clusters", "Iters", "ARI", "NMI")
colnames(avg_output) <- c("Algorithm", "Data", "Clusters", "Iters", "ARI", "NMI", "Max_ARI", "Max_NMI")


write.table(all_output, "Clarans_all_runs.csv", sep = ",", row.names = FALSE)
write.table(avg_output, "Clarans_avg_runs.csv", sep = ",", row.names = FALSE)







