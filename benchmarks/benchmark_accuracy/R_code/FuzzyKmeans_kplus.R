library(mclust)
library(e1071)
library(aricode)

args <- commandArgs(trailingOnly=TRUE)
BASE_PATH <- args[1]


# file_list = c("Pollen.csv", "Darmanis.csv", "Usoskin.csv", "mouse_pan.csv", 
#                           "Muraro.csv", "QSDiaphragm.csv")

# data_list <- c("Pollen", "Darmanis", "Usoskin", "Mouse_pan", "Muraro", "QSDiaphragm")

# labels_list <- c("labels_pollen.csv", "labels_darmanis.csv", "labels_usoskin.csv", 
#               "labels_mouse_pan.csv", "labels_Muraro.csv", "labels_QSDiaphragm.csv")

# centroid_list = c("PollenCentroids", "DarmanisCentroids", "UsoskinCentroids", 
#             "Mouse_panCentroids", "MuraroCentroids", "QSDiaphragmCentroids")


file_list <- c("pollen_raw.csv", "darmanis_raw.csv", "usoskin_raw.csv", "mouse_pan.csv", 
            "Muraro_raw.csv", "QSLimb_raw.csv", "QSTrachea_raw.csv", "QSLung_raw.csv", "QSDiaphragm_raw.csv", 
             "Q10XSpleen_raw.csv")

data_list <- c("Pollen", "Darmanis", "Usoskin", "Mouse_pan", "Muraro", "QSLimb", "QSTrachea", "QSLung", 
             "QSDiaphragm", "Q10XSpleen")

labels_list <- c("labels_pollen.csv", "labels_darmanis.csv", "labels_usoskin.csv", 
              "labels_mouse_pan.csv", "labels_Muraro.csv", "labels_QSLimb.csv", "labels_QSTrachea.csv", 
              "labels_QSLung.csv", "labels_QSDiaphragm.csv",
              "labels_Q10XSpleen.csv")


# file_list <- c("darmanis_raw.csv", "usoskin_raw.csv")
# data_list <- c("Darmanis", "Usoskin")
# labels_list <- c("labels_darmanis.csv", "labels_usoskin.csv")

# centroid_path <- file.path(BASE_PATH, "benchmark_1_centroids")

centroids <- read.csv2(file.path("kpp_centroid_indices.csv"), 
                      sep=",", stringsAsFactors = FALSE)

num_rep <- 19

all_output <- data.frame()
avg_output <- data.frame()


for(i in 1:length(file_list)){
  
  data_name <- data_list[i]
  
  
  if (data_list[i] == "Pollen"){
        num_clusters = 11 
      } else if (data_list[i] == "Darmanis"){
        num_clusters = 8
      } else if (data_list[i] == "Usoskin"){
        num_clusters = 4
      } else if (data_list[i] == "Mouse_pan"){
        num_clusters <- 13
      } else if (data_list[i] == "Muraro"){
        num_clusters <- 9
      } else if (data_list[i] == "QSDiaphragm"){
        num_clusters <- 5
      } else if (data_list[i] == "QSLimb"){
        num_clusters <- 6
      } else if (data_list[i] == "QSTrachea"){
        num_clusters <- 4
      } else if (data_list[i] == "QSLung"){
        num_clusters <- 11
      } else if (data_list[i] == "Q10XSpleen"){
        num_clusters <- 5
      }
  
    print(paste("Processing:", data_name))
  
    # Load the data
    filepath <- file.path(BASE_PATH, "raw_data", file_list[i])
    raw_data <- read.csv2(filepath, sep=",", stringsAsFactors=FALSE)

    labels <- read.csv2(file.path(BASE_PATH, "raw_data", labels_list[i]))


    for (rep in 0:num_rep){
      
      # Read the file with kpp centroids
      # raw_kpp_centroids = read.table(file.path(centroid_path, paste(centroid_list[i], "_raw_" , rep, ".csv" , sep="")), 
      #                       sep=",", stringsAsFactors = FALSE)

      # Extract the indices
      init_cen <- centroids[(centroids$Data == data_name & centroids$Run == rep), 3]
      init_cen <- as.numeric(unlist(strsplit(init_cen, split="+", fixed=TRUE)))

      # RAW Data
      results <- cmeans(raw_data, centers = raw_data[init_cen, ], iter.max = 100, method = "cmeans", m = 1.5)

      
      # Calculate the statistics
      ari <- adjustedRandIndex(labels[, 1], results$cluster)
      nmi <- NMI(labels[, 1], results$cluster)
      
      temp <- c(data_name, "FZKM", num_clusters, ari, nmi)
      all_output <- rbind(all_output, temp)
    }
    
    avg_ari_raw <- mean(as.numeric(all_output[(all_output[, 1] == data_name), 4]))
    ari_std <- sd(as.numeric(all_output[(all_output[, 1] == data_name), 4]))

    avg_nmi_raw <- mean(as.numeric(all_output[(all_output[, 1] == data_name), 5]))
    nmi_std <- sd(as.numeric(all_output[(all_output[, 1] == data_name), 5]))

    # Prepare average results
    temp <- c(data_name, "FZKM", num_clusters, avg_ari_raw, ari_std, avg_nmi_raw, nmi_std) 
    avg_output <- rbind(avg_output, temp)

}

colnames(all_output) <- c("Data", "Algorithm", "Clusters", "ARI", "NMI")
colnames(avg_output) <- c("Data", "Algorithm", "Clusters", "ARI", "ARI_STD", "NMI", "NMI_STD")

write.table(all_output, "FuzzyKM_kpp_all_results.csv", sep = ",", row.names = FALSE)
write.table(avg_output, "FuzzyKM_kpp_avg_results.csv", sep = ",", row.names = FALSE)
