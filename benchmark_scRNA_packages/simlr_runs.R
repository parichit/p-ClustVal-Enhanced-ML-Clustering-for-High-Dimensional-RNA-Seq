library(aricode)
library(mclust)
library(SIMLR)


print(paste("Started execution on:", Sys.time()))


args <- commandArgs(trailingOnly=TRUE)
BASE_PATH <- args[1]



file_list <- c("pollen_raw.csv", "darmanis_raw.csv", "mouse_pan.csv",
            "QSLimb_raw.csv", "QSTrachea_raw.csv", 
            "QSLung_raw.csv", "QSDiaphragm_raw.csv", 
             "Q10XSpleen_raw.csv")

data_list <- c("Pollen", "Darmanis", "Mouse_pan", "QSLimb", 
                "QSTrachea", "QSLung", 
                "QSDiaphragm", "Q10XSpleen")

labels_list <- c("labels_pollen.csv", "labels_darmanis.csv", 
              "labels_mouse_pan.csv", "labels_QSLimb.csv", 
              "labels_QSTrachea.csv", 
              "labels_QSLung.csv", "labels_QSDiaphragm.csv",
              "labels_Q10XSpleen.csv")


# file_list <- c("Muraro_raw.csv")
# data_list <- c("Muraro")
# labels_list <- c("labels_Muraro.csv")


final_output <- data.frame()



for(i in 1:length(file_list)){
  

  data_name <- data_list[i]
  
  if (data_list[i] == "Pollen"){
    num_clusters <- 11 
    coh <- 12.999
  } else if (data_list[i] == "Darmanis"){
    num_clusters <- 8
    coh <- 10.817
  } else if (data_list[i] == "Usoskin"){
    num_clusters <- 4
    coh <- 8.945
  } else if (data_list[i] == "Mouse_pan"){
    num_clusters <- 13
    coh <- 6.882
  } else if (data_list[i] == "Muraro"){
    num_clusters <- 9
    coh <- 10.751
  } else if (data_list[i] == "QSDiaphragm"){
    num_clusters <- 5
    coh <- 8.143
  } else if (data_list[i] == "QSLung"){
    num_clusters <- 11
    coh <- 6.99
  } else if (data_list[i] == "QSTrachea"){
    num_clusters <- 4
    coh <- 8.1377
  } else if (data_list[i] == "Q10XSpleen"){
    num_clusters <- 5
    coh <- 6.5776
  } else if (data_list[i] == "QSLimb"){
    num_clusters <- 6
    coh <- 7.950
  }

  
    print(paste("Processing:", data_name))
  
    # Load the data
    filepath <- file.path(BASE_PATH, "raw_data", file_list[i])
    raw_data <- read.csv(filepath, sep=",", header=FALSE)

    labels <- read.csv(file.path(BASE_PATH, "raw_data" ,labels_list[i]), header=FALSE, sep=",")
    labels = labels[, 1]
      
    raw_data <- t(raw_data)
    
    # Run SIMLR
    if (data_name == "QSTrachea"| data_name == "QSLimb" | data_name == "QSDiaphragm"){
      res = SIMLR(X = raw_data, 
            c = length(unique(labels)), normalize = TRUE)
    } else{
      res = SIMLR_Large_Scale(X = raw_data, 
            c = length(unique(labels)), normalize = TRUE)
    }
    
    # Calculate the statistics
    ari <- round(adjustedRandIndex(labels, res$y$cluster), 2)
    nmi <- round(NMI(labels, res$y$cluster), 2)


    # Append the output
    temp <- c(data_name, "SIMLR", num_clusters, 
    length(unique(res$y$cluster)), ari, nmi)
    
    final_output <- rbind(final_output, temp)

    
    # Transform data
    print("Working on transformed data")

    transformed_data <- raw_data^(coh+1)
    transformed_data <- floor(log(transformed_data, base=coh))
    transformed_data[sapply(transformed_data, simplify = 'matrix', 
    is.infinite)] <- 0

    print("Transformation complete")


    # Run SIMLR
    if (data_name == "QSTrachea"| data_name == "QSLimb" | data_name == "QSDiaphragm"){
      res = SIMLR(X = transformed_data, 
            c = length(unique(labels)), normalize = TRUE)
    } else{
      res = SIMLR_Large_Scale(X = transformed_data, 
            c = length(unique(labels)), normalize = TRUE)
    }

    ari <- round(adjustedRandIndex(labels, res$y$cluster), 2)
    nmi <- round(NMI(labels, res$y$cluster), 2)
    

    # Append the output
    temp <- c(data_name, "T-SIMLR", num_clusters, 
    length(unique(res$y$cluster)), 
    ari, nmi)

    final_output <- rbind(final_output, temp)

}


colnames(final_output) <- c("Data", "Algorithm", "Clusters", "Detected_Clusters", "ARI", "NMI")
write.table(final_output, "SIMLR_results.csv", sep = ",", row.names = FALSE)



print(paste("Completed execution on:", Sys.time()))