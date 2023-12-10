library(mclust)
library(Seurat)
library(SeuratObject)
library(aricode)

args <- commandArgs(trailingOnly=TRUE)
BASE_PATH <- args[1]


# file_list = c("pollen_raw.csv", "darmanis_raw.csv", "usoskin_raw.csv", "mouse_pan.csv", 
#                           "Muraro_raw.csv", "QSDiaphragm_raw.csv")

# data_list <- c("Pollen", "Darmanis", "Usoskin", "Mouse_pan", "Muraro", "QSDiaphragm")

# labels_list <- c("labels_pollen.csv", "labels_darmanis.csv", "labels_usoskin.csv", 
#               "labels_mouse_pan.csv", "labels_Muraro.csv", "labels_QSDiaphragm.csv")


file_list <- c("pollen_raw.csv", "darmanis_raw.csv", "usoskin_raw.csv", "mouse_pan.csv", 
            "Muraro_raw.csv", "QSLimb_raw.csv", "QSTrachea_raw.csv", "QSLung_raw.csv", "QSDiaphragm_raw.csv", 
             "Q10XSpleen_raw.csv")

data_list <- c("Pollen", "Darmanis", "Usoskin", "Mouse_pan", "Muraro", "QSLimb", "QSTrachea", "QSLung", 
             "QSDiaphragm", "Q10XSpleen")

labels_list <- c("labels_pollen.csv", "labels_darmanis.csv", "labels_usoskin.csv", 
              "labels_mouse_pan.csv", "labels_Muraro.csv", "labels_QSLimb.csv", "labels_QSTrachea.csv", 
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
    raw_perplexity <- 4
    transformed_perplexity <- 3.4
    coh <- 12.999
    min_features <- 100
    ndims = 20
  } else if (data_list[i] == "Darmanis"){
    num_clusters <- 8
    raw_perplexity <- 0.9
    transformed_perplexity = 1
    coh <- 10.817
    min_features <- 100
    ndims <- 20 
  } else if (data_list[i] == "Usoskin"){
    num_clusters <- 4
    raw_perplexity <- 0.1
    transformed_perplexity = 0.4
    coh <- 8.945
    min_features <- 100
    ndims <- 20 
  } else if (data_list[i] == "Mouse_pan"){
    num_clusters <- 13
    raw_perplexity <- 0.9
    transformed_perplexity = 0.9
    coh <- 6.882
    min_features <- 100
    ndims <- 20 
  } else if (data_list[i] == "Muraro"){
    num_clusters <- 9
    raw_perplexity <- 0.2
    transformed_perplexity <- 0.25
    coh <- 10.751
    min_features <- 100
    ndims <- 20 
  } else if (data_list[i] == "QSDiaphragm"){
    num_clusters <- 5
    raw_perplexity <- 0.3
    transformed_perplexity = 0.3
    coh <- 8.143
    min_features <- 100
    ndims <- 20 
  } else if (data_list[i] == "QSLung"){
    num_clusters <- 11
    raw_perplexity <- 0.05
    transformed_perplexity = 0.03
    coh <- 6.99
    min_features <- 100
    ndims <- 20 
  } else if (data_list[i] == "QSTrachea"){
    num_clusters <- 4
    raw_perplexity <- 0.04
    transformed_perplexity = 0.025
    coh <- 8.1377
    min_features <- 100
    ndims <- 20 
  } else if (data_list[i] == "Q10XSpleen"){
    num_clusters <- 5
    raw_perplexity <- 0.04
    transformed_perplexity = 0.05
    coh <- 6.5776
    min_features <- 100
    ndims <- 20 
  } else if (data_list[i] == "QSLimb"){
    num_clusters <- 6
    raw_perplexity <- 0.16
    transformed_perplexity = 0.17
    coh <- 7.950
    min_features <- 100
    ndims <- 20 
  } 
  
  print(paste("Processing:", data_name))
  
  # Load the data
  filepath <- file.path(BASE_PATH, "raw_data", file_list[i])
  raw_data <- read.csv(filepath, sep=",", header=FALSE)

    labels <- read.csv(file.path(BASE_PATH, "raw_data" ,labels_list[i]), header=FALSE, sep=",")

      
    raw_data <- t(raw_data)

    colnames(raw_data) <- seq(1, dim(raw_data)[2], by=1)

    Seurat_data <- CreateSeuratObject(counts = raw_data, min.cells = 3, min.features = min_features)

    Seurat_data <- NormalizeData(Seurat_data, normalization.method = "LogNormalize", scale.factor = 10000)
    Seurat_data <- FindVariableFeatures(Seurat_data, selection.method = "vst")

    all.genes <- rownames(Seurat_data)
    Seurat_data <- ScaleData(Seurat_data, features = all.genes)

    Seurat_data <- RunPCA(Seurat_data, features = VariableFeatures(object = Seurat_data))
    Seurat_data <- FindNeighbors(Seurat_data, dims = 1:ndims)

    print(data_name)
    print(paste("Actual clusters: ", length(unique(labels[, 1])), sep=""))
    Seurat_data <- FindClusters(Seurat_data, resolution = raw_perplexity, random.seed = 12)

    
    # Calculate the statistics
    ari <- round(adjustedRandIndex(labels[, 1], Idents(Seurat_data)), 2)
    nmi <- round(NMI(labels[, 1], Idents(Seurat_data)), 2)
    
    temp <- c(data_name, "Seurat", num_clusters, length(unique(Idents(Seurat_data))), ari, nmi)
    final_output <- rbind(final_output, temp)


    ##### Run on transformed data
    transformed_data <- raw_data^(coh+1)
    transformed_data <- floor(log(transformed_data, base=coh))
    transformed_data[sapply(transformed_data, simplify = 'matrix', is.infinite)] <- 0

    Seurat_data <- CreateSeuratObject(counts = transformed_data, min.cells = 3, min.features = min_features)

    Seurat_data <- NormalizeData(Seurat_data, scale.factor = 10000)
    Seurat_data <- FindVariableFeatures(Seurat_data, selection.method = "vst")

    all.genes <- rownames(Seurat_data)
    Seurat_data <- ScaleData(Seurat_data, features = all.genes)

    Seurat_data <- RunPCA(Seurat_data, features = VariableFeatures(object = Seurat_data))
    Seurat_data <- FindNeighbors(Seurat_data, dims = 1:ndims)

    print(data_name)
    print(paste("Actual clusters: ", length(unique(labels[, 1])), sep=""))
    Seurat_data <- FindClusters(Seurat_data, resolution = transformed_perplexity, random.seed = 12)


    # Calculate the statistics
    ari <- round(adjustedRandIndex(labels[, 1], Idents(Seurat_data)), 2)
    nmi <- round(NMI(labels[, 1], Idents(Seurat_data)), 2)
    
    temp <- c(data_name, "Seurat-T", num_clusters, length(unique(Idents(Seurat_data))), ari, nmi)
    final_output <- rbind(final_output, temp)

}

colnames(final_output) <- c("Data", "Algorithm", "Clusters", "Detected_Clusters", "ARI", "NMI")
write.table(final_output, "Seurat_results.csv", sep = ",", row.names = FALSE)
