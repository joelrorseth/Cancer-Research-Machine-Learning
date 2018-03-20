#
# Stage Classification
# Attempt to classify patients' stage of breast cancer given around 
# 24,000 gene expressions and 1400 patient samples.
#
# Joel Rorseth
# February 18, 2018
#

library(class)
library(tree)
library(e1071)


classify_selected <- function(X, y, selected_genes) {
  
  # Get feature selected X
  selected_genes = readLines("genes_rfe_top_all.txt")[0:100]
  new_genetic_data <- X[, colnames(X) %in% selected_genes]
  
  # Split up data to get train / test dataset
  train_filter = sort(sample(nrow(new_genetic_data), nrow(new_genetic_data)*.7))
  train.X <- new_genetic_data[train_filter,]
  test.X <- new_genetic_data[-train_filter,]
  train.y <- y[train_filter]
  test.y <- y[-train_filter]
  
  
  # SECTION 1: KNN
  knn_class <- knn(train.X, test.X, as.factor(train.y), k=3)
  
  # Print confusion matrix and test classification score
  # print( table(knn_class, test.y) )
  print( mean(knn_class == test.y) ) # Mean classification score
  
  
  # SECTION 2: Try SVM
  #svm_clf <- svm(as.factor(patient_stages)~., data=data.matrix(new_genetic_data), 
  #               kernel="polynomial", cost=10, scale=FALSE)
  #plot(svm_clf, data.matrix(new_genetic_data))
  
  
  # SECTION 3: Random Forest
  #tree.patients <- tree(as.factor(patient_stages)~., new_genetic_data)
  
  #plot(tree.patients)
  #text(tree.patients, pretty=0) 
}



# Read preprocessed patient data
setwd("/Users/Joel/Documents/Research/Cancer Stage Diagnosis/Datasets")
patient_data <- read.csv("anon_patient_expressions_stage_all.csv", header=T, fill=T)

# Separate predictors and response, X and y
patient_stages <- patient_data[,ncol(patient_data)]
patient_data <- patient_data[,-ncol(patient_data)]



# Section: RFE
# Try KNN on the top n genes
for (n in c(5, 10, 50, 100, 200, 500, 1000)) {
  
  # Get feature selected X
  selected_genes = readLines("genes_rfe_top_all.txt")[0:n]
  new_genetic_data <- patient_data[, colnames(patient_data) %in% selected_genes]
  
  # Split up data to get train / test dataset
  train_filter = sort(sample(nrow(new_genetic_data), nrow(new_genetic_data)*.7))
  train.X <- new_genetic_data[train_filter,]
  test.X <- new_genetic_data[-train_filter,]
  train.y <- patient_stages[train_filter]
  test.y <- patient_stages[-train_filter]
  
  
  # Try KNN
  knn_class <- knn(train.X, test.X, as.factor(train.y), k=3)
  
  # Print confusion matrix and test classification score
  print( mean(knn_class == test.y) )
}




for (selection in c("genes_tree_top_50.txt", "genes_tree_top_100.txt", 
                    "genes_tree_top_500.txt", "genes_tree_top_1000.txt")) {
  
  # Read in list of feature selected genes
  selected_genes <- readLines(selection)
  
  # Run classifiers on this feature subset
  classify_selected(patient_data, patient_stages, selected_genes)
}
