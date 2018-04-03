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


# Perform classification using feature selection
classify_selected <- function(X, y, filename_selected, limit) {
  
  # Get feature selected X using feature selected genes provided
  selected_genes = readLines(filename_selected)[0:limit]
  new_genetic_data <- X[, colnames(X) %in% selected_genes]
  
  # Split up data to get train / test dataset
  train_filter = sort(sample(nrow(new_genetic_data), nrow(new_genetic_data)*.7))
  train.X <- new_genetic_data[train_filter,]
  test.X <- new_genetic_data[-train_filter,]
  train.y <- y[train_filter]
  test.y <- y[-train_filter]
  
  
  # SECTION 1: KNN
  #knn_class <- knn(train.X, test.X, as.factor(train.y), k=3)
  
  # Print confusion matrix and test classification score
  # print( table(knn_class, test.y) )
  #print( mean(knn_class == test.y) ) # Mean classification score
  
  
  # SECTION 2: Try SVM
  #svm_clf <- svm(as.factor(patient_stages)~., data=data.matrix(new_genetic_data), 
  #               kernel="polynomial", cost=10, scale=FALSE)
  #plot(svm_clf, data.matrix(new_genetic_data))
  
  
  # SECTION 3: Random Forest
  tree.patients <- tree(as.factor(y)~., new_genetic_data, subset=train_filter)
  tree.pred <- predict(tree.patients, test.X, type="class")
  
  tab <- table(tree.pred, test.y)
  print(mean(tree.pred==test.y))
  
  #plot(tree.patients)
  #text(tree.patients, pretty=0) 
}





# Read preprocessed patient data
setwd("/Users/Joel/Documents/Research/Cancer Stage Diagnosis/Datasets")
patient_data <- read.csv("anon_patient_expressions_stage_all.csv", header=T, fill=T)

# Separate predictors and response, X and y
patient_stages <- patient_data[,ncol(patient_data)]
patient_data <- patient_data[,-ncol(patient_data)]


# Perform classification using top 10 selected features from file
# ~48% accuracy
#classify_selected(patient_data, patient_stages, "genes_rfe_top_all.txt", 10)

# ~47% accuracy
#classify_selected(patient_data, patient_stages, "genes_tree_top_50.txt", 10)



# This data has a class imbalance problem
# We will test classification between all possible pairs of the 5 classes


# Generate each unique pair of classes / stages (0...4)
for (i in c(0,1,2,3)) {
  for (j in (i+1):4) {
    
    # Determine rows which are stage i or j
    i_indices <- which(patient_stages==i)
    j_indices <- which(patient_stages==j)
    ij_indices <- which(patient_stages==i | patient_stages==j)
    
    # Form new matrix by reducing to i or j rows only
    patient_data_ij = patient_data[ij_indices,]
    patient_stages_ij <- patient_stages[ij_indices]
    
    # Classify this instead
    s <- sprintf("Classifying Stage %d vs. %d, each having %d and %d samples", 
                 i, j, length(i_indices), length(j_indices))
    print(s)
    classify_selected(patient_data_ij, patient_stages_ij, 
                      "genes_rfe_top_all.txt", 20)
  }
}
