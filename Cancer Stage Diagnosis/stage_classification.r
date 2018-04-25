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


# SVM classification
classify.svm <- function(X, y, train_indices, filename, top_n) {
  
  # Get feature selected gene names
  selected_genes = readLines(filename)[0:top_n]
  
  # Reduce X to only feature selected columns
  selected.X <- X[, colnames(X) %in% selected_genes]
  
  dat <- data.frame(x=selected.X, y=as.factor(y))
  
  # Train radial svm on training data
  # Create the SVM classifier
  svmfit <- svm(y~., data=dat[train_indices,], kernel="radial", gamma=1, cost=1)
  
  # Use CV to determine optimal gamma
  tune.out <- tune(svm, as.factor(y)~., data=dat[train_indices,], kernel="radial",
                   ranges=list(cost=c(0.1,1,10,100,1000), gamma=c(0.5,1,2,3,4)))
  
  # And use to predict
  #table(true=dat[-train_indices,"y"], 
  #      pred=predict(tune.out$best.model, newdata=dat[-train_indices,]))
  
  svm.pred <- predict(tune.out$best.model, newdata=dat[-train_indices,])
  svm.mean <- mean(svm.pred==dat[-train_indices,"y"])
  
  print(sprintf( "SVM CLF: %f", svm.mean ))
}



# KNN classifier -- this is more or less a baseline classification
classif.knn <- function(train.X, train.y, test.X, test.y, filename, top_n) {
  
  # Get feature selected gene names
  selected_genes = readLines(filename)[0:top_n]
  
  # Reduce train and test X to only feature selected columns
  train.X <- train.X[, colnames(train.X) %in% selected_genes]
  test.X <- test.X[, colnames(test.X) %in% selected_genes]
  
  # Classify
  knn.clf <- knn(train.X, test.X, as.factor(train.y), k=3)
  
  # Mean classification score
  knn.mean <- mean(knn.clf == test.y)
  print(sprintf( "KNN CLF: %f", knn.mean ))
  
  #tab <- table(knn.clf, test.y)
}



# Tree classification -- this should be one of the most promising
classif.tree <- function(X, y, train_indices, filename, top_n) {
  
  # Get feature selected gene names
  selected_genes = readLines(filename)[0:top_n]
  
  # Reduce train and test X to only feature selected columns
  selected.X <- X[, colnames(X) %in% selected_genes]
  
  # Calculate the test sets here to reduce # of parameters
  test.X <- X[-train_indices,]
  test.y <- y[-train_indices]
  
  # Create, and predict using, a tree classifier
  tree.clf <- tree(as.factor(y)~., selected.X, subset=train_indices)
  tree.pred <- predict(tree.clf, test.X, type="class")
  
  # Print mean classification score
  # tab <- table(tree.pred, test.y)
  tree.mean <- (mean(tree.pred==test.y))
  
  print(sprintf( "TREE CLF: %f", tree.mean ))
}




# Read preprocessed patient data
setwd("/Users/Joel/Documents/Research/Breast Cancer Stage Diagnosis Ngom/Cancer Stage Diagnosis/Datasets")
patient_data <- read.csv("anon_patient_expressions_stage_all.csv", header=T, fill=T)

# Separate predictors and response, X and y
X <- patient_data[,-ncol(patient_data)]
y <- patient_data[,ncol(patient_data)]



# This data has a class imbalance problem
# Straight classification w/ feature selection yielded ~47% avg accuracy
# We will now try classification between all possible pairs of the 5 classes

# Generate each unique pair of classes / stages (0...4)
for (i in c(0,1,2,3)) {
  for (j in (i+1):4) {
    
    # Determine rows which are stage i or j
    i_indices <- which(y==i)
    j_indices <- which(y==j)
    ij_indices <- which(y==i | y==j)
    
    # Form new matrix by reducing to i or j rows only
    ij.X = X[ij_indices,]
    ij.y <- y[ij_indices]
    
    # Split i, j data into train (80%) and test (20%) sets
    train_indices = sort(sample(nrow(ij.X), nrow(ij.X)*.8))
    train.X <- ij.X[train_indices,]
    test.X <- ij.X[-train_indices,]
    train.y <- ij.y[train_indices]
    test.y <- ij.y[-train_indices]
    

    s <- sprintf("Classifying Stage %d vs. %d, each having %d and %d samples", 
                 i, j, length(i_indices), length(j_indices))
    print(s)
    print(sprintf( "Train size = %d, Test size = %d", length(train.y), length(test.y) ))
    
    
    if (length(train.y) == 0) next
    if (length(test.y) == 0) next
    
    # Perform classification for this 1v1 classification using different models
    classif.tree(ij.X, ij.y, train_indices, "genes_mRMR_top_50.txt", 7)
    classif.knn(train.X, train.y, test.X, test.y, "genes_mRMR_top_50.txt", 7)
    classify.svm(ij.X, ij.y, train_indices, "genes_mRMR_top_50.txt", 10)
    print(" ")
  }
}
