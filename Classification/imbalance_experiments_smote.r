#
# Imbalance Experiments: SMOTE
#
# Investigate the performance of the SMOTE strategy, as a means of 
# accomodating the class imbalance in the dataset.
#
# Joel Rorseth
# August 4, 2018
#

library(caret)
library(glmnetUtils)

# SMOTE experiments
balance_smote <- function(train.X, train.y, test.X, test.y) {
  
  training <- train.X
  training$target = train.y
  training <- as.matrix(training)
  
  # TODO: train() is throwing protection stack overflow error
  # TODO: Figure out a way to increase stack size, using --max-ppsize=500000 still caps at 4gb on macOS
  ctrl <- trainControl(method = "cv", number = 5)
  model <- train(target~., data=training, method = "treebag", trControl = ctrl)
}

# Run several classifiers between all possible pairs of classes
one_vs_one_classify <- function(X, y) {
  
  i <- 1
  j <- 3
  
  # Form new matrix by reducing to i or j rows only
  ij_indices <- which(y==i | y==j)
  ij.X = X[ij_indices,]
  ij.y <- y[ij_indices]
  
  # 0/1 encode the y vector
  ij.y <- ifelse((ij.y)==i, 0, 1)

  # Split into test / train  
  split_index <- createDataPartition(as.factor(ij.y), p=0.5, list=F, times=1)
  train.X <- ij.X[split_index,]
  test.X <- ij.X[-split_index,]
  train.y <- ij.y[split_index]
  test.y <- ij.y[-split_index]

  # Pass off to SMOTE
  balance_smote(train.X, train.y, test.X, test.y)
}

# Main program script
driver <- function() {
  
  # Read preprocessed patient data
  setwd("/Users/Joel/Documents/Research/Breast Cancer Stage Diagnosis/Datasets")
  patient_data <- read.csv("anon_patient_expressions_stage_all.csv", header=T, fill=T)
  
  # Separate predictors and response, X and y
  X <- patient_data[,-ncol(patient_data)]
  y <- factor(patient_data[,ncol(patient_data)])
  
  # Almost 90% fall under stage 1 and 2
  table(y)
  prop.table(table(y))
  
  one_vs_one_classify(X, y)
}

driver()
