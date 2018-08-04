#
# Gene Selection
# Experiment with different feature selection algorithms to identify important
# genes in the breast cancer data, remove unimportant genes in a new dataset.
#
# Joel Rorseth
# March 13, 2018
#

# Read preprocessed patient data
setwd("/Users/Joel/Documents/Research/Breast Cancer Stage Diagnosis/Datasets")
patient_data <- read.csv("anon_patient_expressions_stage_all.csv", header=T, fill=T)

# Using the caret library for feature selection
library(mlbench)
library(caret)
library(mRMRe)


# Perform Recursive Feature Elimination (RFE) feature selection
selection.RFE <- function() {
  
  # Generate a control object to specify the algorithm behind feature selection
  # We will try Random Forest to evaluate each iteration
  control <- rfeControl(functions=rfFuncs, method="cv", number=10)


  # Run Recursive Feature Elimination (RFE) to find 100 best genes
  results <- rfe(patient_data[,-ncol(patient_data)],  # X
               patient_data[,ncol(patient_data)],   # y
               sizes=c(1:100), rfeControl=control)
  # Summary
  #print(results)
  
  # List selected features, write them to file
  predictors(results)
  write(predictors(results), "genes_rfe_top_all.txt", sep="\n")
}



# Perform Minimum Redundancy Maximum Relevance (mRMR) feature selection
selection.mRMR <- function(X, top_n) {
  
  data <- mRMR.data(data = patient_data)
  selected <- mRMR.ensemble(data = data, target_indices = 1, 
                            feature_count = top_n, solution_count = 1)
  
  # Using column indices selected by mRMR, get the gene (column) names
  selected_col_indices <- unlist(solutions(selected))
  gene_names <- colnames(patient_data)[selected_col_indices]
  
  # Write colnames for each index selected to file
  filename <- sprintf("genes_mRMR_top_%d.txt", top_n)
  write(gene_names, filename, sep="\n")  
}



# Run the selection process
selection.mRMR(patient_data, 50)