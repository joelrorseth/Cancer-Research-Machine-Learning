#
# Imbalance Experiments: SMOTE
#
# Investigate the performance of the SMOTE strategy, as a means of 
# accomodating the class imbalance in the dataset.
#
# Joel Rorseth
# August 4, 2018
#

# SMOTE experiments
balance_smote <- function(X, y) {

}

# Run several classifiers between all possible pairs of classes
one_vs_one_classify <- function(X, y, pairings) {
  
  for (p in pairings) {
    i <- p[0]
    j <- p[1]
    
    # Form new matrix by reducing to i or j rows only
    ij_indices <- which(y==i | y==j)
    ij.X = X[ij_indices,]
    ij.y <- y[ij_indices]
    
    # Pass off to SMOTE
    balance_smote(ij.X, ij.y)
  }
}

# Main program script
driver <- function() {
  
  # Read preprocessed patient data
  setwd("/Users/Joel/Documents/Research/Breast Cancer Stage Diagnosis/Datasets")
  patient_data <- read.csv("anon_patient_expressions_stage_all.csv", header=T, fill=T)
  
  # Separate predictors and response, X and y
  X <- patient_data[,-ncol(patient_data)]
  y <- factor(patient_data[,ncol(patient_data)])
  
  # Almost 90% fall under stage 2 and 3
  table(y)
  prop.table(table(y))
  
  pairings <- mapply(c, c(0,1), c(1,2), SIMPLIFY=F)
  one_vs_one_classify(X, y, pairings)
}

driver()