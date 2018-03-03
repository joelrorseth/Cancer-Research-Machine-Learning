#
# Stage Classification
# Attempt to classify patients' stage of breast cancer given around 
# 24,000 gene expressions and 1400 patient samples.
#
# Joel Rorseth
# February 18, 2018
#


# Read preprocessed patient data
setwd("/Users/Joel/Documents/Research/Cancer Stage Diagnosis/Datasets")
patient_data <- read.csv("anon_patient_expressions_stage_all.csv", header=T, fill=T)


# Determine which columns have NA values
# Try to get a picture of how sparse the dataset is
num_complete_cols <- 0
num_incomplete_cols <- 0

# Check all columns for NA (cols are genes, last col is stage, check that also)
for (col in names(patient_data)) {
  
  num_na <- sum( is.na(patient_data[[col]]) )
  if (num_na > 0) { 
    num_incomplete_cols <- num_incomplete_cols + 1
  } else {
    num_complete_cols <- num_complete_cols + 1
  }
}

sprintf("There were %d complete columns, of %d total", 
        num_complete_cols, length(names(patient_data)))
sprintf("There were %d incomplete columns, of %d total", 
        num_incomplete_cols, length(names(patient_data)))



#patient_stages <- patient_data[,ncol(patient_data)]
#patient_data <- patient_data[,-ncol(patient_data)]

# TODO
# SECTION 1: Try SVM
#library(e1071)
#svm_clf <- svm(patient_stages~., data=data.matrix(patient_data), kernel="polynomial", cost=10, scale=FALSE)

#library(tree)
#tree.patients <- tree(patient_stages~., patient_data)