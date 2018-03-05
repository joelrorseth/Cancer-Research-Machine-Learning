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

# Column missing value filtering moved to preprocessing scripts


#patient_stages <- patient_data[,ncol(patient_data)]
#patient_data <- patient_data[,-ncol(patient_data)]

# TODO
# SECTION 1: Try SVM
#library(e1071)
#svm_clf <- svm(patient_stages~., data=data.matrix(patient_data), kernel="polynomial", cost=10, scale=FALSE)

#library(tree)
#tree.patients <- tree(patient_stages~., patient_data)