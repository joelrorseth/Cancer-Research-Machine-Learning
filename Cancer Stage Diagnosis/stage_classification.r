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

# Read in list of feature selected genes
selected_genes <- readLines("genes_tree_top_1000.txt")

# Separate predictors and response, X and y
patient_stages <- patient_data[,ncol(patient_data)]
patient_data <- patient_data[,-ncol(patient_data)]

new_genetic_data <- patient_data[, colnames(patient_data) %in% selected_genes]


# SECTION 1: Try SVM
library(e1071)
svm_clf <- svm(as.factor(patient_stages)~., data=data.matrix(new_genetic_data), 
               kernel="polynomial", cost=10, scale=FALSE)
plot(svm_clf, data.matrix(new_genetic_data))


# SECTION 2: Random Forest
library(tree)
tree.patients <- tree(as.factor(patient_stages)~., new_genetic_data)

plot(tree.patients)
text(tree.patients, pretty=0)


# TODO: KNN