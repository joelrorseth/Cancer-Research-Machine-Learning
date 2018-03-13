#
# Gene Selection
# Experiment with different feature selection algorithms to identify important
# genes in the breast cancer data, remove unimportant genes in a new dataset.
#
# Joel Rorseth
# March 13, 2018
#

# Read preprocessed patient data
setwd("/Users/Joel/Documents/Research/Cancer Stage Diagnosis/Datasets")
patient_data <- read.csv("anon_patient_expressions_stage_all.csv", header=T, fill=T)

# Using the caret library for feature selection
library(mlbench)
library(caret)


# Generate a control object to specify the algorithm behind feature selection
# We will try Random Forest to evaluate each iteration
control <- rfeControl(functions=rfFuncs, method="cv", number=10)


# Run Recursive Feature Elimination (RFE) to find 100 best genes
results <- rfe(patient_data[,-ncol(patient_data)],  # X
               patient_data[,ncol(patient_data)],   # y
               sizes=c(1:100), rfeControl=control)

# Summary
print(results)

# List selected features
predictors(results)

# plot the results
# plot(results, type=c("g", "o"))