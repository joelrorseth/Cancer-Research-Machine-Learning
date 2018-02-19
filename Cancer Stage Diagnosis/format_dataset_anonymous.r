#
# Format Dataset (omit patient ID's)
# Script to read and merge two patient gene expression datasets
#
# Joel Rorseth
# January 30, 2018
#



# Change working directory to current dataset location
setwd("/Users/Joel/Documents/Research/Cancer Stage Diagnosis/Datasets/brca_metabric")


###########################################
# Function to convert MB-0002 to MB.0002
###########################################
dotted_name <- function(name) {
  return(gsub("-", ".", name))
}

###########################################
# Function to convert MB.0002 to MB-0002
###########################################
dashed_name <- function(name) {
  return(gsub(".", "-", name))
}


###########################################
# Function to construct the patient gene expression matrix for given cancer stage
###########################################
construct_summary_matrix <- function(stage, patients) {
  
  # Create new list filtered to samples of current stage only
  current_patients <- patients[patients$TUMOR_STAGE==as.character(stage),]
  
  # Get list of current patient id's (dot format)
  current_patients_ids <- dotted_name(current_patients$PATIENT_ID)
  n_current_patients <- length(current_patients_ids)
  
  # Extract expression data for all patients (each a column) with current cancer stage
  # Now have a matrix of 24,375 rows (expressions) by N patients with this stage
  
  current_expressions <- subset(expressions, select=current_patients_ids)
  
  
  # Add row of stages to bottom
  stages_row <- matrix(stage, nrow=1, ncol=n_current_patients)
  colnames(stages_row) <- colnames(current_expressions)
  summary <- rbind.data.frame(current_expressions, stages_row, stringsAsFactors=F)
  
  
  # For anonymous version, we omit patient ids entirely
  # Manually insert PATIENT_IDs as first row in the matrix (embedded column names)
  #new_row_names <- as.character(current_patients$PATIENT_ID)
  #summary <- rbind.data.frame(new_row_names, summary, stringsAsFactors=F)
  
  
  # Transpose
  summary <- t(summary)
  
  # Omit PATIENT_ID from header, this is not present in anonymous version
  #new_column_names <- c("PATIENT_ID", as.character(gene_names_1), "TUMOR_STAGE")
  
  # Manually insert Gene Names as the first row (embedded column names)
  new_column_names <- c(as.character(gene_names_1), "TUMOR_STAGE")
  summary <- rbind.data.frame(new_column_names, summary, stringsAsFactors=F)
  
  
  # Now #patients+1(colnames) x #expressions+1(rownames)+1(stage)
  #print(dim(summary))
  
  return(summary)
}


###########################################
# Write summary matrix of given stage to file
###########################################
write_matrix_to_file <- function(stage, summary) {
  
  # Write to file
  filename <- paste0("anon_patient_expressions_stage_", as.character(stage), ".csv")
  write.table(summary, file=filename, row.names=F, na="NA", col.names=F, sep=",", quote=F)
}





# SECTION 1: Reading in gene expression data
# Read in expression data, remove duplicate gene measurements
expressions <- read.table("data_expression.txt", header=TRUE, fill=TRUE)
expressions <- expressions[!duplicated(expressions$Hugo_Symbol),]

# SECTION 2: Reading in patient cancer stage data
# Read patient dataset, maps patients to breast cancer stage
patients_1 <- read.table("data_clinical_sample.txt", header=TRUE, fill=TRUE)
dim(patients_1) # 2509


# SECTION 3: Cleaning / filtering data
# Toss out any samples with no recorded tumor stage
valid_patient_indices <- !(is.na(patients_1$TUMOR_STAGE) | 
                             patients_1$TUMOR_STAGE=="" |
                             patients_1$TUMOR_STAGE=="null")

patients_2 <- patients_1[valid_patient_indices,]
dim(patients_2) # 1466

# Filter out patients whom there is no recorded gene expression data
patients_3 <- subset(patients_2, dotted_name(PATIENT_ID) %in% colnames(expressions))




# Keep a list of the actual names of genes recorded for patients
gene_names_1 <- expressions$Hugo_Symbol
gene_names_2 <- expressions$Entrez_Gene_Id


summary_all <- matrix()


# SECTION 4: Cross referencing to create patient gene expression matrices
# For each stage of cancer...
for (stage in c(0:4)) {
  
  # Obtain matrix for this stage
  summary <- construct_summary_matrix(stage, patients_3)
  
  print(dim(summary))
  
  # Omit any potential error rows where NA was written
  summary <- summary[!is.na(summary[,ncol(summary)]), ]
  
  print(dim(summary))
  
  # Write it to file
  write_matrix_to_file(stage, summary)
  
  # Add it to the comprehensive summary matrix of all stages
  summary_all <- rbind.data.frame(summary, summary_all[-1,], stringsAsFactors=F)
}


# SECTION 5: Write a sampled matrix containing all patient expression data
# Shuffle the all-stages matrix (remove header, shuffle, reinsert)
header_row <- summary_all[1,]
summary_all <- summary_all[-1,]
summary_all <- summary_all[sample(nrow(summary_all)),]
summary_all <- rbind.data.frame(header_row, summary_all, stringsAsFactors=F)

print(dim(summary_all))
# Omit any potential error rows where NA was written
summary_all <- summary_all[!is.na(summary_all[,ncol(summary_all)]), ]
print(dim(summary_all))

# Write it to file
write_matrix_to_file("all", summary_all)


