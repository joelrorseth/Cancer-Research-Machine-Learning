#
# Gene Selection
# Use Python feature selection algorithms to attempt to identify important
# genes, and remove unimportant.
#

import pymrmr
from pandas import read_csv
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFE


# Read the patient genetic dataset
def read_dataset(filename):
    return read_csv(filename, header=0)


# Dataset formatting to convert DataFrame to numpy arrays X and y
def separate_dataset(df):

    _, columns = df.shape
    numpy_dataset = df.values

    colnames = df.columns.values[:-1]   # Omit TUMOR_STAGE
    X = numpy_dataset[:, :(columns-1)]  # Omit stage in last col
    y = numpy_dataset[:, columns-1]

    return X,y,colnames


# Write list of (feature selected) genes to txt file
def write_genes(genes, algo):

    genes_file = open('genes_' + algo + '_top_' + str(len(genes)) + '.txt', 'w')

    for gene in genes:
        genes_file.write("%s\n" % gene)


# Run the mRMR algorithm using a DataFrame
# TODO: In progress, waiting on installation bug
def select_n_genes_mRMR(df, num_genes):
    return pymrmr.mRMR(df, 'MIQ', num_genes)


# Find top n genes, produce files using all selection algorithms for several n
def find_top_n_genes(df, X, y, gene_names):

    ranked_tree = select_genes_tree(X, y, gene_names)

    # Write each result to file
    for n in [10, 20, 50, 100, 200, 500, 1000]:
        write_genes(ranked_tree[0:n], "tree")

        ranked_mRMR = select_n_genes_mRMR(df, n)
        write(ranked_mRMR, "mRMR")

        # NOTE: Taking too long to run
        #ranked_rfe = select_n_genes_rfe(X, y, gene_names, n)
        #write(ranked_rfe, "rfe")


# Run RFE selection using logistic regression, determine list ranking genes
def select_n_genes_rfe(X, y, gene_names, n):

    # Run RFE to f
    rfe = RFE( LogisticRegression(), n )
    fit = rfe.fit(X, y)

    return [name for i, name in enumerate(gene_names) if (fit.support_)[i]==True]


# Run decision tree selection algorithm, determine list ranking genes
def select_genes_tree(X, y, gene_names):

    # Use bagged decision tree to decide importance
    model = ExtraTreesClassifier()
    model.fit(X, y)

    importances = model.feature_importances_

    # Zip together gene names with their importances, sort by importance
    ranked = zip(gene_names, importances)
    ranked = sorted(ranked, key=lambda pair: pair[1], reverse=True)

    # Return the top 'num_genes' gene names, sorted by importance
    return [name for name, _ in ranked]


def main():

    dataset = read_dataset("./Datasets/anon_patient_expressions_stage_all.csv")
    X, y, gene_names = separate_dataset(dataset)

    # Try selecting the top n genes, write to separate files
    find_top_n_genes(dataset, X, y, gene_names)

main()
