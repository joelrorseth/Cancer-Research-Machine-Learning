#
# Experiments to address the class imbalance
#
# Joel Rorseth
# August 19, 2018
#

from pandas import read_csv

def read_dataset(filename):
    return read_csv(filename, header=0)

def separate_dataframe(df, i, j):

    # Filter out rows that are not stage 1 or 3 samples, 0/1 encode them
    df = df.drop(df[(df['TUMOR_STAGE']!=i) & (df['TUMOR_STAGE']!=j)].index)
    df['TUMOR_STAGE'] = [1 if stg==i else 0 for stg in df.TUMOR_STAGE]
    df['TUMOR_STAGE'].value_counts()

    y = df.TUMOR_STAGE
    X = df.drop('TUMOR_STAGE', axis=1)
    return X, y


dataframe = read_dataset("../Datasets/anon_patient_expressions_stage_all.csv")
X, y = separate_dataframe(dataframe, 1, 3)
print(y)
