# Task 1: Preprocessing

A brief report of the preprocessing that takes place in `format_dataset.r` and
`format_dataset_anonymous.r`. Many of the original observations from the BRCA dataset are discared
due to several reasons.

## Formatting

The final datasets produced each contain a single header row, most importantly to name each of the
genes listed. There are two near identical formatting scripts, one of which records the patients' IDs in
the first column, and the other (Anon) which leaves this out.


## Pruning Observations

Of the original 2509 patients accounted for in `data_clinical_sample.txt`, many are removed.

### Pass 1 -- Original BRCA breast cancer data
- **Total observations: 2509**

### Pass 2 - Valid patients labelled with stage in `data_clinical_sample.txt`
- Stage 0: 12
- Stage 1: 501
- Stage 2: 825
- Stage 3: 118
- Stage 4: 10
- "null" : 514
- blank : 529
- **Total usable observations: 1466**

### Pass 3 - Patients with data recorded in `data_expression.txt`
- Stage 0: 4
- Stage 1: 475
- Stage 2: 798
- Stage 3: 115
- Stage 4: 9
- **Total recorded: 1401**

## Pruning Genes

Several missing values are present in the raw dataset. The preprocessing script determines that 14
columns (genes) in the matrix have one or more missing values, with some genes simply not having
any value recorded for any patient. Furthermore, the tMDC gene appears twice in the original genetic
dataset. For simplicity, the second column occurrence of this gene is removed. From the 24,375
genes in the original `data_expression.txt`, the following are discarded from the pruned dataset:

- TMPRSS7
- SLC25A19
- CSNK2A1
- BAMBI
- X5.201721816
- X5.140565682
- MRPL24
- X5.908985652
- X5.443175141
- AK127905
- X5.977749751
- C2H2
- X5.369014764
- X5.560773129
- tMDC (the second occurrence of two)

## Output

The final dataset for all cancer stages should be `1402 x 24,361` for the anonymous option, or `1402 x
24,362` when patient ID is included. The anonymous output breaks down as follows:

- 24,374 columns for gene expressions (tMDC gene is recorded twice in expression data, omit second occurrence)
- 1 column for the breast cancer stage (last column)
- N rows for N patients of the selected stage (1401 across all stages)
- 1 row for header (gene name 1, 2 ... 24360, TUMOR_STAGE)
- Matrix size:  `N+1 x 24,361`

### Notes
- When header=T, R will interpret header row and row count will be 1 less
- Right now, 2 samples are generating NA stage, so these are omitted, thus 800 => 798
- These two generated NAs are both for stage 2 samples
