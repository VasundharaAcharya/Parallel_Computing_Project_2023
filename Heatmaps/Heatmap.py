
import pandas as pd
import seaborn as sns
import numpy as np
# Load the gene expression data into a Pandas DataFrame. But here we need to merge the cluster assignment column with the rest of the dataset.
gene_data = pd.read_csv('HEATMAP_DATASET_WITH_ASSIGNMENTS.csv')

numeric_cols = gene_data.select_dtypes(include=np.number).columns.tolist()

matrix = gene_data.pivot_table(index='Gene Description', columns='Cluster_assignments', values='ID')


matrix


sns.heatmap(matrix)