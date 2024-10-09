import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
import os
import sys
import warnings

class PNNComputer:
    def __init__(self, embedding_file, metadata_file, output_directory, n_neighbors=100):
        self.embedding_file = embedding_file
        self.metadata_file = metadata_file
        self.output_directory = output_directory
        self.n_neighbors = n_neighbors
        self.load_data()

    def load_data(self):
        # Load embedding data with header and rename the first column to library_id
        self.df = pd.read_csv(self.embedding_file, sep='\t')
        self.df.rename(columns={self.df.columns[0]: 'library_id'}, inplace=True)
        
        # Load metadata
        self.metadata = pd.read_csv(self.metadata_file, sep='\t')
        
        # Merge metadata with embeddings data on library_id
        self.df = self.df.merge(self.metadata, on='library_id', how='left')

        # Ensure the embedding columns are numeric
        self.df.iloc[:, 1:3] = self.df.iloc[:, 1:3].apply(pd.to_numeric)

    def compute_PNN_for_row(self, row, knn, y, label):
        # Silence the FutureWarning
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            
            # Get the neighbors for the current row
            distances, indices = knn.kneighbors([row[1:3]])
            neighbor_labels = y[indices[0][1:]]
            
            # Calculate the percentage of true labels among the neighbors
            true_label_count = (neighbor_labels == row[label]).sum()
            
        return (true_label_count / self.n_neighbors) * 100

    def compute_PNN(self, label, permute=False):
        """
        Computes the Percentage of Nearest Neighbors (PNN) for a given label.

        Parameters:
        - label (str): The name of the label column in the dataframe to use for PNN calculation.
        - permute (bool): If True, the labels are permuted (randomized) to serve as a baseline.

        Returns:
        - PNN_values (pd.Series): A series containing the PNN values for each row in the dataframe.
        """
        
        # Select the embeddings columns (assuming columns 2 and 3 are the embeddings)
        X = self.df.iloc[:, 1:3].values
        
        # If permute is True, permute the labels
        if permute:
            y = np.random.permutation(self.df[label])
        else:
            y = self.df[label]
        
        # Initialize the Nearest Neighbors model with n_neighbors + 1 (to account for the point itself)
        knn = NearestNeighbors(n_neighbors=self.n_neighbors + 1)
        
        # Fit the KNN model with the embeddings data
        knn.fit(X)
        
        # Compute the PNN values for each row
        PNN_values = self.df.apply(lambda row: self.compute_PNN_for_row(row, knn, y, label), axis=1)
        return PNN_values

    def compute_all_PNN(self):
        # Compute PNN for seq_center_label and its baseline
        self.df['PNN_SC'] = self.compute_PNN('seq_center_label')
        self.df['PNN_SC_baseline'] = self.compute_PNN('seq_center_label', permute=True)
        # Compute PNN for WHO_label and its baseline
        self.df['PNN_lineage'] = self.compute_PNN('WHO_label')
        self.df['PNN_lineage_baseline'] = self.compute_PNN('WHO_label', permute=True)
        # Compute PNN for month and its baseline
        self.df['PNN_month'] = self.compute_PNN('month')
        self.df['PNN_month_baseline'] = self.compute_PNN('month', permute=True)
        
        # Ensure the seq_center_label and WHO_label columns are included in the final output
        self.df['seq_center_label'] = self.df['seq_center_label']
        self.df['WHO_label'] = self.df['WHO_label']
        # Ensure the seq_center_label and WHO_label columns are included in the final output
        self.df['seq_center_label'] = self.df['seq_center_label']
        self.df['WHO_label'] = self.df['WHO_label']

    def save_results(self):
        # Define the output file path
        output_file = os.path.join(self.output_directory, os.path.basename(self.embedding_file).replace('.tsv', '_PNN.tsv'))
        # Reorder columns to match the desired output
        columns_order = ['library_id', self.df.columns[1], self.df.columns[2], 'PNN_SC', 'PNN_SC_baseline', 'PNN_lineage', 'PNN_lineage_baseline', 'WHO_label', 'seq_center_label', 'seq_center_label_acronym', 'month']
        self.df = self.df[columns_order]
        # Save the DataFrame with the computed PNN values
        self.df.to_csv(output_file, sep='\t', index=False)

    def run(self):
        self.compute_all_PNN()  # Compute all PNN values
        self.save_results()  # Save the results

if __name__ == "__main__":
    embedding_file = sys.argv[1]
    metadata_file = sys.argv[2]
    output_directory = sys.argv[3]

    pnn_computer = PNNComputer(embedding_file, metadata_file, output_directory)
    pnn_computer.run()
    print(f"PNN computation completed. Results saved to {output_directory}")