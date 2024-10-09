import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns

class OutlierLibraryIdentifier:
    def __init__(self, input_file, output_directory):
        self.input_file = input_file
        self.output_directory = os.path.join(output_directory, 'outlier_libraries_analysis')
        os.makedirs(self.output_directory, exist_ok=True)
    
    def load_data(self):
        self.data = pd.read_csv(self.input_file, sep='\t')
    
    def identify_outliers(self):
        # Group by library_id and count the number of unique position_ids
        self.library_counts = self.data.groupby('library_id')['position_id'].nunique()
        
        # Calculate the threshold for the top 1% libraries
        self.threshold = self.library_counts.quantile(0.99)
        
        # Identify outlier libraries
        self.outliers = self.library_counts[self.library_counts > self.threshold]
    
    def save_outliers(self):
        # Save the outlier libraries to a file
        outliers_file = os.path.join(self.output_directory, 'outlier_libraries.txt')
        self.outliers.index.to_series().to_csv(outliers_file, sep='\t', header=False, index=False)
        return outliers_file
    
    def save_distribution_data(self):
        # Save the distribution data to a TSV file
        distribution_data_file = os.path.join(self.output_directory, 'library_distribution_data.tsv')
        self.library_counts.to_csv(distribution_data_file, sep='\t', header=True)
        return distribution_data_file
    
    def plot_distribution(self):
        sns.set(rc={'figure.figsize':(10,10)}, style="ticks", context='talk', font_scale=1.5)
        fig, ax = plt.subplots()
        
        sns.histplot(self.library_counts, bins=30, kde=False, label='Library Counts', alpha=0.5, ax=ax)
        ax.axvline(self.threshold, color='red', linestyle='dashed', linewidth=1, label='99th Percentile Threshold')
        
        ax.legend(loc='upper right')
        ax.set_xlabel('Number of Genomic Positions')
        ax.set_ylabel('Number of Libraries')
        ax.set_title('Distribution of Genomic Positions per Library')
        ax.grid(True)
        
        plot_file = os.path.join(self.output_directory, 'library_distribution.png')
        plt.savefig(plot_file)
        plt.close()
        return plot_file
    
    def generate_all(self):
        self.load_data()
        self.identify_outliers()
        outliers_file = self.save_outliers()
        distribution_data_file = self.save_distribution_data()
        plot_file = self.plot_distribution()
        return outliers_file, distribution_data_file, plot_file

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_directory = sys.argv[2]
    
    outlier_identifier = OutlierLibraryIdentifier(input_file, output_directory)
    outliers_file, distribution_data_file, plot_file = outlier_identifier.generate_all()
    print(f"Outlier libraries saved to: {outliers_file}")
    print(f"Distribution data saved to: {distribution_data_file}")
    print(f"Distribution plot saved to: {plot_file}")