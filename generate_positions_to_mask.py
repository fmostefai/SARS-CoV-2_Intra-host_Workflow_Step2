import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy import stats

class PositionMaskGenerator:
    def __init__(self, S_threshold, input_file, output_directory):
        self.S_threshold = S_threshold
        self.input_file = input_file
        self.output_directory = output_directory
        # Create a specific results folder within the output directory
        self.results_folder = os.path.join(output_directory, f'masking_results_S_threshold_{S_threshold}')
        os.makedirs(self.results_folder, exist_ok=True)  # Create results folder if it does not exist
    
    def load_data(self):
        # Load data from the input file into a pandas DataFrame
        self.data = pd.read_csv(self.input_file, sep='\t')
    
    def apply_filter_label(self):
        # Check if 'strand_bias_likelihood' column exists
        if 'strand_bias_likelihood' not in self.data.columns:
            raise KeyError("'strand_bias_likelihood' column not found in input file.")
        # Apply filter label based on the S_threshold
        self.data['filter_label'] = np.where(self.data['strand_bias_likelihood'] >= self.S_threshold, 'kept', 'removed')
    
    def generate_permutations(self):
        # Generate permutations of the filter label for further analysis
        for i in range(1, 10):
            column_name = f'permutation_n{i}'
            self.data[column_name] = np.random.permutation(self.data['filter_label'])
    
    def create_position_counts(self):
        # Create position counts based on the filter labels
        filter_counts = self.data.groupby('position_id')['filter_label'].value_counts().unstack(fill_value=0)
        self.position_counts = filter_counts
        
        # Join permutation counts to the main position counts
        for i in range(1, 10):
            permutation_counts = self.data.groupby('position_id')[f'permutation_n{i}'].value_counts().unstack(fill_value=0)
            permutation_counts.columns = [f'permutation_n{i}_{col}' for col in permutation_counts.columns]
            self.position_counts = self.position_counts.join(permutation_counts, how='outer').fillna(0).astype(int)
    
    def calculate_library_count_threshold(self):
        # Calculate the library count threshold based on the shuffled data
        removed_counts = self.position_counts['removed']
        shuffled_removed_counts = self.position_counts['permutation_n1_removed']
        library_count_threshold = shuffled_removed_counts.quantile(0.99)
        self.library_count_threshold = round(library_count_threshold)
    
    def generate_positions_to_mask(self):
        # Generate positions to mask based on the library count threshold
        genomic_positions_to_mask = self.position_counts[self.position_counts['removed'] > self.library_count_threshold].index
        positions_to_mask_file = os.path.join(self.results_folder, "genomic_positions_to_mask.txt")
        np.savetxt(positions_to_mask_file, genomic_positions_to_mask, fmt='%d')
        return positions_to_mask_file
    
    def generate_distribution_plots(self):
        # Generate distribution plots of the true and shuffled data
        sns.set(rc={'figure.figsize':(10,10)}, style="ticks", context='talk', font_scale=1.5)
        fig, ax = plt.subplots()
        
        true_counts = self.position_counts['removed']
        shuffled_counts = self.position_counts['permutation_n1_removed']
        
        sns.histplot(true_counts, bins=20, kde=False, label='True Distribution', alpha=0.5, ax=ax)
        sns.histplot(shuffled_counts, bins=20, kde=False, label='Shuffled Distribution', color='#dd8453', alpha=0.5, ax=ax)
        
        ax.legend(loc='upper right')
        ax.set_xlabel('Removed iSNV Count')
        ax.set_ylabel('Frequency')
        ax.set_title('Distribution of Removed iSNV Counts')
        
        distribution_plot_file = os.path.join(self.results_folder, "distribution_plot.png")
        plt.savefig(distribution_plot_file)
        
        return distribution_plot_file
    
    def save_distribution_data(self):
        # Save the distribution data to a TSV file
        distribution_data_file = os.path.join(self.results_folder, 'distribution_data.tsv')
        distribution_data = self.position_counts[['removed', 'permutation_n1_removed']].reset_index()
        distribution_data.columns = ['position_id', 'removed', 'permutation_n1_removed']
        distribution_data.to_csv(distribution_data_file, sep='\t', index=False)
        return distribution_data_file
    
    def save_distribution_info(self):
        # Save distribution information to a text file
        distribution_info_file = os.path.join(self.results_folder, "distribution_info.txt")
        with open(distribution_info_file, 'w') as f:
            f.write("Distribution Information:\n")
            f.write(f"S Threshold: {self.S_threshold}\n")
            f.write(f"Library Count Threshold: {self.library_count_threshold}\n")
            f.write("\nTrue Distribution:\n")
            f.write(f"Mean: {self.position_counts['removed'].mean()}\n")
            f.write(f"Median: {self.position_counts['removed'].median()}\n")
            f.write(f"75th Percentile: {np.percentile(self.position_counts['removed'], 75)}\n")
            f.write(f"99th Percentile: {np.percentile(self.position_counts['removed'], 99)}\n")
            f.write(f"Maximum: {self.position_counts['removed'].max()}\n")
            f.write("\nShuffled Distribution:\n")
            f.write(f"Mean: {self.position_counts['permutation_n1_removed'].mean()}\n")
            f.write(f"Median: {self.position_counts['permutation_n1_removed'].median()}\n")
            f.write(f"75th Percentile: {np.percentile(self.position_counts['permutation_n1_removed'], 75)}\n")
            f.write(f"99th Percentile: {np.percentile(self.position_counts['permutation_n1_removed'], 99)}\n")
            f.write(f"Maximum: {self.position_counts['permutation_n1_removed'].max()}\n")
        
        return distribution_info_file

    def generate_all(self):
        # Execute all the steps to generate required files and plots
        self.load_data()
        try:
            self.apply_filter_label()
        except KeyError as e:
            print("Error: Column names in the input file:", self.data.columns)
            raise e
        self.generate_permutations()
        self.create_position_counts()
        self.calculate_library_count_threshold()
        positions_to_mask_file = self.generate_positions_to_mask()
        distribution_plot_file = self.generate_distribution_plots()
        distribution_data_file = self.save_distribution_data()
        distribution_info_file = self.save_distribution_info()
        return positions_to_mask_file, distribution_plot_file, distribution_data_file, distribution_info_file

if __name__ == "__main__":
    import sys
    
    # Parse command line arguments
    S_threshold = float(sys.argv[1])
    input_file = sys.argv[2]
    output_directory = sys.argv[3]
    
    # Instantiate PositionMaskGenerator
    mask_generator = PositionMaskGenerator(S_threshold, input_file, output_directory)
    
    # Generate all required files and plots
    positions_to_mask_file, distribution_plot_file, distribution_data_file, distribution_info_file = mask_generator.generate_all()
    
    print(f"Generated positions to mask file: {positions_to_mask_file}")
    print(f"Generated distribution plot file: {distribution_plot_file}")
    print(f"Generated distribution data file: {distribution_data_file}")
    print(f"Generated distribution info file: {distribution_info_file}")