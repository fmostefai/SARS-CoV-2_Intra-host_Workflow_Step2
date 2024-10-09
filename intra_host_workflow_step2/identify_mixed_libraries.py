import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

class MixedLibraryIdentifier:
    def __init__(self, filtered_consensus_iSNVs_file, filtered_denovo_iSNVs_file, metadata_file, output_directory):
        self.filtered_consensus_iSNVs_file = filtered_consensus_iSNVs_file
        self.filtered_denovo_iSNVs_file = filtered_denovo_iSNVs_file
        self.metadata_file = metadata_file
        self.output_directory = os.path.join(output_directory, 'mixed_libraries_analysis')
        os.makedirs(self.output_directory, exist_ok=True)
    
    def load_data(self):
        # Specify dtype to avoid mixed types warning
        dtype_spec = {
            'library_id': str,
            'position_id': int,
            'num_alleles_observed': int,
            'A_count': int,
            'C_count': int,
            'G_count': int,
            'T_count': int,
            'allele_depth_count': int,
            'forward_strand_count': int,
            'forward_strand_ratio': float,
            'reverse_strand_count': int,
            'reverse_strand_ratio': float,
            'alternative_allele': str,
            'alternative_allele_count': int,
            'alternative_allele_forward_strand_count': int,
            'alternative_allele_forward_strand_ratio': float,
            'alternative_allele_reverse_strand_count': int,
            'alternative_allele_reverse_strand_ratio': float,
            'alternative_allele_frequency': float,
            'strand_bias_likelihood': float
        }
        self.filtered_consensus_iSNVs = pd.read_csv(self.filtered_consensus_iSNVs_file, sep='\t', dtype=str, low_memory=False)
        self.filtered_denovo_iSNVs = pd.read_csv(self.filtered_denovo_iSNVs_file, sep='\t', dtype=str, low_memory=False)

        for col, dtype in dtype_spec.items():
            if dtype == float:
                self.filtered_consensus_iSNVs[col] = self.filtered_consensus_iSNVs[col].replace('None', np.nan).astype(dtype)
                self.filtered_denovo_iSNVs[col] = self.filtered_denovo_iSNVs[col].replace('None', np.nan).astype(dtype)
            else:
                self.filtered_consensus_iSNVs[col] = self.filtered_consensus_iSNVs[col].astype(dtype)
                self.filtered_denovo_iSNVs[col] = self.filtered_denovo_iSNVs[col].astype(dtype)

        self.metadata = pd.read_csv(self.metadata_file, sep='\t', dtype={'library_id': str}, low_memory=False)
        self.metadata['month'] = pd.to_datetime(self.metadata['month'], format='%Y-%m')
    
    def calculate_position_freq(self, threshold=0.1):
        self.filtered_consensus_iSNVs = self.filtered_consensus_iSNVs.merge(self.metadata[['library_id', 'month']], on='library_id')
        position_counts = self.filtered_consensus_iSNVs.groupby(['month', 'position_id']).library_id.count().unstack().fillna(0)
        library_counts = self.filtered_consensus_iSNVs.groupby(['month']).library_id.nunique().to_frame().rename(columns={'library_id': 'library_count'})
        self.position_freq = position_counts.div(library_counts['library_count'], axis=0)
        
        # Save the position frequency table
        self.position_freq.to_csv(os.path.join(self.output_directory, 'position_frequency_over_time.tsv'), sep='\t')
        
        # Plot the frequency over time
        sns.set(rc={'figure.figsize':(16,10)}, style="ticks", context='notebook', font_scale=1.5)
        for position in self.position_freq.columns:
            plt.plot(self.position_freq.index, self.position_freq[position], c='grey')
        plt.title('Frequency of Genomic Positions for Each Month \n')
        plt.xlabel('\n Month')
        plt.ylabel('Derived Allele Frequency \n')
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_directory, 'position_frequency_over_time.png'))
        plt.close()

    def select_positions_based_on_freq(self, threshold=0.1):
        selected_positions = self.position_freq.loc[:, (self.position_freq > threshold).any()]
        self.selected_positions_list = list(selected_positions.columns)
        
        # Save selected positions
        with open(os.path.join(self.output_directory, 'selected_positions.txt'), "w") as f:
            for position in self.selected_positions_list:
                f.write(f"{position}\n")

    def find_possible_mixed_libraries(self, number_of_consensus_iSNVs, max_VAF=0.75):
        possible_mixed = self.filtered_denovo_iSNVs[(self.filtered_denovo_iSNVs.alternative_allele_frequency <= max_VAF) & (self.filtered_denovo_iSNVs.position_id.isin(self.selected_positions_list))]
        possible_mixed = possible_mixed.groupby('library_id').filter(lambda x: len(x) >= number_of_consensus_iSNVs)
        self.possible_mixed_libraries = list(possible_mixed.library_id.unique())
        
        # Save possible mixed libraries
        with open(os.path.join(self.output_directory, 'possible_mixed_libraries.txt'), "w") as f:
            for lib in self.possible_mixed_libraries:
                f.write(f"{lib}\n")

    def generate_all(self):
        self.load_data()
        self.calculate_position_freq()
        self.select_positions_based_on_freq()
        self.find_possible_mixed_libraries(3)
        print(f"Position frequency table and plot saved to: {self.output_directory}")
        print(f"Selected positions saved to: {self.output_directory}/selected_positions.txt")
        print(f"Possible mixed libraries saved to: {self.output_directory}/possible_mixed_libraries.txt")

if __name__ == "__main__":
    filtered_consensus_iSNVs_file = sys.argv[1]
    filtered_denovo_iSNVs_file = sys.argv[2]
    metadata_file = sys.argv[3]
    output_directory = sys.argv[4]
    
    mixed_library_identifier = MixedLibraryIdentifier(filtered_consensus_iSNVs_file, filtered_denovo_iSNVs_file, metadata_file, output_directory)
    mixed_library_identifier.generate_all()
    print(f"Analysis results saved to: {output_directory}")
