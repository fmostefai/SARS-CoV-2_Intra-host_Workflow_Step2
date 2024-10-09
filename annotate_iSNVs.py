import pandas as pd
import numpy as np
import os
import sys
import re

class ISNVAnnotator:
    def __init__(self, filtered_iSNVs_file, reference_file, metadata_file, nt_fitness_file, output_directory):
        self.filtered_iSNVs_file = filtered_iSNVs_file
        self.reference_file = reference_file
        self.metadata_file = metadata_file
        self.nt_fitness_file = nt_fitness_file
        self.output_directory = output_directory
        os.makedirs(self.output_directory, exist_ok=True)
        self.load_data()
    
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

        self.filtered_iSNVs = pd.read_csv(self.filtered_iSNVs_file, sep='\t', dtype=str, low_memory=False)
        
        for col, dtype in dtype_spec.items():
            if dtype == float:
                self.filtered_iSNVs[col] = self.filtered_iSNVs[col].replace('None', np.nan).astype(dtype)
            else:
                self.filtered_iSNVs[col] = self.filtered_iSNVs[col].astype(dtype)

        self.reference = pd.read_csv(self.reference_file, sep='\t', index_col='position_id')
        self.metadata = pd.read_csv(self.metadata_file, sep='\t', dtype={'library_id': str, 'date': str, 'WHO_label': str}, low_memory=False)
        self.metadata['date'] = pd.to_datetime(self.metadata['date'], format='%Y-%m-%d')
        
        if os.path.exists(self.nt_fitness_file):
            self.nt_fitness = pd.read_csv(self.nt_fitness_file, sep=',')
        else:
            self.nt_fitness = pd.DataFrame(columns=['nt_site', 'nt', 'fitness'])

    def annotate_isnvs(self):
        iSNVs = self.filtered_iSNVs.copy()
        
        # Mutation
        iSNVs['mutation'] = iSNVs.apply(lambda x: self.reference.loc[x['position_id'], 'ref'] + '>' + x['alternative_allele'] if x['position_id'] in self.reference.index else 'unknown', axis=1)
        
        # Consequence
        iSNVs['consequence'] = iSNVs.apply(lambda x: self.reference.loc[x['position_id'], f'{x["alternative_allele"]}_consequence'] if x['position_id'] in self.reference.index else 'unknown', axis=1)
        
        # Gene
        iSNVs['gene'] = iSNVs.apply(lambda x: self.reference.loc[x['position_id'], 'gene'] if x['position_id'] in self.reference.index else 'unknown', axis=1)

        # Amino Acid Change
        iSNVs['aa'] = iSNVs.apply(lambda x: self.reference.loc[x['position_id'], f'{x["alternative_allele"]}_aa'] if x['position_id'] in self.reference.index else 'unknown', axis=1)
        
        amino_acid_map = {
            'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
            'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Gly': 'G',
            'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
            'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
            'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
        }

        def convert_format(s):
            if pd.isna(s) or not isinstance(s, str):
                return s
            match = re.match(r'p\.([A-Za-z]+)(\d+)([A-Za-z]+)', s)
            if match:
                return amino_acid_map.get(match.group(1), 'X') + match.group(2) + amino_acid_map.get(match.group(3), 'X')
            return s
        
        iSNVs['aa_letter'] = iSNVs['aa'].apply(convert_format)
        iSNVs['amino_acid_change'] = iSNVs['gene'] + ':' + iSNVs['aa_letter']

        # Merge with metadata
        iSNVs = pd.merge(iSNVs, self.metadata[['library_id', 'date', 'month', 'WHO_label', 'seq_center_label']], on='library_id', how='left')
        
        # Merge with nt_fitness to add fitness values
        if not self.nt_fitness.empty:
            self.nt_fitness['nt_site'] = self.nt_fitness['nt_site'].astype(int)
            iSNVs['position_id'] = iSNVs['position_id'].astype(int)
            annotated_iSNVs = pd.merge(iSNVs, self.nt_fitness[['nt_site', 'nt', 'fitness']], left_on=['position_id', 'alternative_allele'], right_on=['nt_site', 'nt'], how='left')
        else:
            iSNVs['fitness'] = 'unknown'
            annotated_iSNVs = iSNVs

        return annotated_iSNVs

    def save_annotated_isnvs(self, iSNVs, filename):
        output_file = os.path.join(self.output_directory, filename)
        iSNVs.to_csv(output_file, sep='\t', index=False)
        return output_file
    
    def generate_all(self):
        annotated_iSNVs = self.annotate_isnvs()
        filename = os.path.basename(self.filtered_iSNVs_file).replace('.tsv', '_annotated.tsv')
        self.save_annotated_isnvs(annotated_iSNVs, filename)
        print(f"Annotated iSNVs saved to: {os.path.join(self.output_directory, filename)}")

if __name__ == "__main__":
    filtered_iSNVs_file = sys.argv[1]
    reference_file = sys.argv[2]
    metadata_file = sys.argv[3]
    nt_fitness_file = sys.argv[4]
    output_directory = sys.argv[5]
    
    annotator = ISNVAnnotator(filtered_iSNVs_file, reference_file, metadata_file, nt_fitness_file, output_directory)
    annotator.generate_all()