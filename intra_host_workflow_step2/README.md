# README for Bash Pipeline Script

## Overview
This script is designed to process `.base` files, and extract intrahost single nucleotide variations (iSNVs) from sequencing data. It filters the iSNVs based on various criteria, including alternative allele frequency (AAF), strand bias likelihood (S), and coverage breath and depth. The pipeline also identifies and removes outlier libraries.

## Usage
```bash
./script.sh -AAF <AAF_threshold> -S <S_threshold> -i <input_file_paths> -o <output_directory>
```

### Arguments:
- `-AAF <AAF_threshold>`: Threshold for minimum Alternative Allele Frequency (AAF).
- `-S <S_threshold>`: Threshold for minimum Strand Bias Likelihood (S).
- `-i <input_file_paths>`: File containing paths to `.base` input files.
- `-o <output_directory>`: Directory where output files will be saved.

### Example:
```bash
./script.sh -AAF 0.05 -S 0.1 -i base_files.txt -o results
```

## Workflow Breakdown

### 1. **Input Argument Check**
   Verifies that exactly 8 input arguments are passed. If not, it displays usage instructions and exits.

### 2. **Parsing `.base` Files**
   Reads the `.base` files and converts them into `.csv` format. The output `.csv` files are stored in the specified output directory under `parsed_csvs`. Errors encountered during parsing are logged.

### 3. **Check for Missing CSV Files**
   Ensures all `.base` files have corresponding `.csv` files. Missing files are logged in `missing_csv_files.txt`.

### 4. **Compute Coverage Metrics**
   For each parsed `.csv` file, the mean depth and breadth of coverage (positions covered by at least 100X depth) are computed and saved to a summary file `coverage_metrics.tsv`.

### 5. **Mask Low Coverage Libraries**
   Filters out libraries that do not meet the coverage thresholds, saving the results to separate files for valid and invalid libraries.

### 6. **Concatenate CSV Files**
   Combines all the `.csv` files into a single `.tsv` file (`unfiltered_iSNVs_all.tsv`). An optional limit can be specified for testing purposes.

### 7. **Filtering and Masking**
   The script applies a series of filters based on the provided AAF and S thresholds, and it masks positions to exclude based on certain criteria.

### 8. **Generate Matrices**
   Creates matrices of iSNVs with 3 columns (position, library ID, AAF) for different filtering stages (unfiltered, de novo, consensus). The matrices are used for downstream analysis.

### 9. **Outlier Detection and Removal**
   Identifies outlier libraries (libraries with abnormal iSNV distributions) and removes them from the final dataset.

### 10. **Annotation**
   Uses a Python script to annotate iSNVs based on metadata, reference genome, and nucleotide fitness data.

### 11. **Output**
   Final output includes:
   - Filtered iSNV files in `.tsv` and matrix formats.
   - Summary of filtering thresholds and counts of valid and invalid libraries.
   - Final annotated iSNV files.

## Parallelization
The script uses GNU Parallel for processing tasks in parallel, improving the efficiency of handling large datasets.

## Error Handling
All steps are logged, and if an error occurs at any point, the script logs the error and exits to avoid continuing with incomplete or incorrect data.

## Dependencies
- `GNU Parallel`
- `awk`, `sed`
- Python environment for running `annotate_iSNVs.py` and other helper scripts.