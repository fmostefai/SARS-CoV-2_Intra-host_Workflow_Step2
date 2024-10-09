#!/bin/bash

# Function to check if input arguments are provided
check_input_arguments() {
    if [ "$#" -ne 8 ]; then
        echo "Usage: $0 -AAF <AAF_threshold> -S <S_threshold> -i <input_file_paths> -o <output_directory>" | tee -a "$log_file"
        exit 1
    fi
}

# Function to parse .base files
parse_base_files() {
    # Check if parsed_csvs directory already exists
    if [ -d "${output_directory}/parsed_csvs" ]; then
        echo "parsed_csvs directory already exists. Skipping parsing step." | tee -a "$log_file"
        return 0
    fi

    mkdir -p "${output_directory}/parsed_csvs"
    mkdir -p "${output_directory}/logs"

    # Run the parsing step in parallel and log errors
    cat "${input_file_paths}" | time parallel "python ${script_directory}/parse_basefile.py {} ${output_directory}/parsed_csvs" 2> "${output_directory}/logs/parse_errors.log" | tee -a "$log_file"
    if [ $? -ne 0 ]; then
        echo "Error parsing .base files" | tee -a "$log_file"
        exit 1
    fi
}

# Function to check for missing .csv files
check_for_missing_csv_files() {
    mkdir -p "${output_directory}/logs"

    # Create a file to store the list of missing .csv files
    missing_file="${output_directory}/logs/missing_csv_files.txt"
    touch "$missing_file"

    # Read all .base files and store their names without the .base extension in a temporary file
    awk -F '/' '{print $NF}' "${input_file_paths}" | sed 's/\.base$//' > "${output_directory}/logs/all_base_files.txt"

    # Read all generated .csv files and store their names without the .csv extension in a temporary file
    find "${output_directory}/parsed_csvs" -name "*.csv" -exec basename {} .csv \; > "${output_directory}/logs/generated_csv_files.txt"

    # Use comm to find the .base files that don't have corresponding .csv files
    comm -23 <(sort "${output_directory}/logs/all_base_files.txt") <(sort "${output_directory}/logs/generated_csv_files.txt") |
    while IFS= read -r base_file_name; do
        # Reconstruct the full path of the missing .base file
        grep "/${base_file_name}.base$" "${input_file_paths}" >> "$missing_file"
    done

    echo "List of .base files without corresponding .csv files saved to $missing_file" | tee -a "$log_file"
}

# Function to compute the mean depth and breath of coverage for each library
compute_coverage_metrics() {
    # Define the path to the metrics summary file and add the header
    metrics_file="${output_directory}/coverage_metrics.tsv"
    echo -e "library_id\tmean_depth\tbreath_100X" > "$metrics_file"

    # Define an AWK script to calculate mean_depth and breath_100X
    awk_script='
    BEGIN { FS = "," }
    NR > 1 {
        sum += $8
        if ($8 > 100) count++
    }
    END {
        mean_depth = sum / 29903
        print mean_depth, count
    }
    '

    export output_directory metrics_file awk_script

    # Loop through each CSV file in the parsed CSVs directory using GNU Parallel
    find "${output_directory}/parsed_csvs" -name '*.csv' | parallel '
        csv_file={}
        library_id=$(basename "$csv_file" .csv | sed "s/\.base$//")
        read mean_depth breath_100X < <(awk "$awk_script" "$csv_file")
        echo -e "${library_id}\t${mean_depth}\t${breath_100X}" >> "$metrics_file"
    '

    # Check if the previous command executed successfully
    if [ $? -ne 0 ]; then
        echo "Error computing coverage metrics" | tee -a "$log_file"
        exit 1
    fi
}

# Application of the filters
mask_low_coverage_libraries() {
    # Read the coverage metrics file and filter libraries based on the criteria
    metrics_file="${output_directory}/coverage_metrics.tsv"
    valid_metrics_file="${output_directory}/valid_coverage_metrics.tsv"
    invalid_metrics_file="${output_directory}/invalid_coverage_metrics.tsv"

    # Filter out libraries that do not meet the coverage criteria
    awk -F '\t' 'NR==1 || ($2 > 100 && $3 > 20000)' "$metrics_file" > "$valid_metrics_file"
    awk -F '\t' 'NR==1 || !($2 > 100 && $3 > 20000)' "$metrics_file" > "$invalid_metrics_file"
    
    if [ $? -ne 0 ]; then
        echo "Error filtering coverage metrics" | tee -a "$log_file"
        exit 1
    fi

    # Create a list of valid and invalid libraries
    valid_libraries=$(awk -F '\t' 'NR>1 {print $1}' "$valid_metrics_file")
    invalid_libraries=$(awk -F '\t' 'NR>1 {print $1}' "$invalid_metrics_file")

    # Move invalid libraries to another directory
    mkdir -p "${output_directory}/parsed_csvs_invalid"
    for library in $invalid_libraries; do
        mv "${output_directory}/parsed_csvs/${library}.csv" "${output_directory}/parsed_csvs_invalid/"
        
        if [ $? -ne 0 ]; then
            echo "Error moving $library.csv" | tee -a "$log_file"
            exit 1
        fi
    done

    # Log the counts before and after filtering
    total_libraries_before=$(wc -l < "$metrics_file")
    valid_libraries_count=$(wc -l < "$valid_metrics_file")
    invalid_libraries_count=$(wc -l < "$invalid_metrics_file")
    echo "Total libraries before filtering: $((total_libraries_before - 1))" >> "${output_directory}/counts.txt"
    echo "Valid libraries after filtering: $((valid_libraries_count - 1))" >> "${output_directory}/counts.txt"
    echo "Invalid libraries after filtering: $((invalid_libraries_count - 1))" >> "${output_directory}/counts.txt"

    # Echo the valid libraries count
    valid_libraries_count=$(wc -l < "$valid_libraries_file")
    echo "Number of valid libraries: $valid_libraries_count" | tee -a "$log_file"

    # Echo the invalid libraries count
    invalid_libraries_count=$(wc -l < "$invalid_metrics_file")
    echo "Number of invalid libraries: $invalid_libraries_count" | tee -a "$log_file"
}

# Function to concatenate CSV files with an optional limit for testing
concatenate_csv_files() {
    mkdir -p "${output_directory}/unfiltered_data"

    # Define output file path
    output_file="${output_directory}/unfiltered_data/unfiltered_iSNVs_all.tsv"

    # Log the start of the concatenation process
    echo "Starting concatenation of CSV files into a single TSV file..." | tee -a "$log_file"

    # Check if a limit is set for the number of files to process (for testing)
    if [ -n "$1" ]; then
        limit=$1
        echo "Testing mode: Only processing the first $limit files." | tee -a "$log_file"
        find "${output_directory}/parsed_csvs" -name "*.csv" | head -n "$limit" > "${output_directory}/file_list.txt"
    else
        find "${output_directory}/parsed_csvs" -name "*.csv" > "${output_directory}/file_list.txt"
    fi

    # Concatenate all CSV files into a single TSV file
    while IFS= read -r csv_file; do
        awk -F ',' '$19 >= 0.01 && $8 >= 100 && $2 >= 101 && $2 <= 29778 {OFS="\t"; $1=$1; print $0}' "$csv_file" >> "$output_file"
    done < "${output_directory}/file_list.txt"

    if [ $? -ne 0 ]; then
        echo "Error concatenating CSV files" | tee -a "$log_file"
        exit 1
    fi

    # Check if the output file is created
    if [ ! -f "$output_file" ]; then
        echo "Error: Output file $output_file was not created." | tee -a "$log_file"
        exit 1
    fi

    # Add header to unfiltered_iSNVs_all.tsv
    header="library_id\tposition_id\tnum_alleles_observed\tA_count\tC_count\tG_count\tT_count\tallele_depth_count\tforward_strand_count\tforward_strand_ratio\treverse_strand_count\treverse_strand_ratio\talternative_allele\talternative_allele_count\talternative_allele_forward_strand_count\talternative_allele_forward_strand_ratio\talternative_allele_reverse_strand_count\talternative_allele_reverse_strand_ratio\talternative_allele_frequency\tstrand_bias_likelihood"
    temp_file="${output_directory}/unfiltered_data/unfiltered_iSNVs_all_with_header.tsv"
    { echo -e "$header"; cat "$output_file"; } > "$temp_file" && mv "$temp_file" "$output_file"

    # # Check if the header was added correctly
    # first_line=$(head -n 1 "$output_file")
    # echo "First line of the output file: $first_line" | tee -a "$log_file"
    # if [ "$first_line" != "$header" ]; then
    #     echo "Error: Header was not added correctly to $output_file." | tee -a "$log_file"
    #     exit 1
    # fi

    # Echo count of unfiltered_iSNVs_all.tsv
    unfiltered_count=$(wc -l < "$output_file")
    echo "Unfiltered all iSNV count is ${unfiltered_count}" | tee -a "$log_file"
}

# Function to zip parsed_csvs folder and remove the unzipped folder
zip_and_remove_unzipped_folder() {
    # Zip the parsed_csvs folder
    zip -r "${output_directory}/parsed_csvs.zip" "${output_directory}/parsed_csvs" | tee -a "$log_file"
    if [ $? -ne 0 ]; then
        echo "Error zipping parsed_csvs folder" | tee -a "$log_file"
        exit 1
    fi

    # Remove the unzipped parsed_csvs folder
    rm -rf "${output_directory}/parsed_csvs" | tee -a "$log_file"
    if [ $? -ne 0 ]; then
        echo "Error removing unzipped parsed_csvs folder" | tee -a "$log_file"
        exit 1
    fi

    # Zip the parsed_csvs_invalid folder
    if [ -d "${output_directory}/parsed_csvs_invalid" ]; then
        zip -r "${output_directory}/parsed_csvs_invalid.zip" "${output_directory}/parsed_csvs_invalid" | tee -a "$log_file"
        if [ $? -ne 0 ]; then
            echo "Error zipping parsed_csvs_invalid folder" | tee -a "$log_file"
            exit 1
        fi

        # Remove the unzipped parsed_csvs_invalid folder
        rm -rf "${output_directory}/parsed_csvs_invalid" | tee -a "$log_file"
        if [ $? -ne 0 ]; then
            echo "Error removing unzipped parsed_csvs_invalid folder" | tee -a "$log_file"
            exit 1
        fi
    fi
}

# Function to create matrices of unfiltered iSNVs, de novo iSNVs, and consensus iSNVs
create_matrix_unfiltered_iSNVs() {
    # Create matrix for all unfiltered iSNVs
    awk 'BEGIN {FS = "\t"} {print $2 "\t" $1 "\t" $19}' "${output_directory}/unfiltered_data/unfiltered_iSNVs_all.tsv" | \
    sed '1d' > "${output_directory}/unfiltered_data/unfiltered_iSNVs_3cols_all.tsv"
    
    if [ $? -ne 0 ]; then
        echo "Error creating matrix unfiltered iSNVs" | tee -a "$log_file"
        exit 1
    fi
    
    awk 'BEGIN {FS = "\t"} $2!=id && NR!=1 {printf id; for(i=101; i<=29778; i++){printf "\t%s", t[i]} printf "\n"} $2!=id {id=$2; for(i=101; i<=29778; i++){t[i]=0}} $2==id {t[$1]=$3} END {printf id; for(i=101; i<=29778; i++){printf "\t%s", t[i]} printf "\n"}' \
    "${output_directory}/unfiltered_data/unfiltered_iSNVs_3cols_all.tsv" > "${output_directory}/unfiltered_data/unfiltered_iSNVs_all.mat"
    
    if [ $? -ne 0 ]; then
        echo "Error creating matrix for unfiltered_iSNVs_all.mat" | tee -a "$log_file"
        exit 1
    fi

    # Create matrix for unfiltered de novo iSNVs
    awk 'BEGIN {FS = "\t"} {print $2 "\t" $1 "\t" $19}' "${output_directory}/unfiltered_data/unfiltered_iSNVs_denovo.tsv" | \
    sed '1d' > "${output_directory}/unfiltered_data/unfiltered_iSNVs_3cols_denovo.tsv"
    
    if [ $? -ne 0 ]; then
        echo "Error creating matrix unfiltered de novo iSNVs" | tee -a "$log_file"
        exit 1
    fi
    
    awk 'BEGIN {FS = "\t"} $2!=id && NR!=1 {printf id; for(i=101; i<=29778; i++){printf "\t%s", t[i]} printf "\n"} $2!=id {id=$2; for(i=101; i<=29778; i++){t[i]=0}} $2==id {t[$1]=$3} END {printf id; for(i=101; i<=29778; i++){printf "\t%s", t[i]} printf "\n"}' \
    "${output_directory}/unfiltered_data/unfiltered_iSNVs_3cols_denovo.tsv" > "${output_directory}/unfiltered_data/unfiltered_iSNVs_denovo.mat"
    
    if [ $? -ne 0 ]; then
        echo "Error creating matrix for unfiltered_iSNVs_denovo.mat" | tee -a "$log_file"
        exit 1
    fi

    # Create matrix for unfiltered consensus iSNVs
    awk 'BEGIN {FS = "\t"} {print $2 "\t" $1 "\t" $19}' "${output_directory}/unfiltered_data/unfiltered_iSNVs_consensus.tsv" | \
    sed '1d' > "${output_directory}/unfiltered_data/unfiltered_iSNVs_3cols_consensus.tsv"
    
    if [ $? -ne 0 ]; then
        echo "Error creating matrix unfiltered consensus iSNVs" | tee -a "$log_file"
        exit 1
    fi
    
    awk 'BEGIN {FS = "\t"} $2!=id && NR!=1 {printf id; for(i=101; i<=29778; i++){printf "\t%s", t[i]} printf "\n"} $2!=id {id=$2; for(i=101; i<=29778; i++){t[i]=0}} $2==id {t[$1]=$3} END {printf id; for(i=101; i<=29778; i++){printf "\t%s", t[i]} printf "\n"}' \
    "${output_directory}/unfiltered_data/unfiltered_iSNVs_3cols_consensus.tsv" > "${output_directory}/unfiltered_data/unfiltered_iSNVs_consensus.mat"
    
    if [ $? -ne 0 ]; then
        echo "Error creating matrix for unfiltered_iSNVs_consensus.mat" | tee -a "$log_file"
        exit 1
    fi
}

# Function to create filtered_iSNVs_consensus.tsv (AAF > 0.75) and filtered_iSNVs_denovo.tsv (AAF <= 0.75)
create_unfiltered_consensus_and_denovo_files() {
    # Define file paths
    input_file="${output_directory}/unfiltered_data/unfiltered_iSNVs_all.tsv"
    consensus_file="${output_directory}/unfiltered_data/unfiltered_iSNVs_consensus.tsv"
    denovo_file="${output_directory}/unfiltered_data/unfiltered_iSNVs_denovo.tsv"

    # Ensure the input file exists
    if [ ! -f "$input_file" ]; then
        echo "Error: Input file $input_file does not exist." | tee -a "$log_file"
        exit 1
    fi

    # Add headers to output files
    header="library_id\tposition_id\tnum_alleles_observed\tA_count\tC_count\tG_count\tT_count\tallele_depth_count\tforward_strand_count\tforward_strand_ratio\treverse_strand_count\treverse_strand_ratio\talternative_allele\talternative_allele_count\talternative_allele_forward_strand_count\talternative_allele_forward_strand_ratio\talternative_allele_reverse_strand_count\talternative_allele_reverse_strand_ratio\talternative_allele_frequency\tstrand_bias_likelihood"
    echo -e "$header" > "$consensus_file"
    echo -e "$header" > "$denovo_file"

    # Split the input file based on AAF > 0.75 and AAF <= 0.75
    awk -F '\t' 'NR>1 && $19 > 0.75' "$input_file" >> "$consensus_file"
    if [ $? -ne 0 ]; then
        echo "Error filtering consensus iSNVs (AAF > 0.75)." | tee -a "$log_file"
        exit 1
    fi

    awk -F '\t' 'NR>1 && $19 <= 0.75' "$input_file" >> "$denovo_file"
    if [ $? -ne 0 ]; then
        echo "Error filtering denovo iSNVs (AAF <= 0.75)." | tee -a "$log_file"
        exit 1
    fi

    # Echo counts after filtering
    consensus_count=$(wc -l < "$consensus_file")
    denovo_count=$(wc -l < "$denovo_file")
    echo "Consensus iSNVs before filtering: $((consensus_count - 1))" >> "${output_directory}/unfiltered_data/thresholds_used.txt"
    echo "De novo iSNVs before filtering: $((denovo_count - 1))" >> "${output_directory}/unfiltered_data/thresholds_used.txt"
    echo "Unfiltered consensus iSNV count is $((consensus_count - 1))" >> "${output_directory}/counts.txt"
    echo "Unfiltered de novo iSNV count is $((denovo_count - 1))" >> "${output_directory}/counts.txt"
}

# Function to create a file containing iSNVs removed by the filter on strand bias likelihood (S)
create_S_less_than_threshold_file() {
    mkdir -p "${output_directory}/filtered_data"
    awk -F '\t' -v S="$S_threshold" '$19 < 0.75 && $20 < S' \
    "${output_directory}/unfiltered_data/unfiltered_iSNVs_denovo.tsv" > "${output_directory}/filtered_data/S_less_than_${S_threshold}_iSNVs.tsv"
    
    if [ $? -ne 0 ]; then
        echo "Error creating S less than threshold file" | tee -a "$log_file"
        exit 1
    fi
}

# Function to create a file containing genomic positions to be masked based on 
# the iSNVs removed by the filter on strand bias likelihood (S)
run_generate_positions_to_mask() {
    python generate_positions_to_mask.py $S_threshold "${output_directory}/unfiltered_data/unfiltered_iSNVs_denovo.tsv" ${output_directory} | tee -a "$log_file"
    if [ $? -ne 0 ]; then
        echo "Error generating positions to mask" | tee -a "$log_file"
        exit 1
    fi
}

# Function to perform filtering on unfiltered_iSNVs.tsv
# Masking Positions: Rows based on a list of positions to mask are removed. See generate_positions_to_mask.tsv for more detail
# Counts Recording: Records the counts of rows before and after filtering, and the number of unique libraries.
# Error Checking: Includes error checking for each step to ensure smooth execution.
filter_unfiltered_iSNVs() {
    mkdir -p "${output_directory}/filtered_data"

    # Count unfiltered iSNVs and libraries before filtering
    unfiltered_count=$(wc -l < "${output_directory}/unfiltered_data/unfiltered_iSNVs_all.tsv")
    unfiltered_library_count=$(cut -f1 "${output_directory}/unfiltered_data/unfiltered_iSNVs_all.tsv" | sort | uniq | wc -l)
    echo "All iSNVs before filtering: $unfiltered_count" >> "${output_directory}/filtered_data/thresholds_used.txt"
    echo "Number of unique libraries before filtering: $unfiltered_library_count" >> "${output_directory}/filtered_data/thresholds_used.txt"
    echo "All iSNVs before filtering: $unfiltered_count" >> "${output_directory}/counts.txt"
    echo "Number of unique libraries before filtering: $unfiltered_library_count" >> "${output_directory}/counts.txt"

    # Apply filters: AAF, S
    awk -F '\t' -v AAF="$AAF_threshold" -v S="$S_threshold" '{
        if ($19 >= AAF && ($19 > 0.75 || $20 >= S)) {
            print $0;
        }
    }' "${output_directory}/unfiltered_data/unfiltered_iSNVs_all.tsv" > "${output_directory}/filtered_data/temp_filtered_iSNVs_all.tsv"
    
    if [ $? -ne 0 ]; then
        echo "Error filtering unfiltered_iSNVs_all.tsv" | tee -a "$log_file"
        exit 1
    fi

    # Remove positions found in genomic_positions_to_mask.txt
    awk 'NR==FNR { positions[$1]; next } !($2 in positions)' \
    "${output_directory}/masking_results_S_threshold_${S_threshold}/genomic_positions_to_mask.txt" \
    "${output_directory}/filtered_data/temp_filtered_iSNVs_all.tsv" > "${output_directory}/filtered_data/filtered_iSNVs_all.tsv"
    
    if [ $? -ne 0 ]; then
        echo "Error removing positions found in genomic_positions_to_mask.txt" | tee -a "$log_file"
        exit 1
    fi

    # Remove temporary file
    rm "${output_directory}/filtered_data/temp_filtered_iSNVs_all.tsv"
    
    echo "Filtering thresholds used:" > "${output_directory}/filtered_data/thresholds_used.txt"
    echo "depth >= 100X" >> "${output_directory}/filtered_data/thresholds_used.txt"
    echo "AAF >= $AAF_threshold" >> "${output_directory}/filtered_data/thresholds_used.txt"
    echo "S >= $S_threshold" >> "${output_directory}/filtered_data/thresholds_used.txt"
    
    # Count filtered iSNVs and libraries after filtering
    filtered_count=$(wc -l < "${output_directory}/filtered_data/filtered_iSNVs_all.tsv")
    filtered_library_count=$(cut -f1 "${output_directory}/filtered_data/filtered_iSNVs_all.tsv" | sort | uniq | wc -l)
    echo "All iSNVs after filtering: $filtered_count" >> "${output_directory}/filtered_data/thresholds_used.txt"
    echo "Number of unique libraries after filtering: $filtered_library_count" >> "${output_directory}/filtered_data/thresholds_used.txt"
    echo "Filtered all iSNV count is ${filtered_count}" >> "${output_directory}/counts.txt"
    echo "Number of unique libraries after filtering: $filtered_library_count" >> "${output_directory}/counts.txt"
}

# Function to create filtered_iSNVs_consensus.tsv (AAF > 0.75) and filtered_iSNVs_denovo.tsv (AAF <= 0.75)
create_filtered_consensus_and_denovo_files() {
    # Define the input and output file paths
    input_file="${output_directory}/filtered_data/filtered_iSNVs_all.tsv"
    consensus_file="${output_directory}/filtered_data/filtered_iSNVs_consensus.tsv"
    denovo_file="${output_directory}/filtered_data/filtered_iSNVs_denovo.tsv"
    
    # Add headers to the output files
    head -n 1 "$input_file" > "$consensus_file"
    head -n 1 "$input_file" > "$denovo_file"
    
    if [ $? -ne 0 ]; then
        echo "Error creating header for filtered consensus and denovo files" | tee -a "$log_file"
        exit 1
    fi

    # Filter the rows and append to the respective files, skipping the header
    tail -n +2 "$input_file" | awk -F '\t' '$19 > 0.75' >> "$consensus_file"
    tail -n +2 "$input_file" | awk -F '\t' '$19 <= 0.75' >> "$denovo_file"
    
    if [ $? -ne 0 ]; then
        echo "Error filtering consensus and denovo files" | tee -a "$log_file"
        exit 1
    fi

    echo "consensus vs de novo is set at 0.75" >> "${output_directory}/filtered_data/thresholds_used.txt"
    
    consensus_count=$(wc -l < "$consensus_file")
    denovo_count=$(wc -l < "$denovo_file")
    echo "Consensus iSNVs after filtering: $((consensus_count - 1))" >> "${output_directory}/filtered_data/thresholds_used.txt"
    echo "De novo iSNVs after filtering: $((denovo_count - 1))" >> "${output_directory}/filtered_data/thresholds_used.txt"
    echo "Filtered consensus iSNV count is $((consensus_count - 1))" >> "${output_directory}/counts.txt"
    echo "Filtered de novo iSNV count is $((denovo_count - 1))" >> "${output_directory}/counts.txt"
}

# Function to identify outlier libraries
identify_outlier_libraries() {
    python identify_outlier_libraries.py "${output_directory}/filtered_data/filtered_iSNVs_denovo.tsv" "${output_directory}" | tee -a "$log_file"
    if [ $? -ne 0 ]; then
        echo "Error identifying outlier libraries" | tee -a "$log_file"
        exit 1
    fi
}

# Function to remove 1% outlier libraries
filter_outlier_libraries() {
    local outlier_file="${output_directory}/outlier_libraries_analysis/outlier_libraries.txt"
    
    # Check if outlier libraries file exists
    if [ ! -f "$outlier_file" ]; then
        echo "Error: Outlier libraries file '$outlier_file' not found." | tee -a "$log_file"
        exit 1
    fi

    # Define an array of files to process
    local files_to_filter=(
        "${output_directory}/filtered_data/filtered_iSNVs_all.tsv"
        "${output_directory}/filtered_data/filtered_iSNVs_denovo.tsv"
        "${output_directory}/filtered_data/filtered_iSNVs_consensus.tsv"
    )

    # Loop through each file and filter out the outlier libraries
    for file in "${files_to_filter[@]}"; do
        awk 'NR==FNR { outliers[$1]; next } !($1 in outliers)' \
        "$outlier_file" "$file" > "${file}.tmp"

        if [ $? -ne 0 ]; then
            echo "Error filtering outlier libraries in $file" | tee -a "$log_file"
            exit 1
        fi

        # Replace original file with the filtered file
        mv "${file}.tmp" "$file"

        if [ $? -ne 0 ]; then
            echo "Error renaming temporary file for $file" | tee -a "$log_file"
            exit 1
        fi

        # Extract the filename part
        file_name=$(basename "$file")
        
        # Count the number of iSNVs and unique libraries after removing outliers
        filtered_count=$(wc -l < "$file")
        library_count=$(cut -f1 "$file" | sort | uniq | wc -l)
        echo "Filtered iSNVs count after removing outliers in $file_name: $filtered_count" >> "${output_directory}/filtered_data/thresholds_used.txt"
        echo "Number of unique libraries after removing outliers in $file_name: $library_count" >> "${output_directory}/filtered_data/thresholds_used.txt"
        echo "Filtered iSNVs count after removing outliers in $file_name: $filtered_count" >> "${output_directory}/counts.txt"
        echo "Number of unique libraries after removing outliers in $file_name: $library_count" >> "${output_directory}/counts.txt"
    done
}

# Function to create a matrix of filtered iSNVs
create_matrix_filtered_iSNVs() {
    awk 'BEGIN {FS = "\t"} {print $2 "\t" $1 "\t" $19}' "${output_directory}/filtered_data/filtered_iSNVs_all.tsv" | sed '1d' > "${output_directory}/filtered_data/filtered_iSNVs_3cols_all.tsv"
    if [ $? -ne 0 ]; then
        echo "Error creating matrix filtered iSNVs (all)" | tee -a "$log_file"
        exit 1
    fi

    awk 'BEGIN {FS = "\t"} {print $2 "\t" $1 "\t" $19}' "${output_directory}/filtered_data/filtered_iSNVs_denovo.tsv" | sed '1d' > "${output_directory}/filtered_data/filtered_iSNVs_3cols_denovo.tsv"
    if [ $? -ne 0 ]; then
        echo "Error creating matrix filtered iSNVs (denovo)" | tee -a "$log_file"
        exit 1
    fi

    awk 'BEGIN {FS = "\t"} {print $2 "\t" $1 "\t" $19}' "${output_directory}/filtered_data/filtered_iSNVs_consensus.tsv" | sed '1d' > "${output_directory}/filtered_data/filtered_iSNVs_3cols_consensus.tsv"
    if [ $? -ne 0 ]; then
        echo "Error creating matrix filtered iSNVs (consensus)" | tee -a "$log_file"
        exit 1
    fi
    
    awk 'BEGIN {FS = "\t"} $2!=id && NR!=1 {printf id; for(i=101; i<=29778; i++){printf "\t%s", t[i]} printf "\n"} $2!=id {id=$2; for(i=101; i<=29778; i++){t[i]=0}} $2==id {t[$1]=$3} END {printf id; for(i=101; i<=29778; i++){printf "\t%s", t[i]} printf "\n"}' "${output_directory}/filtered_data/filtered_iSNVs_3cols_all.tsv" > "${output_directory}/filtered_data/filtered_iSNVs_all.mat"
    if [ $? -ne 0 ]; then
        echo "Error creating matrix for filtered_iSNVs_all.mat" | tee -a "$log_file"
        exit 1
    fi

    awk 'BEGIN {FS = "\t"} $2!=id && NR!=1 {printf id; for(i=101; i<=29778; i++){printf "\t%s", t[i]} printf "\n"} $2!=id {id=$2; for(i=101; i<=29778; i++){t[i]=0}} $2==id {t[$1]=$3} END {printf id; for(i=101; i<=29778; i++){printf "\t%s", t[i]} printf "\n"}' "${output_directory}/filtered_data/filtered_iSNVs_3cols_denovo.tsv" > "${output_directory}/filtered_data/filtered_iSNVs_denovo.mat"
    if [ $? -ne 0 ]; then
        echo "Error creating matrix for filtered_iSNVs_denovo.mat" | tee -a "$log_file"
        exit 1
    fi

    awk 'BEGIN {FS = "\t"} $2!=id && NR!=1 {printf id; for(i=101; i<=29778; i++){printf "\t%s", t[i]} printf "\n"} $2!=id {id=$2; for(i=101; i<=29778; i++){t[i]=0}} $2==id {t[$1]=$3} END {printf id; for(i=101; i<=29778; i++){printf "\t%s", t[i]} printf "\n"}' "${output_directory}/filtered_data/filtered_iSNVs_3cols_consensus.tsv" > "${output_directory}/filtered_data/filtered_iSNVs_consensus.mat"
    if [ $? -ne 0 ]; then
        echo "Error creating matrix for filtered_iSNVs_consensus.mat" | tee -a "$log_file"
        exit 1
    fi
}

# Function to check if unfiltered data file exists
check_unfiltered_data() {
    if [ -f "${output_directory}/unfiltered_data/unfiltered_iSNVs_all.tsv" ]; then
        echo "Unfiltered data file already exists. Skipping steps before filtering." | tee -a "$log_file"
        return 1
    fi
}

# Function to identify mixed libraries
run_identify_mixed_libraries() {
    echo "Identifying mixed libraries..." | tee -a "$log_file"
    
    # Define paths to input files
    local consensus_file=$1
    local denovo_file=$2
    
    # Run the Python script
    python identify_mixed_libraries.py "$consensus_file" "$denovo_file" "$metadata_file" "$output_directory"
    
    if [ $? -ne 0 ]; then
        echo "Error identifying mixed libraries" | tee -a "$log_file"
        exit 1
    fi
}

# Function to remove mixed libraries
remove_mixed_libraries() {
    mkdir -p "${output_directory}/filtered_data_noMixed"

    # List of mixed libraries
    mixed_libraries_file="${output_directory}/mixed_libraries_analysis/possible_mixed_libraries.txt"
    if [ ! -f "$mixed_libraries_file" ]; then
        echo "Error: Possible mixed libraries file '$mixed_libraries_file' not found." | tee -a "$log_file"
        exit 1
    fi

    # Loop through each file in the filtered_data folder
    for file in "${output_directory}/filtered_data"/*; do
        if [[ $file == *.tsv ]]; then
            output_file="${output_directory}/filtered_data_noMixed/$(basename "$file" .tsv)_noMixed.tsv"
            awk 'NR==FNR { mixed_libraries[$1]; next } !($1 in mixed_libraries)' "$mixed_libraries_file" "$file" > "$output_file"
        fi
        if [[ $file == *.mat ]]; then
            output_file="${output_directory}/filtered_data_noMixed/$(basename "$file" .mat)_noMixed.mat"
            awk 'NR==FNR { mixed_libraries[$1]; next } !($1 in mixed_libraries)' "$mixed_libraries_file" "$file" > "$output_file"
        fi
    done

    if [ $? -ne 0 ]; then
        echo "Error removing mixed libraries" | tee -a "$log_file"
        exit 1
    fi

    echo "Mixed libraries removed and saved to filtered_data_noMixed" | tee -a "$log_file"
}

# Function to run annotate_iSNVs.py on all *iSNVs_noMixed.tsv files in filtered_data_noMixed
annotate_iSNVs() {
    local output_directory=$1
    local test_mode=$2  # Optional argument to enable test mode

    # Ensure the output directory exists
    if [ ! -d "$output_directory" ]; then
        mkdir -p "$output_directory"
        if [ $? -ne 0 ]; then
            echo "Error: Cannot create directory $output_directory. Permission denied or invalid path." | tee -a "$log_file"
            exit 1
        fi
    fi

    for file in "${output_directory}"/*noMixed.tsv; do
        if [[ -f "$file" ]]; then
            if [ "$test_mode" == "test" ]; then
                echo "Annotating first 1000 lines of $file" | tee -a "$log_file"
                head -n 1000 "$file" | python annotate_iSNVs.py /dev/stdin "$reference_file" "$metadata_file" "$nt_fitness_file" "$output_directory" | tee -a "$log_file"
            else
                echo "Annotating $file" | tee -a "$log_file"
                python annotate_iSNVs.py "$file" "$reference_file" "$metadata_file" "$nt_fitness_file" "$output_directory" | tee -a "$log_file"
            fi
            if [ $? -ne 0 ]; then
                echo "Error annotating $file" | tee -a "$log_file"
                exit 1
            fi
        fi
    done
}

# Function to run compute_PNN.py on all *_embedding.tsv files in the specified directory
compute_PNN() {
    local embedding_directory=$1
    
    for file in "${embedding_directory}"/*_embeddings.tsv; do
        if [[ -f "$file" ]]; then
            echo "Computing PNN for $file" | tee -a "$log_file"
            python compute_PNN.py "$file" "$metadata_file" "$embedding_directory" | tee -a "$log_file"
            if [ $? -ne 0 ]; then
                echo "Error computing PNN for $file" | tee -a "$log_file"
                exit 1
            fi
        fi
    done
}

convert_tsv_to_mat() {
    local input_file=$1
    local output_file="${input_file%.tsv}.mat"
    
    awk 'BEGIN {FS = "\t"} {print $2 "\t" $1 "\t" $19}' "$input_file" | sed '1d' > "${input_file%.tsv}_3cols.tsv"
    
    awk 'BEGIN {FS = "\t"} $2!=id && NR!=1 {printf id; for(i=101; i<=29778; i++){printf "\t%s", t[i]} printf "\n"} $2!=id {id=$2; for(i=101; i<=29778; i++){t[i]=0}} $2==id {t[$1]=$3} END {printf id; for(i=101; i<=29778; i++){printf "\t%s", t[i]} printf "\n"}' \
    "${input_file%.tsv}_3cols.tsv" > "$output_file"
    
    rm "${input_file%.tsv}_3cols.tsv"
}

# Main function
main() {
    # Assign input arguments to variables
    AAF_threshold="$2"
    S_threshold="$4"
    input_file_paths="$6"
    output_directory="$8"
    
    # In your output directory make sure you have the following files if you want the following functions to run:
    # run_identify_mixed_libraries
    # remove_mixed_libraries
    # annotate_iSNVs
    reference_file="${output_directory}/reference.tsv"
    metadata_file="${output_directory}/metadata.tsv"
    nt_fitness_file="${output_directory}/nt_fitness.csv"

    # Create log file
    log_file="${output_directory}/pipeline_log.txt"
    touch "$log_file"
    
    # Check input arguments
    check_input_arguments "$@"
    
    # Check if the input file paths and output directory exist
    if [ ! -f "$input_file_paths" ]; then
        echo "Error: Input file paths '$input_file_paths' not found." | tee -a "$log_file"
        exit 1
    fi

    if [ ! -d "$output_directory" ]; then
        echo "Error: Output directory '$output_directory' not found." | tee -a "$log_file"
        exit 1
    fi

    # Determine path to the script directory
    script_directory="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

    # Activate Python environment
    source ~/projects/ctb-hussinju/fmost/virtual_environments/msphate/bin/activate
    if [ $? -ne 0 ]; then
        echo "Error activating Python environment" | tee -a "$log_file"
        exit 1
    fi

    check_input_arguments "$@"
    check_unfiltered_data
    if [ $? -ne 1 ]; then
        parse_base_files
        check_for_missing_csv_files
        compute_coverage_metrics
        mask_low_coverage_libraries
        concatenate_csv_files
        zip_and_remove_unzipped_folder &
    fi
    create_unfiltered_consensus_and_denovo_files
    create_S_less_than_threshold_file
    run_generate_positions_to_mask
    filter_unfiltered_iSNVs
    create_filtered_consensus_and_denovo_files
    identify_outlier_libraries
    filter_outlier_libraries
    create_matrix_unfiltered_iSNVs &
    create_matrix_filtered_iSNVs &
    run_identify_mixed_libraries "${output_directory}/filtered_data/filtered_iSNVs_consensus.tsv" "${output_directory}/filtered_data/filtered_iSNVs_denovo.tsv"
    remove_mixed_libraries
    
    # Annotate iSNVs
    annotate_iSNVs "${output_directory}/filtered_data_noMixed" &
    annotate_iSNVs "${output_directory}/filtered_data_noMixed_no211"
    annotate_iSNVs "${output_directory}/filtered_data_noMixed_no211_noTop1"
    annotate_iSNVs "${output_directory}/filtered_data_noMixed_no211_noTop5"
    
    # Compute PNN
    compute_PNN "${output_directory}/filtered_data_noMixed" &
    
    # Compute PNN
    compute_PNN "${output_directory}/filtered_data" &

    # Compute PNN
    compute_PNN "${output_directory}/unfiltered_data" &

    # Compute PNN
    compute_PNN "${output_directory}/filtered_data_noMixed_no211" &
    compute_PNN "${output_directory}/filtered_data_noMixed_no211_noTop1" &
    compute_PNN "${output_directory}/filtered_data_noMixed_no211_noTop5" &
}

# Run the main function
main "$@"