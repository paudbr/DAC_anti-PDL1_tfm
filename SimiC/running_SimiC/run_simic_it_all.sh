#!/bin/bash
# --------------------------------------------------------------------
# Title: run_simic_it_all.sh
# Author: Paula de Blas Rioja
# Description: Run SimiC analysis for all specified directory/treatment combinations.
# Usage: ./run_simic_it_all.sh lambda1 lambda2
# --------------------------------------------------------------------
# Instructions:
# 1. Load required Python environment and dependencies.
# 2. Set run paths and parameters (lambda1, lambda2, directories, subdirectories).
# 3. Run the SimiC script for each combination of directory and subdirectory.
# 4. Save the output weight matrix and logs.
# --------------------------------------------------------------------



# Check if exactly 2 arguments are provided (lambda1 and lambda2)
if [ "$#" -ne 2 ]; then
    echo "Usage: lambda1 lambda2"
    exit 1
fi

# Capture the first and second arguments
lambda1=$1
lambda2=$2

# Create the logs directory if it does not already exist
mkdir -p ~/workdir/HDAC_tfm/logs

# Define the list of main directories (cell lines or data prefixes)
directories=("K25L" "KUV")

# Define the list of subdirectories (conditions, sample types, etc.)
subdirectories=("Tumor")

# Loop over all combinations of directories and subdirectories
for directory in "${directories[@]}"; do
    for subdirectory in "${subdirectories[@]}"; do
        
        # Define the log file path with parameters in the filename
        log_file="$HOME/workdir/HDAC_tfm/logs/simic_run_${directory}_${subdirectory}_all_L1_${lambda1}_L2_${lambda2}_2.log"
        
        # Print status message
        echo "Running: prefix=$directory, subdirectory=$subdirectory"

        # Run the Python script with nohup so it continues even if terminal closes
        # Arguments passed: directory, subdirectory, lambda1, lambda2
        # Output and errors are redirected to the log file
        nohup python3 ~/workdir/HDAC_tfm/scripts/SimiC_bien/run_SimiC_it_withR_all.py \
            "$directory" "$subdirectory" "$lambda1" "$lambda2" > "$log_file" 2>&1 &
        
        # Store the process ID of the last background command
        pid=$!

        # Wait for the process to finish before starting the next combination
        wait $pid

        # Print completion message
        echo "Finished: prefix=$directory, subdirectory=$subdirectory (all treatments)"
    done
done

