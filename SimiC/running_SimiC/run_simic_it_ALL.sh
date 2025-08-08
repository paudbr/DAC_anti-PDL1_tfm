#!/bin/bash

# --------------------------------------------------------------------
# Title: run_simic_it_ctrl.sh
# Author: Paula de Blas Rioja
# Description: Run SimiC analysis for a given directory across
#              specified subdirectories and treatments. (FOR BOTH CELL LINES MIXED)
# Usage: ./run_simic_it_ctrl.sh ALL <lambda1> <lambda2>
# --------------------------------------------------------------------
# Instructions:
# 1. Load required Python environment and dependencies.
# 2. Set run paths and parameters (directory, lambda1, lambda2, subdirectories, treatments).
# 3. Run the SimiC script for each combination of subdirectory and treatment.
# 4. Save the output weight matrix and logs.
# --------------------------------------------------------------------

# Capture arguments
directory=$1  # First argument: main directory (ALL)
lambda1=$2    # Second argument: lambda1 parameter
lambda2=$3    # Third argument: lambda2 parameter

# Create the logs directory if it does not already exist
mkdir -p ~/workdir/HDAC_tfm/logs

# Define list of subdirectories (cell types or conditions)
subdirectories=("Macrophages" "Immune")

# Define list of treatments to run
treatments=("ALL")

# Loop over all combinations of subdirectories and treatments
for subdirectory in "${subdirectories[@]}"; do
    for treatment in "${treatments[@]}"; do
        
        # Define log file path with parameters in the filename
        log_file="$HOME/workdir/HDAC_tfm/logs/simic_run_${directory}_${subdirectory}_ctrl_${treatment}_L1_${lambda1}_L2_${lambda2}_2.log"
        
        # Print status message
        echo "Running: prefix=$directory, subdirectory=$subdirectory, treatment=$treatment"

        # Run the Python script with nohup so it continues even if the terminal closes
        # Arguments: directory, subdirectory, treatment, lambda1, lambda2
        # Redirect both stdout and stderr to the log file
        nohup python3 ~/workdir/HDAC_tfm/scripts/SimiC_bien/run_SimiC_it_withR_ALL.py \
            "$directory" "$subdirectory" "$treatment" "$lambda1" "$lambda2" > "$log_file" 2>&1 &
        
        # Store the process ID of the last background command
        pid=$!

        # Wait for the process to finish before starting the next combination
        wait $pid

        # Print completion message
        echo "Finished: prefix=$directory, subdirectory=$subdirectory, treatment=$treatment"
    done
done