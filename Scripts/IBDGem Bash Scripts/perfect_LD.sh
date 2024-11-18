#!/bin/bash

# Define the directories containing the subdirectories
main_directories=($(find "/$HOME/copy_number" -mindepth 1 -maxdepth 1 -type d))

# Loop through each main directory
for main_directory in "${main_directories[@]}"; do
    # Get the value of n from the directory name
    n_value=$(basename "$main_directory" | cut -d 'n' -f 2 | cut -d '_' -f 1)
    # Get the value of k from the directory name
    k_value=$(basename "$main_directory" | cut -d 'k' -f 2 | cut -d '_' -f 1)

    # Loop through each subdirectory in the main directory
    for subdir in "$main_directory"/*/*/*; do
        # Check if the subdirectory contains the necessary files
        if [[ -e "$subdir/test.hap" && -e "$subdir/test.legend" && -e "$subdir/test.indv" && -e "$subdir/unknowndna.pileup" ]]; then
            # Determine the value of -N based on the value of n
            sample_value="sample$n_value"
	    m_value=$(basename "$subdir" | cut -d '_' -f 5 | cut -d 'm' -f 2)
	    window_size=$((k_value * m_value))
	    echo "window size is $window_size"
            # Check if the directory is 'no_redundant' or 'redundant'
            if [[ "$subdir" == *"no_redundant"* ]]; then
                # Run the command for no_redundant directories
                ./IBDGem-2.0.2/ibdgem -H "$subdir/test.hap" -L "$subdir/test.legend" -I "$subdir/test.indv" -P "$subdir/unknowndna.pileup" -N "$sample_value" -O "$subdir/output_non_LD" --allele-freqs "$subdir/output.frq" --sample "$sample_value"
            elif [[ "$subdir" == *"redundant"* ]]; then
                # Run the command for redundant directories
                ./IBDGem-2.0.2/ibdgem --LD -H "$subdir/test.hap" -L "$subdir/test.legend" -I "$subdir/test.indv" -P "$subdir/unknowndna.pileup" -N "$sample_value" -O "$subdir/output_LD" --allele-freqs "$subdir/output.frq" --window-size "$window_size" --sample "$sample_value"
		        ./IBDGem-2.0.2/ibdgem -H "$subdir/test.hap" -L "$subdir/test.legend" -I "$subdir/test.indv" -P "$subdir/unknowndna.pileup" -N "$sample_value" -O "$subdir/output_non_LD" --allele-freqs "$subdir/output.frq" --window-size "$window_size" --sample "$sample_value"
            else
                echo "No valid marker type found in $subdir"
            fi
        else
            echo "Files not found in $subdir"
        fi
	echo "done with $subdir"
    done
done

