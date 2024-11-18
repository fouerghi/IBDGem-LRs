#!/bin/bash

# Define the directories containing the subdirectories
main_directories=($(find "/$HOME/complete_LD" -mindepth 1 -maxdepth 1 -type d))


# Loop through each main directory
for main_directory in "${main_directories[@]}"; do
    # Get the value of n from the directory name
    n_value=$(basename "$main_directory" | cut -d 'n' -f 2 | cut -d '_' -f 1)
    # Get the value of k from the directory name
    k_value=$(basename "$main_directory" | cut -d 'k' -f 2 | cut -d '_' -f 1)
    
    for subdir in "$main_directory"/*/*; do
        if [[ "$subdir" == *"downsampled"* ]]; then
            echo "Processing subdirectory: $subdir"
            if [[ -e "$subdir/test.hap" && -e "$subdir/test.legend" && -e "$subdir/test.indv" && -e "$subdir/unknowndna.pileup" ]]; then
                sample_value="sample$n_value"
                ./IBDGem-2.0.2/ibdgem --LD -H "$subdir/test.hap" -L "$subdir/test.legend" -I "$subdir/test.indv" -P "$subdir/unknowndna.pileup" -N "$sample_value" -O "$subdir/output_LD" --allele-freqs "$subdir/output.frq" --window-size "$k_value" --sample "$sample_value"
                ./IBDGem-2.0.2/ibdgem -H "$subdir/test.hap" -L "$subdir/test.legend" -I "$subdir/test.indv" -P "$subdir/unknowndna.pileup" -N "$sample_value" -O "$subdir/output_non_LD" --allele-freqs "$subdir/output.frq" --window-size "$k_value" --sample "$sample_value"
            else
                echo "Files not found in $subdir"
            fi
            echo "Done processing $subdir"
        else
            echo "not downsampled"
        fi
    done
done


for main_directory in "${main_directories[@]}"; do
    # Get the value of n from the directory name
    n_value=$(basename "$main_directory" | cut -d 'n' -f 2 | cut -d '_' -f 1)
    # Get the value of k from the directory name
    k_value=$(basename "$main_directory" | cut -d 'k' -f 2 | cut -d '_' -f 1)
    # the only change from before is now we have /*/* below instead of /* so it goes in the subdirectories
    for subdir in "$main_directory"/*/*; do
        if [[ "$subdir" == *"entire"* ]]; then
            echo "Processing subdirectory: $subdir"
            if [[ -e "$subdir/test.hap" && -e "$subdir/test.legend" && -e "$subdir/test.indv" && -e "$subdir/unknowndna.pileup" ]]; then
                size=$(basename "$subdir" | cut -d 's' -f 2 | cut -d '_' -f 1)
                size=$((size + 1))
                sample_value="sample$size"
                ./IBDGem-2.0.2/ibdgem --LD -H "$subdir/test.hap" -L "$subdir/test.legend" -I "$subdir/test.indv" -P "$subdir/unknowndna.pileup" -N "$sample_value" -O "$subdir/output_LD" --allele-freqs "$subdir/output.frq" --window-size "$k_value" --sample "$sample_value"
            else
                echo "Files not found in $subdir"
            fi
            echo "Done processing $subdir"
        else
            echo "Already processed"
        fi
    done
done
