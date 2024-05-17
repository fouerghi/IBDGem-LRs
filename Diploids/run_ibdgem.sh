#/bin/bash

SCRIPT_DIR="/Users/ouerghi1"

# Change to the script directory
cd "$SCRIPT_DIR" || {
    echo "Failed to change directory to $SCRIPT_DIR. Exiting."
    exit 1
}

echo "Starting to run all scripts in parallel"

./run_ibdgem_panel_size.sh &
PID1=$!

./run_ibdgem_read_number.sh &
PID2=$!

./run_ibdgem_loci_number.sh &
PID3=$!

./run_ibdgem_error_rate.sh &
PID1=$!

wait $PID1
if [ $? -ne 0 ]; then
    echo "run_ibdgem_panel_size.sh failed"
    exit 1
fi

wait $PID2
if [ $? -ne 0 ]; then
    echo "run_ibdgem_panel_read_number.sh failed"
    exit 1
fi


wait $PID3
if [ $? -ne 0 ]; then
    echo "run_ibdgem_panel_loci_number.sh failed"
    exit 1
fi

wait $PID4
if [ $? -ne 0 ]; then
    echo "run_ibdgem_error_rate.sh failed"
    exit 1
fi

echo "All scripts ran successfully in parallel"


