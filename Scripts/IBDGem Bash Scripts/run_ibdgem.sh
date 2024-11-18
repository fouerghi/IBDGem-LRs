#!/bin/bash

SCRIPT_DIR="$HOME"

# Change to the script directory
cd "$SCRIPT_DIR" || {
    echo "Failed to change directory to $SCRIPT_DIR. Exiting."
    exit 1
}

echo "Starting to run all scripts in parallel"

./panel_size.sh &
PID1=$!

./read_number.sh &
PID2=$!

./loci_number.sh &
PID3=$!

./error_rate.sh &
PID4=$!

./perfect_LD.sh &
PID5=$!

./complete_LD.sh &
PID6=$!

wait $PID1
if [ $? -ne 0 ]; then
    echo "panel_size.sh failed"
    exit 1
fi

wait $PID2
if [ $? -ne 0 ]; then
    echo "read_number.sh failed"
    exit 1
fi

wait $PID3
if [ $? -ne 0 ]; then
    echo "loci_number.sh failed"
    exit 1
fi

wait $PID4
if [ $? -ne 0 ]; then
    echo "error_rate.sh failed"
    exit 1
fi

wait $PID5
if [ $? -ne 0 ]; then
    echo "perfect_LD.sh failed"
    exit 1
fi

wait $PID6
if [ $? -ne 0 ]; then
    echo "complete_LD.sh failed"
    exit 1
fi

echo "All scripts ran successfully in parallel"
