#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <angle> <start_tol> <end_tol> <num_tols>"
    exit 1
fi

# Set the angle, tolerances, and number of tolerances from the arguments
angle="$1"  # Fixed angle
start="$2"  # Starting tolerance
end="$3"    # Ending tolerance
num_tols="$4"  # Number of tolerance values
output_file="output.log"  # File to capture output

# Calculate tolerances and run the process
for i in $(seq 0 $((num_tols - 1))); do
    tol=$(echo "scale=4; $start - $i * ($start - $end) / ($num_tols - 1)" | bc)
    # Run the process with the specified angle and tolerance, and capture output
    ./process "$angle" "$tol" >> "$output_file"
done
