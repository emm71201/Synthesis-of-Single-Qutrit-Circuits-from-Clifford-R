#!/bin/bash

# Function to generate a random double between two values
generate_random_double() {
    local min=$1
    local max=$2
    awk -v min="$min" -v max="$max" 'BEGIN { srand(); print min + (max - min) * rand() }'
}

# Function to generate a random integer between two values
generate_random_int() {
    local min=$1
    local max=$2
    echo $((RANDOM % (max - min + 1) + min))
}

# Loop to run the program 100 times
for i in {1..100}; do
    # Generate random values
    angle=$(generate_random_double -1.5708 1.5708) # -pi/2 to pi/2
    tol=$(generate_random_double 0.01 0.5)         # 0.01 to 0.5
    k=$(generate_random_int 5 17)                  # 5 to 17

    # Run the program with the generated values
    ./process "$angle" "$tol" "$k"
done