#!/bin/bash

# Fixed contraction parameter
CONTRACTION=0.3

# Tolerance values
TOL_VALUES=(0.1 0.01 0.001 0.0001 0.00001 0.000001)

# Generate 10 random angles between -pi/2 and pi/2
for i in {1..10}; do
    # Random number between 0 and 1
    RANDOM_NUM=$(awk 'BEGIN{srand(); print rand()}' 2>/dev/null)
    
    # Scale to range between -pi/2 and pi/2
    ANGLE=$(awk -v r="$RANDOM_NUM" 'BEGIN{print (r * 3.14159265359) - (3.14159265359/2)}' 2>/dev/null)
    
    # Run the program for each tolerance value
    for TOL in "${TOL_VALUES[@]}"; do
        ./main "$ANGLE" "$TOL" "$CONTRACTION" 
    done
done
