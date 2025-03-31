#!/bin/bash

# Number of iterations (you can change this value)
iterations=10

for ((i=1; i<=iterations; i++))
do
    # Generate a random angle between 0 and pi/2
    angle=$(echo "scale=6; ($RANDOM/32767) * (a(1)*4)" | bc -l)
    
    # Generate a random error between 0.001 and 0.1
    error=$(echo "scale=6; 0.001 + ($RANDOM/32767) * (0.1 - 0.001)" | bc -l)
    
    # Run the program with the generated angle and error
    echo "Running iteration $i with angle=$angle and error=$error"
    ./process $angle $error
done
