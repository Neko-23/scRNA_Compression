#!/bin/bash

# Define arrays for the outer and inner loops
outer_elements=(4 5 6 7 8)
inner_numbers=(2 5 10 20 30)

# Outer loop
for sample in "${outer_elements[@]}"; do
  # Inner loop
  for k in "${inner_numbers[@]}"; do
    echo "Sample $sample with $k clusters"
    python compression_pipeline.py $sample $k
  done
done