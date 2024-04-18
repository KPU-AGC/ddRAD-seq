#!/bin/bash

# Reading the list from the file into an array
mapfile -t items < $1

# Get the number of items in the array
num_items=${#items[@]}

# Generate all unique combinations of two items
for (( i=0; i<$num_items; i++ )); do
    for (( j=i+1; j<$num_items; j++ )); do
        echo -e "${items[i]}\t${items[j]}"
    done
done