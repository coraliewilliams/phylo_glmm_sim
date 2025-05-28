#!/bin/bash

# Ask for path to results directory
read -p "Enter the path to the results folder: " result_dir

# Check if folder exists
if [ ! -d "$result_dir" ]; then
  echo "Error: Directory not found."
  exit 1
fi

# Ask for expected job ID range
read -p "Enter the starting job ID: " start_id
read -p "Enter the ending job ID: " end_id

# Extract existing job IDs from files like res_88598, res_90000, etc.
existing_ids=()
for file in "$result_dir"/res_*; do
  basename=$(basename "$file")
  jobid=${basename#res_}
  if [[ $jobid =~ ^[0-9]+$ ]]; then
    existing_ids+=("$jobid")
  fi
done

# Sort and deduplicate
existing_ids_sorted=($(printf "%s\n" "${existing_ids[@]}" | sort -n | uniq))

# Identify missing IDs
echo "Missing job IDs:"
for (( id=$start_id; id<=$end_id; id++ )); do
  if ! printf "%s\n" "${existing_ids_sorted[@]}" | grep -q -w "$id"; then
    echo "$id"
  fi
done
