#!/bin/bash

file=$1

while read -r vcf_line; do
  # Skip header lines in VCF file
  [[ $vcf_line =~ ^# ]] && continue

  # Extract quality and INFO fields
  quality=$(echo "$vcf_line" | cut -f 6)
  info_field=$(echo "$vcf_line" | cut -f 8)

  # Extract DP value from the INFO field
  DP=$(echo "$info_field" | awk -F';' '{for(i=1;i<=NF;i++) if($i ~ /^DP=/) print substr($i,4)}')

  # Use bc for floating-point comparisons
  quality_check=$(echo "$quality > 30" | bc)
  DP_check=$(echo "$DP > 5" | bc)

  if [[ $quality_check -eq 1 && $DP_check -eq 1 ]]; then
    chromosome=$(echo "$vcf_line" | cut -f 1)
    position=$(echo "$vcf_line" | cut -f 2)
    reference=$(echo "$vcf_line" | cut -f 4)
    alternative=$(echo "$vcf_line" | cut -f 5)
    echo -e "$chromosome\t$position\t$reference\t$alternative"
  fi
done < "$file" > "$file.hc_calls.list"
