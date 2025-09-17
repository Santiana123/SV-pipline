#!/bin/bash
# Filter SVIM VCF by SUPPORT field and PASS status

if [ $# -ne 4 ]; then
    echo "Usage: $0 <input.vcf> <output.vcf> <depth_min> <depth_max>"
    exit 1
fi

in_vcf=$1
out_vcf=$2
depth_min=$3
depth_max=$4

awk -F'\t' -v min=$depth_min -v max=$depth_max '!/^#/ && $7 == "PASS" {
    match($8, /SUPPORT=([0-9]+)/, a);
    if (a[1] >= min && a[1] <= max) print
}' "$in_vcf" > "$out_vcf"
