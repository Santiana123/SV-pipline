#!/bin/bash
# Filter PBSV VCF by DP field

if [ $# -ne 4 ]; then
    echo "Usage: $0 <input.vcf> <output.vcf> <depth_min> <depth_max>"
    exit 1
fi

in_vcf=$1
out_vcf=$2
depth_min=$3
depth_max=$4

awk -F'\t' -v min=$depth_min -v max=$depth_max '$0 ~ "^#" {print; next}
{
    split($9,f,":"); for(i=1;i<=length(f);i++) if(f[i]=="DP") d=i;
    split($10,v,":");
    if(v[d] != "." && v[d] + 0 >= min && v[d] + 0 <= max) print
}' "$in_vcf" > "$out_vcf"
