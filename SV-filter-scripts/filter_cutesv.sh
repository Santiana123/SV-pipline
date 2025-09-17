#!/bin/bash
# Filter CuteSV VCF by depth range

if [ $# -ne 4 ]; then
    echo "Usage: $0 <input.vcf> <output.vcf> <depth_min> <depth_max>"
    exit 1
fi

in_vcf=$1
out_vcf=$2
depth_min=$3
depth_max=$4

awk -F'\t' -v min=$depth_min -v max=$depth_max '!/^#/ && $6 != "." {
    split($9,f,":"); split($10,v,":");
    dr=dv=0;
    for(i=1;i<=length(f);i++) {
        if(f[i]=="DR") dr=v[i];
        if(f[i]=="DV") dv=v[i];
    }
    if(dr+dv >= min && dr+dv <= max) print
}' "$in_vcf" > "$out_vcf"
