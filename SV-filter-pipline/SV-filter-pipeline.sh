#!/bin/bash
# SV merging and post-analysis pipeline
# Author: Yuejingjing
# Description: Filter individual caller VCFs, merge results, adjust BND, and generate summary statistics.

# Usage:
#   bash SV-merge-pipeline.sh <depth_min> <depth_max>
#
# Example:
#   bash SV-merge-pipeline.sh 9 71

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <depth_min> <depth_max>"
    exit 1
fi

depth_min=$1
depth_max=$2

echo "=== SV Merging Pipeline Started ==="
echo "Depth range: $depth_min - $depth_max"

# Step 1: Filtering
echo "[1/6] Filtering VCFs ..."
bash scripts/filter_cutesv.sh cutesv.vcf cutesv.f.vcf $depth_min $depth_max
bash scripts/filter_sniffles2.sh sniffles2.vcf sniffles2.f.vcf $depth_min $depth_max
bash scripts/filter_svim.sh svim.vcf svim.f.vcf $depth_min $depth_max
bash scripts/filter_pbsv.sh pbsv.vcf pbsv.f.vcf $depth_min $depth_max

# Step 2: Merge with SURVIVOR
echo "[2/6] Merging VCFs with SURVIVOR ..."
ls SyRI.vcf svim-asm.vcf pbsv.f.vcf svim.f.vcf cutesv.f.vcf sniffles2.f.vcf > sample
SURVIVOR merge sample 1000 4 1 1 0 50 merge.vcf

# Step 3: Extract INV breakpoints
echo "[3/6] Extracting INV breakpoints ..."
bash scripts/INV-breakpoints-finder.sh merge.vcf INV.breakpoint.txt

# Step 4: Adjust BND with Python script
echo "[4/6] Adjusting BND (filter misclassified INV) ..."
python3 scripts/BND-adjust.py
grep -v -Ff deleted_ids.txt merge.vcf > merge.filtered.vcf

# Step 5: Statistics
echo "[5/6] Generating statistics ..."
bash scripts/stat_sv.sh merge.filtered.vcf statistic.number+length.txt

# Step 6: Extract DUP info
echo "[6/6] Extracting DUP info ..."
bash scripts/extract_dup.sh merge.vcf DUP.location.txt

echo "=== SV Merging Pipeline Finished Successfully ==="
echo "Results:"
echo " - merge.filtered.vcf"
echo " - statistic.number+length.txt"
echo " - DUP.location.txt"
echo " - INV.breakpoint.txt / INV.filtered.txt"
