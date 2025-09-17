#!/bin/bash
# Extract and filter INV breakpoints from merged VCF

if [ $# -ne 2 ]; then
    echo "Usage: $0 <merge.vcf> "
    exit 1
fi

merge_vcf=$1
out_txt=$2

awk -F'\t' '
BEGIN {
    OFS = "\t";
    print "CHR", "POS", "ID", "CHR2", "END", "SVLEN", "SOURCE", "Counted";
}
!/^#/ {
    svtype = ""; chr2 = "NA"; end = "NA"; svlen = "NA"; counted = "YES"; source = "unknown";
    split($8, info, ";");

    for (i in info) {
        if (info[i] ~ /^SVTYPE=/) { split(info[i], a, "="); svtype = a[2]; }
        if (info[i] ~ /^CHR2=/)   { split(info[i], b, "="); chr2 = b[2]; }
        if (info[i] ~ /^END=/)    { split(info[i], c, "="); end = c[2]; }
        if (info[i] ~ /^SVLEN=/)  {
            split(info[i], d, "=");
            svlen = d[2];
            gsub(/[^0-9.-]/, "", svlen);
            svlen = (svlen < 0) ? -svlen : svlen;
        }
    }

    if (svtype == "INV") {
        if ($3 ~ /^svim_asm/) source = "svim";
        else if ($3 ~ /^pbsv/) source = "pbsv";
        else if ($3 ~ /^Sniffles2/) source = "Sniffles2";

        if ($3 ~ /BND/ && $1 == chr2 && end != "NA" && ($2 + 0) > (end + 0)) {
            counted = "NO";
        }

        print $1, $2, $3, chr2, end, svlen, source, counted;
    }
}' "$merge_vcf" > INV.breakpoint.txt
