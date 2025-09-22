#!/bin/bash
# Extract DUP info from merged VCF

if [ $# -ne 2 ]; then
    echo "Usage: $0 <merge.vcf> <output.txt>"
    exit 1
fi

in_vcf=$1
out_txt=$2

awk -F'\t' '
BEGIN {
    OFS = "\t";
    print "CHR", "POS", "ID", "END", "SVLEN", "TYPE_HINT", "SOURCE";
}
!/^#/ {
    svtype = ""; svlen = "NA"; end = "NA"; type_hint = "DUP"; source = "NA";
    split($8, info, ";");

    for (i in info) {
        if (info[i] ~ /^SVTYPE=DUP/) { svtype = "DUP"; }
        if (info[i] ~ /^END=/)       { split(info[i], a, "="); end = a[2]; }
        if (info[i] ~ /^SVLEN=/)     {
            split(info[i], b, "=");
            svlen = b[2];
            gsub(/[^0-9.-]/, "", svlen);
            svlen = (svlen < 0) ? -svlen : svlen;
        }
    }

    if (svtype == "DUP") {
        if ($3 ~ /^pbsv/) source = "pbsv";
        else if ($3 ~ /^Sniffles2/) source = "Sniffles2";
        else if ($3 ~ /^svim/) source = "svim";

        if ($3 ~ /DUP_TANDEM/) type_hint = "DUP_TANDEM";
        else if ($3 ~ /INS\.DUP/) type_hint = "INS+DUP";

        print $1, $2, $3, end, svlen, type_hint, source;
    }
}' "$in_vcf" | column -t > "$out_txt"
