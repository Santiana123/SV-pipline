#!/bin/bash
# Count and summarize SV types and lengths from filtered VCF

if [ $# -ne 2 ]; then
    echo "Usage: $0 <merge.filtered.vcf> <statistic.txt>"
    exit 1
fi

in_vcf=$1
out_txt=$2

awk -F'\t' '
BEGIN {
    OFS = "\t";
    INS_n = DEL_n = INV_n = DUP_n = TRA_n = 0;
    INS_len = DEL_len = INV_len = DUP_len = 0;
}
!/^#/ {
    svtype = ""; svlen = ""; end = ""; chr2 = "";
    split($8, info, ";");
    for (i in info) {
        if (info[i] ~ /^SVTYPE=/) { split(info[i], a, "="); svtype = a[2]; }
        if (info[i] ~ /^SVLEN=/)  {
            split(info[i], b, "=");
            svlen = b[2];
            gsub(/[^0-9.-]/, "", svlen);
            svlen = (svlen < 0) ? -svlen : svlen;
        }
        if (info[i] ~ /^END=/)    { split(info[i], c, "="); end = c[2]; }
        if (info[i] ~ /^CHR2=/)   { split(info[i], d, "="); chr2 = d[2]; }
    }

    if (svtype == "INS") { INS_n++; if (svlen != "") INS_len += svlen; }
    else if (svtype == "DEL") { DEL_n++; if (svlen != "") DEL_len += svlen; }
    else if (svtype == "DUP") { DUP_n++; if (svlen != "") DUP_len += svlen; }
    else if (svtype == "INV") {
        count_this = 1;
        if ($3 ~ /BND/ && $1 == chr2 && end != "" && ($2 + 0) > (end + 0)) count_this = 0;
        if (count_this) { INV_n++; if (svlen != "") INV_len += svlen; }
    }
    else if (svtype == "TRA") { TRA_n++; }
}
END {
    print "Type\tCount\tTotalLength\tAvgLength"
    print "INS", INS_n, INS_len, (INS_n > 0 ? INS_len / INS_n : 0)
    print "DEL", DEL_n, DEL_len, (DEL_n > 0 ? DEL_len / DEL_n : 0)
    print "INV", INV_n, INV_len, (INV_n > 0 ? INV_len / INV_n : 0)
    print "DUP", DUP_n, DUP_len, (DUP_n > 0 ? DUP_len / DUP_n : 0)
    print "TRA", TRA_n, "-", "-"
}' "$in_vcf" > "$out_txt"
