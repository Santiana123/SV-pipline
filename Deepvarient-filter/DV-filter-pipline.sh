#!/bin/bash
#PBS -j oe 
#PBS -q FAFU1
#PBS -V 
#PBS -l nodes=1:ppn=4

cd $PBS_O_WORKDIR
threads=4

# 输入参数
name="AU9"
in_vcf="${name}.ccs.concat.vcf.gz"
depth_min=9
depth_max=66

# 输出文件
filtered_vcf="${name}.filtered.vcf"
snp_vcf="${name}.snp.vcf"
indel_vcf="${name}.indel.vcf"
snp_count="${name}.snp.count.txt"
indel_stats="${name}.indel.statistic.txt"

# 激活环境
source /public/home/yuejingjing/.bashrc
source activate /public/home/yuejingjing/mambaforge/envs/test

# Step 1: DP过滤
bash filter_vcf_by_dp.sh "$in_vcf" "$filtered_vcf" $depth_min $depth_max

# Step 2: 拆分SNP/Indel
bash split_snp_indel.sh "$filtered_vcf" "$snp_vcf" "$indel_vcf"

# Step 3: 统计SNP数
bash count_snp.sh "$snp_vcf" "$snp_count"

# Step 4: 统计Indel指标
bash indel_stats.sh "$indel_vcf" "$indel_stats"
