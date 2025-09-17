# SV-pipline

Based on the benchmarking study of Liu et al., four optimal combinations were selected for HiFi read–based SV calling: **pbmm2 + cuteSV2**, **minimap2 + Sniffles2**, **ngmlr + SVIM**, and **Winnowmap + pbsv**.  
In addition, two assembly-based approaches were applied, namely **minimap2 + SVIM-asm** and **nucmer + SyRI**.


## Workflow

The overall workflow of SV-pipline is illustrated in the following diagram:

![SV Pipeline Flowchart](pictures/SV-pipline.png)


## Installation

```bash
# Clone the repository
git clone https://github.com/Santiana123/SV-pipline.git
cd SV-pipline

# Create a conda environment using the provided YAML file
conda env create -f SV-pipline.yaml

# Activate the environment
conda activate SV
```

## Usage on HPC

This pipeline is designed to run on an HPC cluster using PBS job scheduling.

```bash
# Submit the pipeline script
qsub SV-pipline.sh
```

The job submission script (`SV-pipline.sh`) includes typical PBS directives:

```bash
#PBS -j oe                 # merge stdout and stderr
#PBS -q FAFU3              # queue name
#PBS -V                    # export environment variables
#PBS -l nodes=1:ppn=40     # resources: 1 node, 40 CPU cores
```



### Configure your sample

Before submission, edit **SV-pipline.sh** to specify:

```bash
name=sample1
ref=/path/to/reference.fasta
query=/path/to/query.fasta
reads=/path/to/sample1.fastq.gz
```

Other variables:

- `threads=40` → number of CPU threads  
- `query` → path to query genome (if assembly-based pipeline is used)  
- `read_group` → read group information (auto-generated from `name`)  

---

## Example

```bash
# Submit a job
qsub SV-pipline.sh
```

Results will be written to directories such as:

- `1.pbmm2-cutesv2/`  
- `2.minimap2-sniffles2/`  
- `3.ngmlr-svim/`  
- `4.winnowmap-pbsv/`  
- `5.minimap2-SVIM-asm/`
- `6.nucmer-SyRI/`
- `7.merge/`
Each folder contains the intermediate alignment files and SV calls.


## Output files description

Typical outputs include:

- **BAM files** → aligned reads (`*.bam`, `*.bam.bai`)  
- **VCF files** → structural variant calls (`*.vcf`, `*.vcf.gz`)  
- **Log files** → runtime logs from each step (`*.log`)  
- **Intermediate index files** → mapping indexes (`*.mmi`, `*.fai`)  

Each pipeline directory is self-contained and can be inspected independently.

---

## Citation

If you use this pipeline, please cite the benchmarking study:

Liu Z, Xie Z, Li M. Comprehensive and deep evaluation of structural variation detection pipelines with third-generation sequencing data[J]. Genome Biology, 2024, 25(1): 188.
