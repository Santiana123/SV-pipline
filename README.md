# SV-pipline

Based on the benchmarking study of Liu et al., four optimal combinations were selected for HiFi read–based SV calling: **pbmm2 + cuteSV2**, **minimap2 + Sniffles2**, **ngmlr + SVIM**, and **Winnowmap + pbsv**.  
In addition, two assembly-based approaches were applied, namely **minimap2 + SVIM-asm** and **nucmer + SyRI**.

---

## Workflow

The overall workflow of SV-pipline is illustrated in the following diagram:

![SV Pipeline Flowchart](pictures/SV-pipline.png)

---

## Installation

```bash
# Clone the repository
git clone https://github.com/你的用户名/你的仓库名.git
cd 你的仓库名

# Create a conda environment using the provided YAML file
conda env create -f SV-pipline.yaml

# Activate the environment
conda activate sv-pipline
