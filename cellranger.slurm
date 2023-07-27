#!/bin/bash

#SBATCH --job-name=cellranger
#SBATCH --partition=64c512g
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=wcj342403669@sjtu.edu.cn
#SBATCH -N 1
#SBATCH --exclusive

export PATH=/dssg/home/acct-medcl1/medcl1/opt/cellranger-3.1.0:$PATH

cellranger count --id=cellranger \
--fastqs=/dssg/home/acct-medcl1/medcl1/Ovarian_carcinoma_GSE184880/HGSOC4 \
--transcriptome=/dssg/home/acct-medcl1/medcl1/opt/refdata-gex-GRCh38-2020-A \
--sample=SRR16093365,SRR16093366,SRR16093367,SRR16093368