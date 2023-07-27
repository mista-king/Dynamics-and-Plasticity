#!/bin/bash

#SBATCH --job-name=cellranger
#SBATCH --partition=64c512g
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=wcj342403669@sjtu.edu.cn
#SBATCH -N 1
#SBATCH --exclusive

export PATH=/dssg/home/acct-medcl1/medcl1/opt/cellranger-7.0.0:$PATH

cellranger vdj --id=cellranger \
--fastqs=/dssg/home/acct-medcl1/medcl1/HCC_TCR_GSE206325/R1006/1T \
--reference=/dssg/home/acct-medcl1/medcl1/opt/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0 \
--sample=SRR19718628,SRR19718629,SRR19718630,SRR19718631