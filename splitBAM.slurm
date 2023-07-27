#!/bin/bash

#SBATCH --job-name=splitBAM
#SBATCH --partition=64c512g
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=wcj342403669@sjtu.edu.cn
#SBATCH -n 5
#SBATCH --array=1-1990

module load miniconda3
source activate mypython3
module load samtools
id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /dssg/home/acct-medcl1/medcl1/Ovarian_carcinoma_GSE184880/HGSOC4/bam/xac)
export BAM_FILE='/dssg/home/acct-medcl1/medcl1/Ovarian_carcinoma_GSE184880/HGSOC4/cellranger/outs/possorted_genome_bam.bam'
#samtools view -@ 16 -H $BAM_FILE > SAM_header
samtools view -@ 16 $BAM_FILE | LC_ALL=C grep -F ${id} > ${id}_SAM_body
cat /dssg/home/acct-medcl1/medcl1/Ovarian_carcinoma_GSE184880/HGSOC4/bam/SAM_header ${id}_SAM_body > ${id}.sam
samtools view -@ 16 -b ${id}.sam > ${id}.bam
rm ${id}_SAM_body
rm ${id}.sam