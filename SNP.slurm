#!/bin/bash

#SBATCH --job-name=SNP
#SBATCH --partition=64c512g
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=wcj342403669@sjtu.edu.cn
#SBATCH -n 10
#SBATCH --array=1-1990

module load miniconda3
source activate mypython3
module load samtools
id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /dssg/home/acct-medcl1/medcl1/Ovarian_carcinoma_GSE184880/HGSOC4/bam/xaa)
samtools mpileup -gSDf /dssg/home/acct-medcl1/medcl1/opt/refdata-gex-GRCh38-2020-A/fasta/genome.fa /dssg/home/acct-medcl1/medcl1/Ovarian_carcinoma_GSE184880/HGSOC4/bam/bam/${id}.bam > ${id}.bcf
bcftools call -vmO z -o ${id}.vcf.gz ${id}.bcf
bcftools filter -O v -o ${id}_filter.vcf -s LOWQUAL -e 'QUAL<10 || FMT/DP <5' --SnpGap 5 --set-GTs . ${id}.vcf.gz
bcftools view -v snps ${id}_filter.vcf > ${id}_snp_filter.vcf
sed -i '1,230d' ${id}_snp_filter.vcf
awk '{print $0"\t'${id}'"}' ${id}_snp_filter.vcf > ${id}_SNP_filter.vcf
rm ${id}.bcf
rm ${id}.vcf.gz
rm ${id}_filter.vcf
rm ${id}_snp_filter.vcf

