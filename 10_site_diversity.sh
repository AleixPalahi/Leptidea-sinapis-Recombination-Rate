module load bioinfo-tools
module load vcftools/0.1.15
zcat /proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/10c_pyrho/28-filtered.vcf.gz | \
 vcftools --vcf - --site-pi --positions SNPs_chr_28_masked.txt --out chr_28_per_site_masked
