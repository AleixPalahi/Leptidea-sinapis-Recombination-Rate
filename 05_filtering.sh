module load bioinfo-tools
module load vcftools/0.1.15

# We set the input and output files we want.

SCAFF=1
VCF_IN=/proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/5_VCF_files/new-scaff-$SCAFF.vcf.gz
VCF_OUT=/proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/6_filtering_SNPs/scaff-$SCAFF-miss0.85.vcf.gz

# Set filters

#MAF=0.02
MISS=0.85
QUAL=30
MIN_DEPTH=8
MAX_DEPTH=100


vcftools --gzvcf $VCF_IN \
--remove-indels \
--mac 2 \
--max-missing $MISS \
--minQ $QUAL \
--min-meanDP $MIN_DEPTH \
--max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH \
--maxDP $MAX_DEPTH \
--recode \
--recode-INFO-all \
--stdout | gzip -c > $VCF_OUT
