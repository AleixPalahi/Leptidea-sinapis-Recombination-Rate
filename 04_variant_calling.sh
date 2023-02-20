# We create a variable to indicate the sample we want to work with, and another for the path to the reference genome and the scaffold.

REFERENCE=/proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/4_SNP_calling/P14502_103.FINAL-deduped-nuc.filtered.fasta
SAMPLE=196
SCAFFOLD=29

# We load the required packages

module load bioinfo-tools
module load GATK/4.1.4.1

# And run the program to eliminate the duplicate 
gatk --java-options "-Xmx6g" HaplotypeCaller \
	-R $REFERENCE \
	-I /proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/3_duplicates_elimination/new_assembly/$SAMPLE-duplicates.bam \
	-O /proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/4_SNP_calling/new_assembly/new_scaffold_$SCAFFOLD/$SAMPLE-scaff-$SCAFFOLD.g.vcf.gz \
	-ERC GVCF \
	-L HiC_scaffold_$SCAFFOLD \
	--add-output-vcf-command-line false
