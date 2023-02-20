# We create a variable to indicate the sample we want to work with.

SAMPLE=195

# We calculate the stats on the new generated file.

module load bioinfo-tools
module load GATK/4.1.4.1
module load samtools/1.10

# And run the program to eliminate the duplicate 

gatk MarkDuplicatesSpark \
	-I /proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/2_Mapping_bwa/new_assembly/$SAMPLE.pairs.bam \
	-O /proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/3_duplicates_elimination/new_assembly/${SAMPLE}-duplicates.bam \
	--remove-sequencing-duplicates \
	--create-output-bam-index \
	--create-output-bam-splitting-index

# Finally, we do flagstat on the file without duplicates to check how many have been eliminated.

samtools flagstat /proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/3_duplicates_elimination/new_assembly/$SAMPLE-duplicates.bam > \
	/proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/3_duplicates_elimination/new_assembly/$SAMPLE-duplicates.flagstat
