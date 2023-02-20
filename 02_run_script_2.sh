## load modules
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.10

READ_1=${1:?msg}                              #read name of files containing paired reads
READ_2=${2:?msg}
SAMPLE=196
#LOCATION=SWE
REFERENCE=

## commands:
# Go from fastQ to BAM #

bwa mem -t 19 -M $REFERENCE -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" $READ_1 $READ_2 > $SAMPLE.sam

samtools flagstat $SAMPLE.sam > $SAMPLE.sam.flagstat &
samtools view -Sb -f 3 $SAMPLE.sam > $SAMPLE.pairs.bam
samtools flagstat $SAMPLE.pairs.bam > $SAMPLE.pairs.bam.flagstat &
#samtools sort $SAMPLE.pairs.bam > $SAMPLE.pairs.sorted.bam
#samtools index $SAMPLE.pairs.sorted.bam

