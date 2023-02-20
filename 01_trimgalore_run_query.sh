# load modules
module load bioinfo-tools
module load TrimGalore/0.4.1
module load cutadapt/1.9.1
module load FastQC/0.11.5
#module load MultiQC/0.9



READ1=${1:?msg}                              #read name of files containing paired reads
READ2=${2:?msg}
SRCDIR=                            #location of output folder


#cd $SNIC_TMP                                # use local scratch disk


# USAGE: trim_galore [options] <filename(s)>

trim_galore \
            --quality 30 \
            --paired \
            --illumina \
            --phred33 \
            --stringency 1 \
            -e 0.1 \
            --clip_R1 12 \
            --clip_R2 12 \
            --three_prime_clip_R1 3 \
            --three_prime_clip_R2 3 \
            --length 30 \
            --gzip \
            --output_dir $SRCDIR \
            --fastqc \
            $READ1 \
            $READ2


#multiqc $SRCDIR/                  # creates single report

