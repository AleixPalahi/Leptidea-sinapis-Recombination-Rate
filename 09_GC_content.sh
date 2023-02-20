module load bioinfo-tools
module load BEDTools/2.26.0

bedtools nuc -fi P14502_103.FINAL-deduped-nuc.filtered.fasta -bed right_flanks_hotspots.bed > GC_right_flanks.txt
