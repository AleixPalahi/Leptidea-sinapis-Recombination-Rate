#We activate the virtual environment in which we execute pyrho.
source pyrho-virtual-env/bin/activate

SCAFF=17

#We run optimize to obtain the recombination landscape for a chromosome.
pyrho optimize --tablefile scaff-27_lookuptable.hdf \
	--vcffile /proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/16_sims/Recomb_simulation_vcfs/concat_76.vcf \
	--outfile sim-long-trial.rmap \
	--ploidy 2 --numthreads 20 \
	--blockpenalty 25 --windowsize 100 \
	--fast_missing --logfile .
