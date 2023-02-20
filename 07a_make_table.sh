#We create a variable for the scaffold.
SCAFF=29

#We activate the virtual environment in which we execute pyrho.
source pyrho-virtual-env/bin/activate

#We create the table
pyrho make_table -n 168 -N 210 --mu 2.9e-9 --logfile scaff-$SCAFF.out --outfile scaff-29_lookuptable.hdf \
	--approx --smcpp_file /proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/9_SMC++/scaff-$SCAFF-reduced.csv \
	--numthreads 20 --decimate_rel_tol 0.1
