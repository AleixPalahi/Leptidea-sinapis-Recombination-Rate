#We activate the virtual environment in which we execute pyrho.
source pyrho-virtual-env/bin/activate

#We test the parameter values
pyrho hyperparam -n 168 --mu 2.9e-9 --blockpenalty 10,25,50,100,150,200 \
	--windowsize 10,25,50,75,100 --logfile . --tablefile scaff-29_lookuptable.hdf \
	--num_sims 4 --ploidy 2 --numthreads 20 \
	--smcpp_file ../9_SMC++/scaff-29-reduced.csv --outfile scaff-29-full_hyperparam_results.txt
