SCAFFOLD=19

cd /proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/9_SMC++/scaff-$SCAFFOLD/definitive
singularity exec /proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/9_SMC++/scaff-$SCAFFOLD/definitive/smcpp_latest.sif echo "starting SMC++"
singularity run smcpp_latest.sif plot -c scaff-$SCAFFOLD-definitive.png /proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/9_SMC++/scaff-$SCAFFOLD/definitive/model.final.json
