echo
echo "Starting Uppmax jobs ..."
echo

#input data files
#PATH_RAW=                               #path to where raw files are located (top folder)
#PATH_MAIN=                          #path to project folder (my analysis)
RAW_FN=       #file with names of files containing raw reads

N_LINES=$(wc -l < "$RAW_FN") 
N_SAMPLES=$((N_LINES / 2))
                                       #number of samples
echo $N_SAMPLES

declare -A READS                                   #declare variable as array
i=1
j=1

while read -r LINE                                 #reads each filename at a time, stores names (paired) in array 
do
   if [ "$j" -eq 1 ]
   then
      READS[$i,1]=$LINE
      let "j=2"

   else
      READS[$i,2]=$LINE
      let "i+=1"
      let "j=1"
   fi

#   echo $j#${READS[$i,$j]}
done < ${RAW_FN}
#echo $READS	


# prepare output folder
#SRCDIR=$(pwd)                                      #remember current path
cd /proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/1_Trimming/trimmed_files/                               #cd to project folder


#cd $SRCDIR

for i in `seq 1 1 $N_SAMPLES`; do
   RNAFILE_A=${READS[$i,1]}
   RNAFILE_B=${READS[$i,2]}
   echo $RNAFILE_A                                          #display fastq file names (paired)
   echo $RNAFILE_B
   sbatch /proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/1_Trimming/trimgalore_run_query.sh $RNAFILE_A $RNAFILE_B /proj/uppstore2017185/b2014034_nobackup/Aleix/LD_recombination/1_Trimming/trimmed_files/   #starts job
   sleep 1                                                  #pauses for 1 sec
   echo
done


squeue -u $USER                                    #check job status
echo

