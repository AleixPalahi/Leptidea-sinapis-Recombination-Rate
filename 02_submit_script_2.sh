## CONTENTS ##
RAW_FN=/proj/uppstore2017185/b2014034/nobackup/Aleix/LD_recombination/2_Mapping_bwa/sample_list_2.txt        #file with names of files co$

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


# navigate to output folder
cd /proj/uppstore2017185/b2014034/nobackup/Aleix/LD_recombination/2_Mapping_bwa/new_assembly                               #cd to project folder


for i in `seq 1 1 $N_SAMPLES`; do
   DNAFILE_A=${READS[$i,1]}
   DNAFILE_B=${READS[$i,2]}
   echo $DNAFILE_A                                          #display fastq file names (paired)
   echo $DNAFILE_B
   sbatch /proj/uppstore2017185/b2014034/nobackup/Aleix/LD_recombination/2_Mapping_bwa/run_script_2.sh $DNAFILE_A  $DNAFILE_B   #sta$
   sleep 1                                                  #pauses for 1 sec
   echo
done

