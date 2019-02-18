
### run the trimming scripts.
### usage ./run_trimming read_input_dir output_dir batch


batch=${3}

let START=84*$(expr ${batch} - 1)+1 ## change the 84 here depending how many genomes you had in each batch 
let STOP=84*$(expr ${batch} - 1)+84

## change here depending what batch I'm on
for (( c=$START; c<=$STOP; c++ ))
do
	echo $c
	arg3=${c}
	job_name=trimming_${c}
	bsub -J ${job_name} -G team216 -o ${job_name}.o -e ${job_name}.e -R"select[mem>2000] rusage[mem=2000]" -M2000 -n5 -R"span[hosts=1]" python3.6 X_trimming.py ${1} ${2} ${arg3}
done

exit

