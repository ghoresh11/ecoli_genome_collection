
for i in $(seq 10 51); do
if [ $i -eq 50 ]; then
continue
fi
job_name=${i}_combine
bsub -J ${job_name} -G team216 -o ${job_name}.o -e ${job_name}.e -R"select[mem>2000] rusage[mem=2000]" -M2000 python 2.0_build_pairwise_network.py -cluster1 ${i} --cluster2 ${j}
done


for i in $(seq 10 49); do
for j in $(seq $((${i}+1)) 51); do
if [ $j -eq 50 ]; then
continue
fi
job_name=${i}_${j}_combine
bsub -J ${job_name} -G team216 -o ${job_name}.o -e ${job_name}.e -R"select[mem>2000] rusage[mem=2000]" -M2000 python 2.0_build_pairwise_network.py --cluster1 ${i} --cluster2 ${j}
done
done

## combine all the pairwise comparisons to a single network
bsub -J ${job_name} -G team216 -o ${job_name}.o -e ${job_name}.e -R"select[mem>5000] rusage[mem=5000]" -M5000 python 2.1_combine_pairwise_networks.py
