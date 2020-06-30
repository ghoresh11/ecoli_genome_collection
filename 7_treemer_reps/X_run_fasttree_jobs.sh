

for d in */ ; do
if [[ $d == *_155* ]]; then
var=${d%?}
echo ${var}
bsub -o ${var}.o -e ${var}.e -J ${var} -n1 -R"span[hosts=1]" bash run_fasttree.sh $d
fi
done


while read d; do
var=${d%?}
echo ${var}
bsub  -R"select[mem>3000] rusage[mem=3000]" -M3000 -o ${var}.o -e ${var}.e -J ${var} -n1 -R"span[hosts=1]" bash run_fasttree.sh $d
done <failed.txt

for d in */ ; do
if [[ $d == 1_1558534808 ]];then
continue
fi
if [[ $d == *_155* ]]; then
var=${d%?}
echo ${var}
bsub -o ${var}.o -e ${var}.e -J ${var} -n1 -R"span[hosts=1]" bash run_fasttree.sh $d
fi
done
