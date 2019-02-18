


while read -r line
do
    echo $line
    fastq-dump.2.9.2 --gzip --split-files -O downloaded -A $line
done < "${1}"


