
cd ${1}

snp_sites core_gene_alignment.aln -m -o core_gene_alignment_snps.aln

FastTreeMP -mlnni 4 -nt core_gene_alignment_snps.aln > tree

## This will run treemer and keep only 10 leaves
runjob5  -e more_cpus2.e -o more_cpus2.o -n10 -R"span[hosts=1]" python ~/Treemmer/Treemmer.py -r 4 --cpu 10 -np -X 10 tree
