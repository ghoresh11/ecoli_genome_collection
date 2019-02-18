
# bin_genomes
Bin and QC a collection of genome sequences

usage: `python bin_genomes.py [options] <job_id> <genomes_file>`

Bin input genomes according to sequence identity using specified method.

### Positional arguments:
  `<job_id> STR `                  Name used to describe this run
 
  `<genomes_file> FILE`               File with list of genome files (Default format: FASTA, set `--gff` for GFF)

### Optional arguments:

  `-h, --help`            show this help message and exit
  
  `--gff`                 Set if input in GFF format [Default: FASTA]
  
  `--species_cutoff FLOAT` Maximum distance between species to remove contaminents [0.04]
  
  `--distance FLOAT `     Maximum distance between two genomes to be considered in same bin [0.005]
  
  `--max_contigs INT`   Skip genomes with more than num_contigs [600]
  
  `--min_length FLOAT` Skip genomes shorter than this length, in MBP [4]
  
  `--max_length FLOAT` Skip genomes longer than this length, in MBP [6]
  
  `--cpu INT`             Number of CPUs to use [16]
  
  `--keep_temp`           Keep temporary files
  
  `--verbose`             Verbose output while run
  
  `--debug `              Set for Debug mode (doesnt run MASH)
  
  `--method STR `         Method to bin the genomes. Options: (MASH), [MASH]
  
  
# create_mash_tree

Create a minimum spanning tree of the MASH clusters

usage: `python create_mash_tree.py <bin_file> <distances_file>`

### Positional arguments:

 ` <bin_file> FILE   `     Output of `bin_genomes.py` (`[job_id]_binned.txt`)
  
  ` <distances_file>FILE    `    MASH distances output from `bin_genomes.py` (`[job_id]_distances.tab`)

### Optional arguments:

  `-h, --help`  show this help message and exit