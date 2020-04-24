### spaTyper.py: Get spa types


Version: 0.2.0

License: GPLv3

USAGE: python3 spaTyper.py [-h] [-r REPEAT_FILE] [-o REPEAT_ORDER_FILE]
                   [-f FASTA [FASTA ...]] [-g GLOB] 
                   [-e] [-c] [--version]
                   fasta_file.fa

Prints spa type to stdout - egenomics letter combination and then the ridom spa type.

If multiple pcr products are found will print spa types for each product.

Will download sparepeats.fasta and spatypes.txt from the ridom server to repository directory if files not provided or already in directory.
```
optional arguments:
  -h, --help            show this help message and exit
  -r REPEAT_FILE, --repeat_file REPEAT_FILE
                        List of spa repeats
                        (http://spa.ridom.de/dynamic/sparepeats.fasta)
  -o REPEAT_ORDER_FILE, --repeat_order_file REPEAT_ORDER_FILE
                        List spa types and order of repeats
                        (http://spa.ridom.de/dynamic/spatypes.txt)
  -f FASTA [FASTA ...], --fasta FASTA [FASTA ...]
                        List of one or more fasta files.
  -g GLOB, --glob GLOB  Uses unix style pathname expansion to run spa typing
                        on all files. If your shell autoexpands wildcards use
                        -f.
  -e, --do_enrich       Do PCR product enrichment. [Default: False]
  -c, --clean_output    Make output clean
  --version             show program's version number and exit

```

## Installation
Clone the repository or download the python script. In a future, pip download will be available. 

Requires python 3

## How it works

Given a fasta file or multiple fasta files, this script identifies the repeats and the order and generates a 
spa type.

The repeat sequences and repeat orders found on http://spaserver2.ridom.de/ are used to identify the spa type of each enriched sequence.

Ridom spa type and the egenomics repeat sequence are then reported back to the user.

If enriched option provided, the script searches for 50bp to 5000bp sequences produced by the following primer sets
```
TAAAGACGATCCTTCGGTGAG, CAGCAGTAGTGCCGTTTGCTT
AGACGATCCTTCGGTGAGC, GCTTTTGCAATGTCATTTACTG
ATAGCGTGATTTTGCGGTT, CTAAATATAAATAATGTTGTCACTTGGA
CAACGCAATGGTTTCATCCA, GCTTTTGCAATGTCATTTACTG
```

If an enriched sequence is found by a primer set, subsequent primer sets are not used.


## Copyright
Original code written by mjsull.

Jose F. Sanchez-Herrero updated the code, change to python3, and set to use it as a module
