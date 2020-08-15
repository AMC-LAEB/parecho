# User instructions for using phyloplace.py  

phyloplace.py is a wrapper script written in Python 3 to perform phylogenetic placement of new query sequences using RAxML-EPA to a phylogenetic tree of reference sequences that were previously clustered by PhyCLIP.  

phyloplace.py categorize the phylogenetically placed query sequences into 3 types:  
 1. Placed on a branch subtending OR within an existing reference cluster.  
 2. Placed on the edge of an outlying unclustered reference sequence.  
 3. Placed on the truck edge between two existing clusters.  

For more details, please check out our paper:  
Parker E., Han A. X., Brouwer L., Wolthers K., Benschop K., Russell C. (2020) Genotypic diversity and dynamic nomenclature of Parechovirus A. bioRxiv 2020.08.14.251231; doi: https://doi.org/10.1101/2020.08.14.251231

## Installation  

phyloplace.py requires the following python libraries to be pre-installed, which can be installed using ```pip```:  

```
pip install -U numpy scipy pandas statsmodels ete3
```

phyloplace.py also needs RAxML to be pre-installed. You can find installation instructions of RAxML here: https://cme.h-its.org/exelixis/web/software/raxml/hands_on.html  

Move the phyloplace.py file in the file directory that you are working in.

## Running phyloplace.py 

phyloplace.py requires minimally the output tree file from PhyCLIP and a sequence alignment FASTA file with all of the reference sequences found in the tree as well as the query sequences to be typed. Note that the headers of the reference sequences must be identical to their corresponding tip names in the phylogenetic tree.  

To run phyloplace.py, open a terminal/command prompt and go to your working directory where the phyloplace.py file had been placed. Type the following minimal command:  

```
phyloplace.py -t <path/to/PhyCLIP_output_tree> -a <path/to/aln_fasta>
```

Make sure that you have changed the permission of phyloplace.py to be executed as a program before:  

```
chmod +x phyloplace.py
```

