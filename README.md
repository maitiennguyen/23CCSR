# 23CCSR Tool
This is a tool that aims to streamline the data mining process and facilitate the analysis of biological sequences. 

It automates the process of downloading data from NCBI databases, executing BLAST searches, performing basic annotation, and conducting Clustal analysis. It can be a helpful workflow for exploring and analyzing biological sequence data in the context of evolutionary relationships.

By providing an amino acid or nucleotide query sequence along with a search set, this tool communicates with NCBI databases to download protein and genome datasets for all species in the search set. It performs BLAST searches, including reciprocal BLAST, to identify sequences in the search set that are similar to the query sequence. Additionally, the tool conducts basic annotation to determine start and stop codons for nucleotide sequence hits. For further analysis, protein and nucleotide sequence hits are inputted into Clustal Omega, allowing meaningful multiple sequence alignments for divergent sequences in the search set.


## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
  

## Installation
<a name="installation"></a>
To download the zip file, click [here](http://github.com/maitiennguyen/23CCSR/zipball/master/) or select ‘Download ZIP’ under ‘<>Code’ drop-down options.

After unzipping, make sure to have all files in the working directory.

Before running the script, make sure you have Python on your computer and the following required Python packages and tools installed (using conda environment and conda for installation is recommended) :

#### Python Packages
-	Biopython (https://biopython.org/wiki/Download)

#### Command-Line Tools
-	NCBI BLAST+ (https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata)
-	NCBI datasets and dataformat (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
-	Clustal Omega (http://www.clustal.org/omega/)
-	ETE Toolkit (http://etetoolkit.org)


## Usage
<a name="usage"></a>
#### Required Arguments
- -qseq <query_sequence_file>: Specifies the file containing the query sequence in FASTA format, can be sequence of nucleotides or amino acids. 

- -qdb <query_species_dataset_file>: Specifies the genome or protein dataset file of the query species in FASTA format. If query is an amino acid sequence, this must be a protein dataset. If query is a nucleotide sequence, this must be a genome dataset

- -qname <query_species_name>: Specifies the scientific name of the query species. The input should be in the format 'Genus,species' (e.g., 'Saccharomyces,cerevisiae').

- -qtype <nucl_or_prot>: Specifies the biological sequence type of the query. Use 'nucl' for nucleotide sequence or 'prot' for protein sequence.

- -sset <search_set_taxon_ids>: Specifies the taxon ID(s) of the search set. The ID(s) can be of any taxonomic rank and there can be multiple taxon IDs (e.g, '000' or '000,1111').
  
- -download <yes_or_no>: Specifies whether the search set data has already been downloaded with the script previously. Use 'yes' if the necessary files are in the working directory (see File Documentation for the specifics of the files), or 'no' if the script needs to download the search set data.

#### Optional Arguments
- -evalue <evalue_threshold>: Specifies the E-value threshold for the blast and reciprocal blast searches. If not specified, the default value is '1e-01'.

