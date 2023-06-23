# 23CCSR Tool
This is a tool that aims to streamline the data mining process and facilitate the analysis of biological sequences. 

It automates the process of downloading data from NCBI databases, executing BLAST searches, performing basic annotation, and conducting Clustal analysis. It can be a helpful workflow for exploring and analyzing biological sequence data in the context of evolutionary relationships.

By providing an amino acid or nucleotide query sequence along with a search set, this tool communicates with NCBI databases to download protein and genome datasets for all species in the search set. It performs BLAST searches, including reciprocal BLAST, to identify sequences in the search set that are similar to the query sequence. Additionally, the tool conducts basic annotation to determine start and stop codons for nucleotide sequence hits. For further analysis, protein and nucleotide sequence hits are inputted into Clustal Omega, allowing meaningful multiple sequence alignments for divergent sequences in the search set.


## Table of Contents
- [Workflow](#workflow)
- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [File Documentation](#filedoc)
- [Troubleshooting](#troubleshoot)


## Workflow
<a name="workflow"></a>

Coming Soon...

## Installation
<a name="installation"></a>

To download the zip file, click [here](http://github.com/maitiennguyen/23CCSR/zipball/master/) or select 'Download ZIP' under '<>Code' drop-down options.

After unzipping, make sure to have all files in the working directory.

Before running the script, make sure you have Python on your computer and the following required Python packages and tools installed (using conda environment and conda for installation is recommended) :

### Python Packages
-	Biopython (https://biopython.org/wiki/Download)

### Command-Line Tools
-	NCBI BLAST+ (https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata)
-	NCBI datasets and dataformat (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
-	Clustal Omega (http://www.clustal.org/omega/)
-	ETE Toolkit (http://etetoolkit.org)


## Usage
<a name="usage"></a>

### Required Arguments
```
-qseq <query_sequence_file>
  ```
Specifies the file containing the query sequence in FASTA format, can be sequence of nucleotides or amino acids. If file is not in the working directory, path of the file is required.

```
-qdb <query_species_dataset_file>
```
Specifies the genome or protein dataset file of the query species in FASTA format. If query is an amino acid sequence, this must be a protein dataset. If query is a nucleotide sequence, this must be a genome dataset. If file is not in the working directory, path of the file is required.

```
-qname <query_species_name>
```
Specifies the scientific name of the query species. The input should be in the format 'Genus,species' (e.g., 'Saccharomyces,cerevisiae').

```
-qtype <nucl_or_prot>
```
Specifies the biological sequence type of the query. Use 'nucl' for nucleotide sequence or 'prot' for protein sequence.

```
-sset <search_set_taxon_ids>
```
Specifies the taxon ID(s) of the search set. The ID(s) can be of any taxonomic rank and there can be multiple taxon IDs (e.g, '000' or '000,1111').

```
-download <yes_or_no>
```
Specifies whether the search set data has already been downloaded with the script previously. Use 'yes' if the necessary files are in the working directory (see File Documentation for the specifics of the files), or 'no' if the script needs to download the search set data.

### Optional Arguments
```
-evalue <evalue_threshold>
```
Specifies the E-value threshold for the blast and reciprocal blast searches. If not specified, the default value is '1e-01'.

### Add-on Command
```
run_clustal.py <fasta_file>
```
This command can be used to perform Clustal analysis of sequences in FASTA format individually.


## Examples
<a name="examples"></a>

```
python main.py -qseq cnp20_spombe.fasta -qdb spombe.faa -qname Schizosaccharomyces,pombe -qtype prot -sset 147537 -download no 
```

```
python main.py -qseq sir2_gene.fna -qdb scereviseae.fna -qname Saccharomyces,cerevisiae -qtype nucl -sset 4930,34365,2916678 -download yes -evalue 1e-50
```

```
python run_clustal.py manual_anno_seqs.fasta
```


## File Documentation
<a name="filedoc"></a>

### Python Files

The script for the tool consists of the following files:

- `main.py`: This is the main script that should be run to execute the tool. It integrates all the other Python files and orchestrates the workflow.

- `user_input.py`: This contains methods to process user terminal arguments and handle any argument errors. It ensures that the user inputs are correctly parsed and validated before further processing.

- `blast.py`: This contains methods to perform BLAST (Basic Local Alignment Search Tool) and reciprocal BLAST searches. It also handles the processing and validation of the search results.

- `annotation.py`: This contains methods to sort the BLAST results into two lists: one that requires manual annotation and another that can undergo automated annotation. This file helps streamline the annotation process by categorizing the results accordingly.

- `clustal.py`: This contains methods to perform Clustal analysis. It facilitates multiple sequence alignment, allowing for the creation of biologically meaningful alignments of divergent sequences in the search set.

- `process_specs.py`: This contains methods to sort species by family. It aids in the organization and classification of species based on their taxonomic families.

- `run_clustal.py`: This script is an add-on tool and can be used if manual Clustal analysis is required.

### Output Files

After running the tool, the following output files will be generated:

- `nucl.fna`: The genome dataset of all species in the search set in FASTA format.

- `prot.faa`: The protein dataset of all species in the search set in FASTA format.

- `all_spec.txt`: A text file containing the names of all species in the search set.

- `prot_data_specs.txt`: A text file containing the names of species in the search set with a protein dataset.

- `prot_files_all_dict.txt`: A text file containing the file paths to the protein files in FASTA format for each species with a protein dataset in the search set.

- Folders named after the input taxon IDs: Each folder contains the genome and protein dataset files for each species in that taxon.

- `[blast_type]1_[#]_dict.txt`: A text file containing the BLAST results in Python dictionary format before any validation. The `#` represents the query round number, and `blast_type` can be either `blastp` or `tblastx` depending on the query sequence type.

- `[blast_type]2_[#]_dict.txt`: A text file containing the BLAST results in Python dictionary format after being validated by reciprocal BLAST. The `#` represents the query round number, and `blast_type` can be either `blastp` or `tblastx` depending on the query sequence type.

- `[blast_type]_[#]_dict.txt`: A text file containing the BLAST results in Python dictionary format before or after reciprocal BLAST. The `#` represents the query round number, and `blast_type` can be either `tblastn` or `blastx`.

- `[blast_type]_[#]_summary_report.txt`: A text file containing a summary of the valid BLAST hits. The `#` represents the query round number.

- `complete_blast_summary_report.txt`: A text file containing an overview summary of the BLAST results from all query rounds.

- `annotated_[species_name]_dict.txt`: A text file containing information about annotated nucleotide hits for a specific query species and round. The `[species_name]` represents the query species for that round.

- `[species_name]_man_anno_seqs.txt`: A text file containing information about nucleotide hit sequences that need to be manually annotated and excluded from the script's Clustal analysis for a specific query species and round.

- `all_man_anno.txt`: A text file containing information about all nucleotide hit sequences that require manual annotation from all query rounds.

- `complete_manual_anno_summary.txt`: A text file containing a summary overview of all sequences that require manual annotation from all query rounds.

- `complete_auto_anno_summary.txt`: A text file containing a summary overview of all sequences that have been automatically annotated by the script from all query rounds.

- `auto_anno_seqs.fasta`: A file containing all script-annotated nucleotide sequences in FASTA format from all query rounds.

- `auto_algn.clustal`: The Clustal analysis results of all protein and nucleotide hit alignments from all query rounds in Clustal format.

- `auto_algn.fasta`: The Clustal analysis results of all protein and nucleotide hit alignments from all query rounds in FASTA format.

These output files provide valuable information and results from the tool's data mining and analysis workflow, any files not mentioned are miscellaneous and can be discarded.

### Required Files for Re-using Search Set

If you want to re-run the script and the search set remains the same as the previous run (when the argument for -download is 'no'), then the following files are required to be in the working directory for the script to run:

- `nucl.fna`
- `prot.faa`
- `all_spec.txt`
- `prot_data_specs.txt`
- `prot_files_all_dict.txt`
- All folders named after the input taxon IDs



## Troubleshooting
<a name="troubleshoot"></a>

If you encounter any issues while using the tool, you can refer to the following troubleshooting steps:

- **Error message indicating duplicate protein or nucleotide sequences when making BLAST database**
  - This usually occurs when the script attempts to download the search set data, but there are existing 'prot.faa' and 'nucl.fna' files from a previous run of the script in the working directory. To resolve this, delete the existing files and run the script again.

- **Error message indicating invalid amino acid or nucleotide when making BLAST database**
  - This error can occur due to a poor internet connection, causing interruptions during the download process from NCBI databases to your local computer. This causes the FASTA file to be formatted incorrectly. To resolve this, try running the script again with a better internet connection and ensure the 'prot.faa' and 'nucl.fna' files are not broken (delete them if necessary).

- **Required search set database and/or files not found in the directory**
  - Make sure all the required files are present in the working directory before running the script. These files include:
    - `nucl.fna`
    - `prot.faa`
    - `all_spec.txt`
    - `prot_data_specs.txt`
    - `prot_files_all_dict.txt`
    - All folders named after the input taxon IDs.

- **No BLAST result text files generated**
  - If there are no blastp results for a protein query or no blastx/tblastn results for a nucleotide query, it could be due to either none of the species in the search set having a protein dataset. In any blast types, another reason could be there being no blast hits for the given query.

- **Query sequence file and/or Query nucleotide/protein dataset file not in valid format**
  - Ensure that the query sequence file and/or query nucleotide/protein dataset are in FASTA format and contain valid amino acid or nucleotide bases. Check for any formatting errors or invalid characters in the files.

