# mign
mign is a tool that aims to streamline the data mining process and facilitate the analysis of biological sequences. 

mign automates the process of downloading data from NCBI databases, executing BLAST searches, performing basic annotation, and conducting Clustal analysis. It can be a helpful workflow for exploring and analyzing biological sequence data in the context of evolutionary relationships.

By providing an amino acid or nucleotide query sequence along with a search set, mign communicates with NCBI databases to download protein and genome datasets for all species in the search set. It performs BLAST searches, including reciprocal BLAST, to identify sequences in the search set that are similar to the query sequence. Additionally, mign conducts basic annotation to determine the start and stop codons for nucleotide sequence hits. For further analysis, protein and nucleotide sequence hits are inputted into Clustal Omega, allowing meaningful multiple sequence alignments for divergent sequences in the search set.


## Table of Contents
- [Workflow](#workflow)
- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [File Documentation](#filedoc)
- [Re-using Search Set](#rerun)
- [Troubleshooting](#troubleshoot)


## Workflow
<a name="workflow"></a>

#### Homologous Sequences Search using BLAST+:
- Users provide a biological sequence of interest and a species search set.
- mign performs a BLAST+ search using the provided biological sequence as the query and the species' biological dataset the search set.
- Homologous sequences to the query are identified based on sequence similarity.

#### Protein Dataset Search:
- Species with annotated genome in the NCBI database are added to the protein search set.
- The biological sequence of interest is used as the query in BLAST+ search to identify homologous sequences in the protein dataset.

#### Reciprocal Best Hit Validation (Protein):
- The tool performs reciprocal BLAST searches for each identified homologous sequence.
- Only sequences where the initial biological sequence of interest is the top hit in the reciprocal BLAST results are considered valid and included in the final set of homologous sequences.

#### Nucleotide Dataset Search:
- Species without annotated genome or without hits from the protein BLAST+ search are added to the nucleotide search set.
- The biological sequence of interest is used as the query in another BLAST+ search to identify homologous sequences in the nucleotide dataset.

#### Reciprocal Best Hit Validation (Nucleotide):
- Reciprocal BLAST searches are performed to confirm the identified nucleotide homologous sequences.
- Only sequences where the initial biological sequence of interest is the top hit in the reciprocal BLAST results are considered valid and included in the final set of homologous sequences.

#### Iterative Search:
- The tool selects a random species in the protein hit results, the species must be from a different family rank than the query species', and uses its homologous sequence as a new query.
- The search for homologous sequences continues iteratively until no new sequences are found, there are no species in other family ranks to select from, or there are no species left in the protein search results to select from.

#### Annotation:
- Nucleotide hits requiring automated annotation and those needing manual annotation are identified.
- For sequences with more than one non-overlapping hit at the same scaffold, manual annotation is required, and the user is provided with an output file to work with.
- Automated annotation of homologous sequences is performed for the remaining hits.
- Start and stop codons are identified based on proximity to the sequence ends and in consideration of the scaffold's start and end positions.
- For the stop codon, it searches for the closest occurrence in the 3' direction of the sequence.
- For the start codon, it searches for the nearest stop codon in the 5' direction of the sequence and proceeds towards the 3' direction until it encounters the first start codon.
- If no stop codon and/or start codon is found, the start and/or end of the scaffold becomes the start and end of the sequence.
- The annotated nucleotide sequences are translated into corresponding amino acid sequences.

#### Multiple Sequence Alignment using Clustal Omega:
- The protein hits and the translated nucleotide sequences with automated annotation are compiled into a FASTA file.
- The FASTA file is inputted into Clustal Omega for multiple sequence alignment.


## Installation
<a name="installation"></a>

To download the zip file, click [here](http://github.com/maitiennguyen/23CCSR/zipball/master/) or select 'Download ZIP' under '<>Code' drop-down options.

After unzipping, make sure to have all files in the working directory.

Before running the script, make sure you have Python on your computer and the following required Python packages and tools installed (using conda environment and conda for installation is recommended):

### Python Version
- 3.10 or later

### Python Packages
-	Biopython (https://anaconda.org/anaconda/biopython)

### Command-Line Tools
-	NCBI BLAST+ (https://anaconda.org/bioconda/blast)
    -	If running on Apple M1/M2, NCBI BLAST+ should be installed using source code (https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
-	NCBI datasets and dataformat (https://anaconda.org/conda-forge/ncbi-datasets-cli)
-	Clustal Omega (https://anaconda.org/bioconda/clustalo)
    -	If running on Apple M1/M2, Clustal Omega should be installed using precompiled binary (http://www.clustal.org/omega/)
-	ETE Toolkit (https://anaconda.org/bioconda/ete3)
    -	If installation fails, use pip:
    	-	```pip install six```
        -	```pip install ete3```



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
Specifies whether the search set data has already been downloaded with the script previously. Use 'yes' if the necessary files are in the working directory (see [Re-using Search Set](#rerun) for the specifics of the files), or 'no' if the script needs to download the search set data.

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
python main.py -qseq sir2_gene.fna -qdb scerevisiae.fna -qname Saccharomyces,cerevisiae -qtype nucl -sset 4930,34365,2916678 -download yes -evalue 1e-50
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

- `organize_files.py`: This contains methods to organize output files into folders and remove miscellaneous files.

- `run_clustal.py`: This script is an add-on tool and can be used if manual Clustal analysis is required.

### Output Files

After running the tool, the following output folders and files will be generated:

- `InputFiles`: A folder containing the user input files.

  
- `SearchSetFiles`: A folder containing files that are required for re-using search set for future runs.

    - `nucl.fna`: The genome dataset of all species in the search set in FASTA format.
    
    - `prot.faa`: The protein dataset of all species in the search set in FASTA format.
    
    - `all_spec.txt`: A text file containing the names of all species in the search set.
    
    - `prot_data_specs.txt`: A text file containing the names of species in the search set with a protein dataset.
    
    - `prot_files_all_dict.txt`: A text file containing the file paths to the protein files in FASTA format for each species with a protein dataset in the search set.
    
    - Folders named after the input taxon IDs: Each folder contains the genome and protein dataset files for each species in that taxon.

 
- `MainFiles`: A folder containing primary output files.

    - `complete_blast_summary_report.txt`: A text file containing an overview summary of the BLAST results from all query rounds.
 
    - `complete_manual_anno_summary.txt`: A text file containing a summary overview of all sequences that require manual annotation from all query rounds.

    - `complete_auto_anno_summary.txt`: A text file containing a summary overview of all sequences that have been automatically annotated by the script from all query rounds.
   
    - `auto_algn.clustal`: The Clustal analysis results of all protein and nucleotide hit alignments from all query rounds in Clustal format.

    - `auto_algn.fasta`: The Clustal analysis results of all protein and nucleotide hit alignments from all query rounds in FASTA format.


- `DictReports`: A folder containing BLAST results in Python dictionary format.

    - `[blast_type]1_[#]_dict.txt`: A text file containing the BLAST results in Python dictionary format before any validation. The `#` represents the query round number, and `blast_type` can be either `blastp` or `tblastx` depending on the query sequence type.
    
    - `[blast_type]2_[#]_dict.txt`: A text file containing the BLAST results in Python dictionary format after being validated by reciprocal BLAST. The `#` represents the query round number, and `blast_type` can be either `blastp` or `tblastx` depending on the query sequence type.
    
    - `[blast_type]_[#]_dict.txt`: A text file containing the BLAST results in Python dictionary format before or after reciprocal BLAST. The `#` represents the query round number, and `blast_type` can be either `tblastn` or `blastx`.

 
- `SummaryReports`: A folder containing summary reports of the valid BLAST hits.

    - `[blast_type]_[#]_summary_report.txt`: A text file containing a summary of the valid BLAST hits. The `#` represents the query round number.


- `AnnotationFiles`: A folder containing files pertaining auto and manual annotation.
  
    - `annotated_[species_name]_dict.txt`: A text file containing information about annotated nucleotide hits for a specific query species and round. The `[species_name]` represents the query species for that round.
    
    - `[species_name]_man_anno_seqs.txt`: A text file containing information about nucleotide hit sequences that need to be manually annotated and excluded from the script's Clustal analysis for a specific query species and round.

    - `all_man_anno.txt`: A text file containing information about all nucleotide hit sequences that require manual annotation from all query rounds.
      
    - `auto_anno_seqs.fasta`: A file containing all script-annotated nucleotide sequences in FASTA format from all query rounds.


- Miscellaneous folders and files:
    - `QuerySeqs`: A folder containing query sequences in fasta format.


## Re-using Search Set
<a name="rerun"></a>
If you want to re-run the script and the search set remains the same as the previous run (when the argument for -download is 'no'), then the following files, which can be found in the `SearchSetFiles` folder, are **required to be moved into the working directory** for the script to run:

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
  - Make sure all the required files are present in the working directory before running the script. These files can be found in the `SearchSetFiles` folder, which include:
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

