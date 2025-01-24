##  Program Name
clustalo_alignment_using_multiple_threads.rb

## Requirements
1. Ruby 3.0 or higher
2. Dependencies: Ensure Clustal Omega is installed and available in your system PATH.

## Usage
   clustalo_alignment_using_multiple_threads.rb -o OUTPUT_DIRECTORY bacteria1.fasta bacteria2.fasta ...

## Description
This script utilizes Clustal Omega to create consensus sequences from 16S ribosomal RNA sequence data for specified bacteria. It supports multi-threaded processing to handle multiple sequence alignments efficiently.

## Options
`-o`: (Required) Specifies the output directory where output files for alignment and consensus sequences are stored.

## Parameters
`bacteria1.fasta bacteria2.fasta ...` : Input FASTA files containing 16S rRNA sequences for the target bacteria.

## Output
The alignment results for each bacterium are stored in the specified directory as genus_species_number_Clustalo.fasta. The consensus sequences for each bacterium are compiled and output as 16S_rRNA_consensus_for_each_bacteria.txt. The consensus sequence for all bacteria, obtained by further aligning the consensus sequences of each bacterium, is output as 16S_rRNA_consensus_for_bacterias.txt

## License
This script is distributed under the MIT License. Feel free to use and modify it as needed.

## Contact
For any questions or issues, please contact:
**Koichi Hagiwara**  
Email: hagiwark@me.com
