## Program Name
cheryblast_plus_ob.rb

## Requirements
1. Ruby 3.0 or higher (https://www.ruby-lang.org).
2. BLAST Tools:
   - `blastn` and `makeblastdb` from the `ncbi-blast-2.8.1+` package.
   - Download from: [NCBI BLAST+ Executables](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/)
   - Set the tools in `/blast_directory/bin/` on your local computer as:
     - `/blast_directory/bin/blastn`
     - `/blast_directory/bin/makeblastdb`

## Usage
For FASTQ Files:
cheryblast_plus_ob.rb -o OUTPUT_DIRECTORY -b BLAST_DIRECTORY -d database 16S_rRNA_file1.fasta 16S_ribosomal_RNA_file2.fasta ...

For Pair-End Read FASTQ Files:
cheryblast_plus_ob.rb -o OUTPUT_DIRECTORY -b BLAST_DIRECTORY -d database 16S_ribosomal_RNA_R1_file1.fastq 16S_ribosomal_RNA_R2_file1.fastq 16S_ribosomal_RNA_R1_file2.fastq 16S_ribosomal_RNA_R2_file2.fastq ...

FASTA files and pair-end read FASTQ files can be processed together.

## Description
This script utilizes BLAST (blastn and makeblastdb) to perform homology searches between input sequences (FASTA or FASTQ) and a local database. It processes sequences to classify them based on bit score thresholds and provides detailed and summary outputs.

- Each sequence (read) in the input file (`file1`) is used as a query.
- Each entry in the local database (`dbase`) is used as a subject for homology searches using `blastn`.
- Reads with bit scores exceeding the threshold are classified and recorded in `file1_local_dbase_result.txt`.

## Output
### (A) `file1_local_dbase_result.txt`
This file contains:
1. **Summary Section:**
   - Example:
     ```
     Streptococcus_oralis_Region_I	48
     Undetermined	97
     ```
     This shows that among 145 reads:
     - 48 reads were classified as `Streptococcus_oralis_Region_I`.
     - 97 reads were classified as `Undetermined` (bit score below threshold).

2. **Details for Each Read:**
   - Example:
     
     ID	Diagnosis	bit_score
     Read1	Streptococcus_oralis_Region_I	907
     
     - `Read1` has a bit score of 907 for `Streptococcus_oralis_Region_I`.

   - For undetermined reads:
     
     Read2	Undetermined		Streptococcus_oralis_Region_I	902	...
     
     - `Read2` did not exceed the threshold for any database entry. Detailed scores for each match are listed.

Threshold values for database entries are specified in line 8 of `cheryblast_plus_ob.rb`.

### (B) `summary.txt`
- Example:
  
  16S_rRNA_file1_local_database_result.txt	Streptococcus_oralis_Region_I	48	Undetermined	97
  
- When processing multiple input files, this file summarizes their results for easy reference.

## License
This script is distributed under the MIT License. Feel free to use and modify it as needed.

## Contact
For any questions or issues, please contact:
**Koichi Hagiwara**  
Email: hagiwark@me.com
