## Program Name
download_16S_rRNA_data_from_genbank.rb

## Requirements
1. Ruby 3.0 or higher (https://www.ruby-lang.org).
2. Ruby gem `nokogiri` must be installed. If it is not installed, you can install it using the following command:

   gem install nokogiri

## Usage
   download_16S_rRNA_data_from_genbank.rb -m MAX_NUMBER_RETRIEVED -o OUTPUT_DIRECTORY bacteria1 bacteria2 ...

## Description
This script downloads 16S ribosomal RNA sequence data for specified bacteria from whole genome sequences available in GenBank. Each bacterium should be specified using its Genus_species name, with the genus and species names separated by an underscore (`_`).

## Options
`-m`: (Required) Specifies the maximum number of sequences to retrieve per bacterium.
`-o`: (Required) Specifies the output directory where the retrieved data will be saved.

## Example
To retrieve data for `Escherichia coli` and `Staphylococcus aureus`, limiting the maximum number of sequences to 100 and saving the output in a directory named `result_dir`, use the following command:

   download_16S_rRNA_data_from_genbank.rb -m 100 -o result_dir Escherichia_coli Staphylococcus_aureus ...

## Output
The output will be saved in the specified directory with the following naming convention:

outputDirectory/_genus_species_from_genbank_numberOfEntriesRetrieved.txt

Each file will contain the retrieved sequence data for the specified bacterium, with up to the maximum number of entries defined by the `-m` option.

## License  
This program is released under the MIT License.  

## Contact  
**Koichi Hagiwara**  
- E-mail: hagiwark@me.com  
- GitHub: [@hagiwark1957](https://github.com/hagiwark1957)  
