#!/usr/bin/env ruby
require 'net/http'
require 'nokogiri'
require 'optparse'

# Usage description and help message
def print_usage_and_exit
  puts <<~USAGE
    \nUsage: download_16S_rRNA_data_from_genbank.rb -m MAX_NUMBER_RETRIEVED -o OUTPUT_DIRECTORY bacteria1 bacteria2, ...
    Description:
      This script downloads 16S ribosomal RNA sequence data for specified bacteria from GenBank.
      Each bacterium should be specified by its Genus_species name separated by an underscore.
      The '-m' option is required and specifies the maximum number of sequences to retrieve per bacterium.
      The '-o' option is required and specifies the output directory.
    Example:
      download_16S_rRNA_data_from_genbank.rb -m 100 -o result_dir Escherichia_coli Staphylococcus_aureus, ...\n
  USAGE
  exit(1)
end

def rc(seq)
  seq.upcase.reverse.tr('ACGTNWRMKYSHBVD', 'TGCANSYKMRWDVBH')
end

#NCBI API endpoint of data search.
SEARCH_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
#NCBI API endpoint of data retrieval.
FETCH_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Initialize options hash to store command line arguments
options = {}

# Method to make API requests to NCBI.
def fetch_from_ncbi(base_url, params)
  uri = URI(base_url)
  uri.query = URI.encode_www_form(params)
  response = Net::HTTP.get_response(uri)
  sleep(0.5)  # Throttling to comply with API usage guidelines
  raise "HTTP Error: #{response.code}" unless response.is_a?(Net::HTTPSuccess)
  response.body
end

def fetch_rRNA_sequence(base_url, accession_number, start, stop)
  params = {
    'db' => 'nuccore',
    'id' => accession_number,
    'rettype' => 'fasta',
    'retmode' => 'text',
    'seq_start' => [start, stop].min,
    'seq_stop' => [start, stop].max
  }
  fetch_from_ncbi(base_url, params)
end

# Method to extract 16S rRNA start and end positions.
def find_16S_rRNA_positions(feature_lines)
  positions = []
  feature_lines.each_cons(2) do |line, next_line|
    if line =~ /rRNA/ && next_line.include?("16S ribosomal RNA") && line =~ /^(\d+)\t(\d+)/
      positions << [$1.to_i, $2.to_i]
    end
  end
  positions
end

#Parse command line options
OptionParser.new do |opts|
  opts.banner = "Usage: download_16S_rRNA_data_from_genbank.rb -m max_number_retrieved bacteria1 bacteria2, ..."  # Optional: Provide a usage banner
  #maximal number of entries retrieved
  opts.on("-m", "--max MAX_NUMBER", Integer, "Set the maximal number of entries retrieved") do |number|
    options["max_number_retrieved"] = number
  end
  opts.on("-o", "--output_directory OUTPUT_DIRECTORY", String, "Set the output directory for saving data") do |dir|
    options["output_directory"] = dir
  end
  opts.on("-h", "--help", "Displays this help message") do
    print_usage_and_exit
  end  
end.parse!

#Check if the maximum number of entries was set to a valid number
if options.fetch("max_number_retrieved", 0) <= 0
  $stderr.puts("\nError: Invalid or missing value for maximum entries. Please specify a positive integer.")
  print_usage_and_exit
end
if options["output_directory"].nil? || !File.directory?(options["output_directory"]) || !File.writable?(options["output_directory"])
  $stderr.puts "\nError: The specified output directory '#{path}' does not exist or is not writable."
  print_usage_and_exit
end

#Set max number of retrieval.
max_retrieval_permitted = options["max_number_retrieved"]
max_retrieval_tried = max_retrieval_permitted + 100

#Set bacterial for data retrieval.
bacterias = []
($stderr.puts("Set bacteria as 'Genus_species', e.g., 'Escherichia_coli'."); print_usage_and_exit) if ARGV.empty?

ARGV.each do |arg|
  bacteria = arg.split('_')
  # Check if the bacteria name follows the required 'Genus_species' format
  if bacteria.length != 2 || bacteria.any?(&:empty?)
    $stderr.puts("Each bacteria must be specified in 'Genus_species' format, e.g., 'Escherichia_coli'.")
    print_usage_and_exit
  end
  bacterias << bacteria #bacterias = [["Bacillus","subtilis"], ["Proteus","vulgaris"], ["Shigella","flexneri"]]
end

remove_ids = [] #List entries that are not investigated.
wrote_ids = []
id_header_hash = {}

bacterias.each do |(genus,species)|
  query = "#{genus}[title] AND #{species}[title] AND complete[title] AND 16S"
  accession_numbers = nil
  begin
    params = {
      'db' => 'nuccore',
      'term' => query,
      'retmax' => max_retrieval_tried
    }
    xml_data = fetch_from_ncbi(SEARCH_BASE_URL, params)
    doc = Nokogiri::XML(xml_data)
    accession_numbers = doc.xpath('//IdList/Id').map(&:text)
    if accession_numbers.empty?
      $stderr.puts "No accessions found."
    else
      $stderr.puts "\n#{genus} #{species}: #{accession_numbers.size} accessions were temporarily retrieved."
    end
  rescue StandardError => e
    $stderr.puts "Error: #{e.message}"
  end

  entries_written = 0
  accession_numbers.each do |accession_number|
    next if remove_ids.include?(accession_number)
    # API request
    begin
      params = {
        'db' => 'nuccore',
        'id' => accession_number,
        'rettype' => 'ft',  #feature table
        'retmode' => 'text',  #text format.
      }
      feature_lines = fetch_from_ncbi(FETCH_BASE_URL, params).split("\n")
      rRNA_positions = find_16S_rRNA_positions(feature_lines)
    
      raise "No 16S ribosomal RNA sequences found in accession #{accession_number}." if rRNA_positions.empty?
    
      rRNA_positions.each do |(rRNA_16S_start,rRNA_16S_end)|
        strand = rRNA_16S_start <= rRNA_16S_end ? "sense" : "anti-sense"
        sequence_data = fetch_rRNA_sequence(FETCH_BASE_URL, accession_number, rRNA_16S_start, rRNA_16S_end)
        header, *sequence_lines = sequence_data.split("\n")
        sequence = sequence_lines.join("")
        File.open(File.join("#{options["output_directory"]}","#{genus}_#{species}_from_genbank_#{entries_written}.fasta"),'a') do |file|
          file.puts header
          file.puts strand == "sense" ? sequence : rc(sequence)
          $stderr.puts "A 16S ribosomal RNA sequence from genbank #{accession_number} is written."
          wrote_ids << accession_number
          id_header_hash[accession_number] = header
          entries_written += 1
        end
        File.rename(File.join("#{options["output_directory"]}","#{genus}_#{species}_from_genbank_#{entries_written-1}.fasta"), File.join("#{options["output_directory"]}","#{genus}_#{species}_from_genbank_#{entries_written}.fasta"))
        if entries_written >= max_retrieval_permitted
          break
        end
      end
      if entries_written >= max_retrieval_permitted
        break
      end
    rescue StandardError => e
      $stderr.puts "Error: #{e.message}"
    end
  end
  $stderr.puts "#{entries_written} entries were written for #{genus} #{species}."
end      
