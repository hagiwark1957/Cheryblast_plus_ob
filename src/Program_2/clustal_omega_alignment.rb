#!/usr/bin/env ruby
require 'optparse'

# Parses command-line options and returns a hash of options
def parse_options
  options = {}
  OptionParser.new do |opts|
    opts.banner = "Usage: clustalo_alignment_using_multiple_threads.rb [options] FASTA1 FASTA2 ..."
    opts.on("-o", "--output-directory DIRECTORY", "Output directory") do |dir|
      options[:output_directory] = dir
    end
  end.parse!
  options
end

# Validates the presence and correctness of required options
def validate_options(options)
  raise OptionParser::MissingArgument, "Output directory (-o) is not specified." if options[:output_directory].nil?
  raise OptionParser::MissingArgument, "Output directory (-o) is not a directory." if !File.directory?(options[:output_directory])
end

# Parses filename to extract bacteria name and number of entries
def parse_filename(file)
  File.basename(file).match(/(^[A-Za-z]+_[A-Za-z]+)_[^\d]+(\d+)/).captures
end

# Conducts clustalo alignment for each FASTA file
def perform_alignment(file, options, output_stream)
  bacteria_name, n_entries = parse_filename(file)
  options[:genus_species] = "#{bacteria_name}_#{n_entries}"
  output_file_path = File.join(options[:output_directory], "#{bacteria_name}_#{n_entries}_Clustalo.txt")
  execute_clustalo(file, output_file_path)
  consensus = generate_consensus(output_file_path)
  write_consensus(bacteria_name, consensus, output_stream) unless consensus.empty?
end

# Executes clustalo command line tool
def execute_clustalo(input_file, output_file_path)
  $stderr.puts "Processing #{File.basename(input_file)}..."
  system("clustalo -i #{input_file} -o #{output_file_path} --iter=4 --outfmt=Clustal --wrap=10000")
end

# Generates a consensus sequence from the alignment file
def generate_consensus(file_path)
  lines = IO.readlines(file_path)[3..-2] # Strip header and footer
  consensus = lines.each_with_object([]) do |line, strings|
    match = line.match(/\s+([AGCT\-]+)/) #Sequences have already pudded by '-' by clustal omega.
    strings << match[1] unless match.nil?
  end
  consensus.empty? ? nil : calculate_consensus(consensus)
end

# Calculates the consensus sequence from an array of strings
def calculate_consensus(strings)
  max_length = strings.max_by(&:length).length
  padded_strings = strings.map { |s| s.ljust(max_length) }
  consensus = padded_strings.map(&:chars).transpose.map do |chars|
    determine_consensus(chars)
  end.join
end

def calculate_frequency(chars)
  chars.each_with_object(Hash.new(0)) { |char, counts| counts[char] += 1 }
end

def determine_nucleotide(most_frequent, second_frequent)
  pairs = {
    'AG' => 'R', 'CT' => 'Y', 'GC' => 'S',
    'AT' => 'W', 'GT' => 'K', 'AC' => 'M'
  }
  sorted_pair = [most_frequent, second_frequent].sort.join
  pairs[sorted_pair] || 'N'
end

def determine_consensus(chars)
  frequency = calculate_frequency(chars)
  sorted_by_frequency = frequency.sort_by { |_, count| -count }
  most_frequent, highest_count = sorted_by_frequency.first
  second_frequent, second_count = sorted_by_frequency[1] || [nil, nil]
  return '-' if most_frequent == '-'
  total_valid = chars.count { |char| char != '-' }
  if total_valid * 0.5 < highest_count
    most_frequent
  elsif total_valid * 0.35 < highest_count && total_valid * 0.35 < second_count
    determine_nucleotide(most_frequent, second_frequent)
  else
    'N'
  end
end
  
# Writes the consensus sequence to the output stream
def write_consensus(bacteria_name, consensus, output_stream)
  output_stream.puts ">#{bacteria_name}\n#{consensus.gsub("-", "").gsub(/\AN+|N+\z/, '')}"
end

def clustalo_alignment_of_consensus_seqs(options)
  `clustalo -i #{options[:consensus_for_each_bacteria]} -o #{options[:consensus_for_bacterias]} --iter=4 --outfmt=Clustal --wrap=10000`
end

#Main routine
options = parse_options
validate_options(options)

threads = []
max_concurrent_processes = 16
#Output *_Clustalo.txt
ARGV.each do |file|
  #Call a method asynchronously.
  threads << Thread.new { perform_alignment(file, options, $stderr) }
  while threads.count >= max_concurrent_processes
    threads.reject!{ |thread| !thread.alive? }
    sleep(0.5)
  end
end
while !threads.empty?
  threads.reject!{ |thread| !thread.alive? }
  sleep(0.5)
end

#Output 16S_ribosomal_RNAs_for_all_bacteria.fasta
consensus_for_each_bacteria = "16S_rRNA_consensus_for_each_bacteria.fasta"
$stderr.puts "\nOutput #{consensus_for_each_bacteria}\n"
options[:consensus_for_each_bacteria] = File.join(options[:output_directory], consensus_for_each_bacteria)
output_stream = File.open(options[:consensus_for_each_bacteria], 'w')
ARGV.each do |file|
  perform_alignment(file, options, output_stream)
end

# Align the final consensus sequences
consensus_for_bacterias = "16S_rRNA_consensus_for_bacterias.txt"
$stderr.puts "\n#{consensus_for_bacterias}\n"
options[:consensus_for_bacterias] = File.join(options[:output_directory], consensus_for_bacterias)
clustalo_alignment_of_consensus_seqs(options)
