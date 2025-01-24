#! /usr/bin/ruby
require 'fileutils'
require 'tempfile'
require 'set'
require 'optparse'
require 'zlib'

THRESHOLD_FOR_SPECIFIC_SPECIES_H = {
  "Acinetobacter_baumannii_Region_I" => 559, 
  "Aerococcus_viridans_Region_I" => 569, 
  "Bacillus_anthracis_Region_I" => 687, 
  "Bordetella_pertussis_Region_I" => 469, 
  "Chlamydia_trachomatis_Region_I" => 687, 
  "Chlamydophila_pneumoniae_Region_I" => 792, 
  "Chlamydophila_psittaci_Region_I" => 803, 
  "Clostridium_botulinum_Region_I" => 369, 
  "Corynebacterium_jeikeium_Region_I" => 602, 
  "Cutibacterium_acnes_Region_I" => 521, 
  "Enterococcus_faecalis_Region_I" => 738, 
  "Enterococcus_faecium_Region_I" => 738, 
  "Escherichia_coli_Region_I" => 781, 
  "Francisella_tularensis_Region_I" => 423, 
  "Haemophilus_influenzae_Region_I" => 694, 
  "Haemophilus_parainfluenzae_Region_I" => 700, 
  "Klebsiella_pneumoniae_Region_I" => 834, 
  "Legionella_pneumophila_Region_I" => 454, 
  "Moraxella_catarrhalis_Region_I" => 563, 
  "Mycobacterium_avium_Region_I" => 840, 
  "Mycobacterium_gordonae_Region_I" => 788, 
  "Mycobacterium_intracellulare_Region_I" => 834, 
  "Mycobacterium_kansasii_Region_I" => 812, 
  "Mycobacterium_tuberculosis_Region_I" => 803, 
  "Mycolicibacterium_fortuitum_Region_I" => 676, 
  "Mycoplasma_pneumoniae_Region_I" => 220, 
  "Nocardia_brasiliensis_Region_I" => 676, 
  "Proteus_mirabilis_Region_I" => 677, 
  "Pseudomonas_aeruginosa_Region_I" => 559, 
  "Rickettsia_rickettsii_Region_I" => 325, 
  "Staphylococcus_aureus_Region_I" => 857, 
  "Staphylococcus_epidermidis_Region_I" => 857, 
  "Stenotrophomonas_maltophilia_Region_I" => 504, 
  "Streptococcus_mitis_Region_I" => 903, 
  "Streptococcus_oralis_Region_I" => 903, 
  "Streptococcus_pneumoniae_Region_I" => 903, 
  "Streptococcus_pyogenes_Region_I" => 674, 
  "Streptococcus_salivarius_Region_I" => 790, 
  "Streptococcus_sanguinis_Region_I" => 862, 
  "Trueperella_pyogenes_Region_I" => 489, 
  "Yersinia_pestis_Region_I" => 720, 
  "Human_HSD11B2" => 0
}

options = {}

#Obtain reverse complement.
def rc(seq)
  seq.upcase.reverse.tr('ACGTNWRMKYSHBVD', 'TGCANSYKMRWDVBH')
end

#Extract a single match from a blastn output.
def extract_substring_by_pattern(input_str, regex, word, stop_patterns)
  match_data = input_str.match(regex) #str matches regex?
  return nil unless match_data #If not, return.
  post_regex_str = input_str[match_data.end(0)..] #str from match to the end.
  start_index = post_regex_str.index(word) #Search match with word.
  return nil unless start_index #If not, return.
  post_word_str = post_regex_str[start_index..] #str from match to the end.
  extracted_lines = []
  post_word_str.each_line do |line|
    if stop_patterns.any? { |pattern| line.include?(pattern) } #If any of pattern appears,
      return extracted_lines.join("\n") #Return lines before the pattern.
    else
      extracted_lines << line
    end
  end
  nil #Stop pattern never appears.
end

def exit_with_error(message)
  $stderr.puts message
  $stderr.puts "Terminating program."
  exit
end

def judgement(sequence, bacteria_data_ar_ar) #seq, [["Chlamydophila_psittaci_part_I", 944],[Chlamydophila_…, 823],...]  
  bacteria_data_ar_ar.each_with_index do |(bacteria, score), index|
    (index == 0 && score.to_f >= THRESHOLD_FOR_SPECIFIC_SPECIES_H[bacteria]) ? (return bacteria) : break
  end
  streptococcus_bacteria_set = bacteria_data_ar_ar.take(5).map{|ar| ar[0]}.to_set
  _, streptococcus_bacteria_score = bacteria_data_ar_ar[4]
  if streptococcus_bacteria_set == ["Streptococcus_pneumoniae_part_I", "Streptococcus_mitis_part_I", "Streptococcus_oralis_part_I", "Streptococcus_sanguinis_part_I", "Streptococcus_salivarius_part_I"].to_set
    if streptococcus_bacteria_score.to_i >= 674
      return "Streptococcus_pneumoniae_part_I" if sequence.include?("TGCACTTGCA")
      return "alpha-Streptococcuses"
    end
  end
  five_bacteria_set = bacteria_data_ar_ar.take(5).map{|ar| ar[0]}.to_set
  _, fifth_bacteria_score = bacteria_data_ar_ar[4]
  if five_bacteria_set == ["Mycobacterium_gordonae_part_I", "Mycobacterium_kansasii_part_I", "Mycobacterium_avium_part_I", "Mycobacterium_intracellulare_part_I", "Mycobacterium_tuberculosis_part_I"].to_set
    if fifth_bacteria_score.to_i >= 646
      return "Mycobacterium_avium_part_I" if sequence.include?("CCTCTTCGGA") && sequence.include?("ATAGGACCTC")
      return "Mycobacterium_gordonae_part_I" if sequence.include?("TGTCCTGTGGT") && sequence.include?("GGTGATGG") #M.intracellulare also has the latter.
      return "Mycobacterium_intracellulare_part_I" if sequence.include?("CGGGGGTACT") && sequence.include?("ATAGGACCTT")
      return "Mycobacterium_kansasii_part_I" if sequence.include?("GACACTCGAG")
      return "Mycobacterium_tuberculosis_part_I" if sequence.include?("AAAGCGCTTT") && sequence.include?("GCTTTAGCGG")
    end
  end
  return "Undetermined"
end

class Blast_search
  attr_reader :database_h
  
  def initialize(options)
    @options = options
    blast_dir = @options["blast_directory"]
    makeblastdb_path = blast_dir.nil? ? "makeblastdb" : "#{blast_dir}/bin/makeblastdb"
    @blastn_path = blast_dir.nil? ? "blastn" : "#{blast_dir}/bin/blastn"
    #Copy database.
    database_name_without_extention = File.basename(@options["database"], ".*")
    @working_database_path = File.join(@options["output_directory"], "#{Process.pid}_#{File.basename(@options["database"])}")
    @working_database_path_without_extension = File.join(File.dirname(@working_database_path), File.basename(@working_database_path, ".*"))
    working_database_name_without_extention = File.basename(@working_database_path_without_extension) #The path of working .fasta file minus .fasta extention. This is used by the makeblastdb program.
    FileUtils.cp(@options["database"], @working_database_path)
    #Make database files for blast.
    begin
      `#{makeblastdb_path} -in #{@working_database_path} -dbtype nucl -parse_seqids -out #{File.join(@options["output_directory"], working_database_name_without_extention)}`
    rescue
      exit_with_error("*****Blast is not installed or the path of blast is not set.*****")
    end
      #Read datafase into @database_h
    @database_h = read_fasta_h(content_ar(options["database"]))#["S.pneumoniae_NC_003098_part_I"] = "AGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCCTAATAC..."
  end
  
  def execute(id_with_underscore, sequence)
    temp_file = Tempfile.new('my_temp_file')
    begin
      temp_file.open
      temp_file.puts ">#{id_with_underscore}"
      temp_file.puts sequence
      temp_file.close
      blastn_result = `#{@blastn_path} -db #{@working_database_path_without_extension} -query #{temp_file.path}`
    rescue
      exit_with_error("*****Blast is not installed or the path of blast is not set.*****")
    ensure
      temp_file.unlink
    end
    blastn_result.include?("No hits found") ? nil : blastn_result
  end
  
  def content_ar(filename) #[line1, line2, ...]
    if File.extname(filename) == '.gz'
      tempfile = Tempfile.new('decompressed_file')
      begin
        Zlib::GzipReader.open(filename) do |gz|
          tempfile.write(gz.read)
        end
        tempfile.rewind
        File.readlines(tempfile.path).map(&:chomp)
      ensure
        tempfile.close
        tempfile.unlink
      end
    else
      File.readlines(filename, chomp: true)
    end
  end
  
  #Take the common part of the IDs of #R1 and #R2 as long as possible.
  def common_prefix(str1, str2)
    min_length = [str1.length, str2.length].min
    result = ""

    (0...min_length).each do |i|
      if str1[i] == str2[i]
        result += str1[i]
      else
        break
      end
    end

    result
  end

  #Read fasta file. Return Hash h[header] = seq
  def read_fasta_h(content_ar)
    sequences_h = {}
    current_sequence = ''
    content_ar.each do |line|
      line.chomp!
      if line.start_with?('>')
        #Start of a new sequence.
        current_sequence = line[1..-1] #Remove '>' from the title.
        sequences_h[current_sequence] = [] unless current_sequence.empty?
      elsif !line.empty?
        #The sequence continues.
        sequences_h[current_sequence] << line
      end
    end
    sequences_h.each do |k,v|
      sequences_h[k] = v.join('')
    end
    sequences_h
  end
  
  def read_fastq_h(content_ar)
    content_ar.each_slice(4).each_with_object({}) do |(id, seq, _, _), sequences_h|
      sequences_h[id[1..-1]] = seq
    end
  end
    
  def read_sequences_h(pass_ar) #return hash[id] = seq1
    raise "Please provide at least one file" if pass_ar.empty?
    if pass_ar.size == 1
      content_ar = content_ar(pass_ar[0])
      first_char = content_ar[0][0]
      case first_char
      when '>' #fasta
        [read_fasta_h(content_ar), 'fasta']
      when '@' #fastq
        [read_fastq_h(content_ar), 'fastq']
      else
        raise "Unknown file format."
      end
    elsif pass_ar.size == 2 #paired read
      print "\nCombining #{File.basename(pass_ar[0])} and #{File.basename(pass_ar[1])} "
      r1_content = content_ar(pass_ar[0])
      r2_content = content_ar(pass_ar[1])
      combined_seq = combine_r1_and_r2(r1_content, r2_content)
      [combined_seq, 'fastq']
    else
      raise "Unknown file format."
    end
  end
  
  #Read R1 and R2 .fastq files or .fastq.gz files. Paired reads are combined into a single read.
  def combine_r1_and_r2(fastq_ar1, fastq_ar2)
    combined_id_seq_h = {}
    reads_total, reads_discarded, reads_utilized = 0, 0, 0
    (0...(fastq_ar1.size / 4)).each do |k|
      reads_total += 1
      print "*" if reads_total % 10 == 1
      read1_seq = fastq_ar1[k * 4 + 1]
      read2_seq = rc(fastq_ar2[k * 4 + 1]) #R2 file is read after converting to its reverse complement.
      if (read1_seq.size < 50 || read1_seq.count("N") > 2) || (read2_seq.size < 50 || read2_seq.count("N") > 2) #Remove low quality data.
        reads_discarded += 1
        next
      end
      #fastq_ar1[k * 4] =~ /^(\S+)/ #Obtain ID without any space in it.
      id = common_prefix(fastq_ar1[k * 4], fastq_ar2[k * 4]).strip[1..-1] #Take as long a common part of #ID as possible with r1 and r2. Exclude the initial @.
      id_under_bar = id.gsub(/[: ,\/\.\(\)]+/, "_") #Replace ":" in the id string with "_".
    
      st_r1 = execute(id, read1_seq)
      if st_r1.nil? #Remove a sequence without any homology.
        reads_discarded += 1
        next
      end
      st_r2 = execute(id, read2_seq)
      if st_r2.nil? #Remove a sequence without any homology.
        reads_discarded += 1
        next
      end    
      st_r1 =~ /Sequences producing significant alignments:.+\n\n([^>]+)\n\n>/
      hit_list1_ar_ar = $1.nil? ? [] : $1.scan(/(\S+)\s+([\d\.]+)/) #[["S.pneumoniae_NC_003098_Part_I", 520],...]
      st_r2 =~ /Sequences producing significant alignments:.+\n\n([^>]+)\n\n>/
      hit_list2_ar_ar = $1.nil? ? [] : $1.scan(/(\S+)\s+([\d\.]+)/) #[["S.pneumoniae_NC_003098_Part_I", 411],...]
      hit_list_h = {}
      hit_list1_ar_ar.each {|(amplicon, score)|
        hit_list_h[amplicon] = score.to_f #hit_list_h["S.pneumoniae_NC_003098_Part_I"] = 520
      }
      #Select entries that exist in both hit_list1_ar_ar and hit_list2_ar_ar.
      hit_list2_ar_ar.reject! { |(amplicon, score)| hit_list_h[amplicon].nil? }
      hit_list2_ar_ar.each {|(amplicon, score)|
        hit_list_h[amplicon] += score.to_f
      }
      sorted_hit_list_ar_ar = hit_list_h.to_a.sort{|a,b| b[1] <=> a[1]}  #sort by score in descending order. #[["S.pneumoniae_NC_003098_Part_I", 931],...]
      if sorted_hit_list_ar_ar.size >= 2 && sorted_hit_list_ar_ar[0][1] == sorted_hit_list_ar_ar[1][1]
        reads_discarded += 1
        next #If 2 items have the same score, discard the entry.
      end
      top_hit_entry = sorted_hit_list_ar_ar[0][0] #Top amplicon like "S.pneumoniae_NC_003098_Part_I"
          
      t_1 = extract_substring_by_pattern(st_r1, /^>#{top_hit_entry}.+\nLength/, "Query", ["Length", ">", "Score", "Lambda"])
      if t_1.nil?
        reads_discarded += 1
        next
      end
      t_2 = extract_substring_by_pattern(st_r2, /^>#{top_hit_entry}.+\nLength/, "Query", ["Length", ">", "Score", "Lambda"])
      if t_2.nil?
        reads_discarded += 1
        next
      end
      query_ar_1 = t_1.scan(/Query\s+(\d+)\s+([AGCT\-]+)\s+(\d+)/)
      sbjct_ar_1 = t_1.scan(/Sbjct\s+(\d+)\s+([AGCT\-]+)\s+(\d+)/)
      query_ar_2 = t_2.scan(/Query\s+(\d+)\s+([AGCT\-]+)\s+(\d+)/)
      sbjct_ar_2 = t_2.scan(/Sbjct\s+(\d+)\s+([AGCT\-]+)\s+(\d+)/)
      query_seq_1 = query_ar_1.map { |ar| ar[1].gsub('-', '') }.join
      sbjct_start_1 = sbjct_ar_1[0][0].to_i-1
      query_seq_2 = query_ar_2.map { |ar| ar[1].gsub('-', '') }.join
      sbjct_start_2 = sbjct_ar_2[0][0].to_i-1
    
      database_sequence = @database_h[top_hit_entry] #Obtain sequence of the top hit entry.
      database_sequence_size = database_sequence.size
    
      r1_aligned_sequence = ' ' * database_sequence_size
      r2_aligned_sequence = ' ' * database_sequence_size
      r1_aligned_sequence[sbjct_start_1, query_seq_1.size] = query_seq_1.upcase
      r2_aligned_sequence[sbjct_start_2, query_seq_2.size] = query_seq_2.upcase
      aligned_sequence = ' ' * database_sequence.size
      (0...(database_sequence.size)).each do |i|
        if r1_aligned_sequence[i] == r2_aligned_sequence[i] && r1_aligned_sequence[i] != ' '
          aligned_sequence[i] = r1_aligned_sequence[i]
        elsif database_sequence[i] == r1_aligned_sequence[i]
          aligned_sequence[i] = r1_aligned_sequence[i]
        elsif database_sequence[i] == r2_aligned_sequence[i]
          aligned_sequence[i] = r2_aligned_sequence[i]
        elsif r1_aligned_sequence[i] != ' '
          aligned_sequence[i] = r1_aligned_sequence[i]
        elsif r2_aligned_sequence[i] != ' '
          aligned_sequence[i] = r2_aligned_sequence[i]
        else
          aligned_sequence[i] = database_sequence[i]
        end
      end
      stripped_sequence = aligned_sequence.strip
      if stripped_sequence.index(" ") #Exclude any string having " ".
        reads_discarded += 1
        next
      end
      combined_id_seq_h[id] = stripped_sequence
      reads_utilized += 1
    end
    $stderr.puts "\nRead_pair_examined =\t#{reads_total}, \tRead_pair_utilized =\t#{reads_utilized}, \tRead_pair_discarded =\t#{reads_discarded}."
    combined_id_seq_h #combined_id_seq_h[id] = combined_seq
  end

  def exit_with_error(message)
    $stderr.puts message
    $stderr.puts "Terminating program."
    exit
  end
  
  def cleanup()
    File.delete("#{@working_database_path}")
    Dir.glob("#{@working_database_path_without_extension}.*").each do |file|
      File.delete(file)
    end
  end
end

def calculate_string_frequencies(ar_ar, index)
  frequencies_h = ar_ar.each_with_object(Hash.new(0)) { |ar, counts| counts[ar[index]] += 1 }
  max_frequency_ar = frequencies_h.max_by { |k, v| v }
  [max_frequency_ar[0], max_frequency_ar[1], ar_ar.count] #Return the mode string and its occurrence count.
end

def ret_score(rows_ar_ar, identified_species)
  rows_ar_ar.each do |(species, score, _)|
    if species == identified_species
      return score
    end
  end
  return ""
end

#Main program
#1. Directory information.
option_parser = OptionParser.new do |opts|
  # output directory option
  opts.on("-o DIR", "--output-directory DIR", "Set the output directory") do |dir|
    options["output_directory"] = dir
  end
  # database option
  opts.on("-d DATABASE", "--database DATABASE", "Specify the database") do |db|
    options["database"] = db
  end
  # blast option
  opts.on("-b BLAST", "--blast BLAST", "Specify the database") do |blast|
    options["blast_directory"] = blast
  end
  # help message
  opts.on("-h", "--help", "Prints this help") do
    puts opts
    exit
  end
end
option_parser.parse!

master_database_path = options["database"] #master_database_path
exit_with_error("\n***** Database file is not set or does not exist.*****") if master_database_path.nil? || !File.exist?(master_database_path)

$stderr.puts "\nDatabase is\t#{options["database"]}"
output_dir = options["output_directory"]
exit_with_error("\n***** Output directory is not set.*****") if output_dir.nil?
exit_with_error("\n***** Output directory is not found or is not writable.*****") if !Dir.exist?(output_dir) || !File.writable?(output_dir)
$stderr.puts "\nOutput directory is\t#{options["output_directory"]}"
output_dir = options["output_directory"]
exit_with_error("Output directory not found or not writable: #{output_dir}.") unless Dir.exist?(output_dir) || File.writable?(output_dir)
master_database_path = options["database"] #master_database_path
exit_with_error("Database file is not set or does not exist.\nThe command line should look like: >ruby top_hits_for_each_read.rb bacteria_database.fasta") if master_database_path.nil? || !File.exist?(master_database_path)
output_files = []

files_ar = ARGV.sort_by do |filename|
  match_data = filename.match(/S(\d+)_/)
  if filename =~ /S(\d+).+(R1|R2)/
    [0, $1.to_i, $2] #If S+number is found, sort using that.
  else
    [1, filename] #If not, sort as a string.
  end
end

#2. Make a Blast_search instance.
blast_search = Blast_search.new(options)

#3.	Read database into database_h["S.pneumoniae_NC_003098_part_I"] = "AGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCCTAATAC..."
database_h = blast_search.database_h

#4-1. Analyze fasta or fastq files.
while files_ar.size != 0
  if files_ar.size >= 2 && files_ar[0].include?('_R1_') && files_ar[1].include?('_R2_') #paired read.
    pass_ar = [files_ar.shift, files_ar.shift]
  else
    pass_ar = [files_ar.shift]
  end
  
  file_basename = File.basename(pass_ar.first).split('.')[0]
  file_stemname = file_basename.include?('_R1_') ? file_basename.match(/(.*?)_R1_/)[1] : file_basename

#4-2. Read fasta file.
  begin
    id_seq_ar_ar, file_format = blast_search.read_sequences_h(pass_ar) #id_seq_ar_ar = [[id, seq1], [id2, seq2], ....]
  rescue RuntimeError => e
    puts "An unexpected error has occurred. #{e.message}．"
    exit 1
  end

#4-3. Blast search.
  identifiedSpecies_count_h = Hash.new(0)
  id_identifiedSpeciesScoreAr_h = {}
  
  $stderr.print "Processing #{file_stemname} "
  output_file = File.join(output_dir, "#{file_stemname}_#{File.basename(master_database_path,".*")}_result.txt")
  output_files << output_file
  File.open(output_file, "w") do |file|
    data_ar_ar = []
    id_seq_ar_ar.each_with_index do |(id, seq), index| #id_seq_ar = [ID, seq]
      $stderr.print "#{index + 1} " if (index + 1) % 10 == 0 #"#{id}\n" #if index % 10 == 0 #Indicating progress.
      id_with_underscore = id.gsub(/[: ,\/\.\(\)]+/, "_") #Replace ": ,()/" or "__" in the id string with "_".
      blastn_result = blast_search.execute(id_with_underscore, seq)
      if blastn_result.nil? #a sequence without any homology.
        identifiedSpecies_count_h["Undetermined"] += 1
        id_identifiedSpeciesScoreAr_h[id] = ["Undetermined",""]
        next
      end
      blastn_result =~ /Sequences producing significant alignments:.+\n\n/
      intermed_st = $'
      rows_ar_ar = intermed_st.split("\n\n").first.split("\n").map {|line| line.split(/\s+/)} # [["Chlamydophila_psittaci_part_I", "944", "0.0"],[Chlamydophila_…, "823", "0.0"],...]
      identified_species = judgement(seq, rows_ar_ar) #"Chlamydophila_psittaci_part_I"
      identifiedSpecies_count_h[identified_species] += 1
      id_identifiedSpeciesScoreAr_h[id] = [identified_species, identified_species == "Undetermined" ? [""] + rows_ar_ar.map{|ar| ar[0..1]}.flatten : ret_score(rows_ar_ar, identified_species)]
    end
    identifiedSpeciesCount_ar_ar = identifiedSpecies_count_h.sort_by {|_,v| -v} #[["Acinetobacter_baumannii_part_I", 995], ["Undetermined", 5]]
    determined_ar_ar = identifiedSpeciesCount_ar_ar.reject{|(identifiedSpecies,_)| identifiedSpecies == "Undetermined"}
    determined_ar_ar.each{|identifiedSpecies_count_ar| file.puts identifiedSpecies_count_ar.join("\t")}
    undetermined_ar_ar = identifiedSpeciesCount_ar_ar.select{|(identifiedSpecies,_)| identifiedSpecies == "Undetermined"}
    undetermined_ar_ar.each{|undetermined_ar| file.puts undetermined_ar.join("\t")}

    file.puts "\nID\tDiagnosis\tbit_score"
    id_identifiedSpeciesScoreAr_h.each do |id_with_underscore, identifiedSpecies_score_ar|
      file.puts "#{id_with_underscore}\t#{identifiedSpecies_score_ar.join("\t")}"
    end
  end
  $stderr.print "\n"
end

#5. Summarize
summary_file = File.join(output_dir, "summary.txt")
File.open(summary_file, 'w') do |file|
  output_files.each do |filename|
    basename = File.basename(filename)
    lines = File.readlines(filename).map(&:strip).take_while { |line| !line.empty? }
    formatted_line = "#{basename}\t#{lines.map { |line| line.split("\t") }.flatten.join("\t")}"
    file.puts formatted_line
  end
end

#6. Clean up.
blast_search.cleanup
