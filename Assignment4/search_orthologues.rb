=begin
Assignment 4. Searching for Orthologues
Author: Gema Castillo García
Date: 21/12/2021

FILES I HAVE USED:
      TAIR.fa
      SPOM.fa

FILES I HAVE CREATED:
  - The main script that makes BLAST files and the report with orthologues is:
      search_orthologues.rb --> creates the following three files

  - The BLAST files are:
      best_tair.txt
      best_spom.txt

  - The report file with pairs of possible orthologues between "A. thaliana" and "S. pombe" is:
      orthologues_list.txt


### Final program to execute: $ ruby search_orthologues.rb


################################ Homologues and Orthologues ################################
Two genes are homologous when their DNA sequence derives from a common origin, and may or may not have the same function. The task of finding homologs (hits) to a sequence of interest (query) in a database containing many other sequences (targets) can be conceptualized as getting the best possible alignment of the query against all the targets, scoring each of these alignments, and choosing those whose scores surpass a given threshold such as E-value (the lower it is, the more significant is the alignment). The Basic Local Alignment Search Tool (BLAST) is probably the most common heuristic algorithm used to find homologs.

As rates of identity between orthologs correlate directly with time to the last common ancestor (https://www.sciencedirect.com/science/article/abs/pii/S1046202309001327?via%3Dihub), two orthologs may be very different. Consequently, I am not considering the identity percent as a BLAST filter in this assignment (I only got 1 pair of possible orthologs trying to filter by >=25%). Thus, to do this assignment, I have run BLAST using the following parameters: E-value <= 1e-6 and coverage >= 50%, as used in https://academic.oup.com/bioinformatics/article/24/3/319/252715.


########################################## BONUS ##########################################
***Orthologs*** are homologous sequences that diverged from speciation events, so they are expected to conserve the same function and are useful for tracing species evolution. The first step to find orthologs is to search for **Best Reciprocal Hits (BRH)** using BLAST (i.e., comparison of two genomes in a pairwise manner using amino acid sequences of each CDS). However, this method may miss true orthologs if it returns a paralog as the best hit, so it would be better to consider more than one hit between high-scoring hits. In addition, to complement this analysis, it can be used the **open reading frame (ORF)-independent method**, which uses nucleotide sequences of CDSs instead of amino acid sequences to find orthologous DNA regions.

After comparing the possible pairs of orthologues obtained using the BRH method with those of the ORF-independent method, the next step would be to perform a **functional analysis** looking at the GO ontolgies and the KEGG pathways (using Web APIs as we did in the Assigment 2) to determine whether both genes perform the same function and/or whether they are involved in the same pathway. In this respect, a more detailed **analysis of the domains** that they share could also be performed, and even a **phylogenetic tree** could be constructed to see if these sequences tend to cluster together (i.e., they have the same origin).


=end



#I have created a new folder "Databases" to store databases and the multifasta files of "Arabidopsis thaliana" and "Schizosaccharomyces pombe"

require 'bio'
require 'stringio'


#indexing FASTA files using makeblastdb to create BLAST databases for "A. thaliana" and "S. pombe" sequences
system("makeblastdb -in './Databases/TAIR.fa' -dbtype 'nucl' -out ./Databases/TAIR")
system("makeblastdb -in './Databases/SPOM.fa' -dbtype 'prot' -out ./Databases/SPOM")


def read_file(file) #defining a function to read and store the content of multifasta files
  abort("Please, check the name of #{file}, I can't find it!\n") unless (File.exists?(file)) #aborting if the file does not exist
  puts "Reading the file #{file}\n"
  @fasta = Bio::FlatFile.auto(file) #storing the content of .fa files as FASTA sequences
  puts "DONE!\n\n"
  return @fasta
end


tair_queries = {} #creating a dictionary to store each sequence (key=ID, value=sequence) of "A. thaliana"
spom_queries = {} #creating a dictionary to store each sequence (key=ID, value=sequence) of "S. pombe"

#reading "A. thaliana" multifasta file to store each sequence in a hash entry
read_file("./Databases/TAIR.fa").each do |query|
  tair_queries[query.entry_id] = query.seq
end

#reading "S. pombe" multifasta file to store each sequence in a hash entry
read_file("./Databases/SPOM.fa").each do |query|
  spom_queries[query.entry_id] = query.seq
end


#If a nucleotide sequence is translated before the search, it is more likely to find better and more accurate hits than just a blastn search.
#One of the reasons for this is that protein sequences are evolutionarily more conserved than nucleotide sequences. 

#to run the BLAST analysis:
#   1. create a BLAST factory
    factory_tair = Bio::Blast.local('tblastn', './Databases/TAIR') #the database is translated into protein
    factory_spom = Bio::Blast.local('blastx', './Databases/SPOM') #the queries are translated into proteins

#   2. run the actual BLAST by querying the factory
    tair_report = {} #creating a dictionary to store BLAST report of each "S. pombe" sequence (key=queryID, value=BLAST report) of "A. thaliana"
    spom_report = {} #creating a dictionary to store BLAST report of each "A. thaliana" sequence (key=queryID, value=BLAST report) of "S. pombe"
    

    #BLAST of "S. pombe" sequences (queries) vs "A. thaliana" database
    spom_queries.each_pair do |key_spom, value_spom|
      tair_report[key_spom] = factory_tair.query(">#{key_spom}\n#{value_spom}") #storing the BLAST report of each query (key) in a hash entry
    end


#defining a function to retrieve the best hits of each BLAST report to store them in an array
def retrieve_best_hits(blast_report,best_hits_array,output_file)
  File.open(output_file, "w") { |f| f.write("Query\tHit\tE-value\tCoverage (%)\tIdentity (%)\n") } #creating a new file to store E-value, coverage and identity of best hits
  blast_report.each do |qual| #for each entry of the hash with BLAST entries
    unless qual[1].iterations[0].message == "No hits found" #unless there is no hits
      #puts hits[1].hits[0].hsps
      qual[1].hits[0].hsps.each do |hsp| #extracting hits information
        if hsp.evalue <= 1e-6 #filtering hits by E-value because the smaller it is, the better is the hit
          ############ definitions of coverage and identity percent in  https://usuhs.libguides.com/c.php?g=468091&p=3260303 ############
          segment = hsp.query_to.to_f + 1 - hsp.query_from.to_f #length of the aligned segment of the query
          
          #percent of the query length included in the aligned segment of the query
          #(coverage = alignment segment/query length)
          coverage = (100*(segment/qual[1].hits[0].query_len.to_f)).to_i #calculating coverage as integer in percent

          #how similar the query is to the aligned sequences
          #(number of identities/alignment segment)
          identity = (100*(hsp.identity.to_f/segment)).to_i #calculating identity as integer in percent

          if coverage >= 50 #and identity >= 25 ##I am not filtering by identity (see reasons above)
            ############ create a report with the best hits ############
            definition = qual[1].hits[0].definition.split("|") #retrieving only the first element of hit's definition
            #puts "QUERY: #{qual[1].hits[0].query_def} | HIT: #{definition[0]} => #{qual[1].hits[0].hit_id} evalue: #{hsp.evalue}\tcoverage: #{coverage}%\tidentity: #{identity}%\n"   
            best_hits_array << qual[1].hits[0].definition[0..10] #adding the filtered best hit of each query to 'best_hits_array'
            File.open(output_file, "a") { |f| f.write("#{qual[1].hits[0].query_def}\t#{definition[0]}\t#{hsp.evalue}\t#{coverage}\t#{identity}\n") } #
            break #to avoid selecting other alignments of the same hit with a higher E-value (the first alignment always has the lower E-value)
          end
        end
      end
    end
  end
  best_hits_array = best_hits_array.uniq #removing possible duplicated genes from the list with the best hits
  return best_hits_array
end


@best_tair = [] #creating an array to return the best hits of the BLAST of "S. pombe" sequences vs "A. thaliana" database
best_tair = retrieve_best_hits(tair_report,@best_tair,"best_tair.txt") #storing the best hits of the BLAST of "S. pombe" sequences vs "A. thaliana" database in an array   


#BLAST of "A. thaliana" best hits vs "S. pombe" database
tair_queries.each_pair do |key_tair, value_tair|
  if best_tair.include? key_tair #only selecting the "A. thaliana" genes included in the list with the best hits to optimise the process
    spom_report[key_tair] = factory_spom.query(">#{key_tair}\n#{value_tair}") #storing the BLAST report of each query (key) in a hash entry
  end
end


@best_spom = [] #creating an array to return the best hits of the BLAST of "A. thaliana" sequences vs "S. pombe" database
best_spom = retrieve_best_hits(spom_report,@best_spom,"best_spom.txt") #storing the best hits of the BLAST of "A. thaliana" sequences vs "S. pombe" database in an array   


#I converted this box to 'Raw NBConvert' by accident after runing everything and this is why it seems that it has not been run, but it has
#(runing everything again so it appears run would take a lot of time!)

def find_orthologues(tair_homologues_file,spom_homologues_file) #defining a function to create a report file with pairs of possible orthologues between "A. thaliana" and "S. pombe"  
  File.open("orthologues_list.txt", "w") { |f| f.write("Pairs of possible orthologues between A. thaliana and S. pombe:\n") } #creating a new file to store pairs of possible orthologues
  
  tair_homologues = [] #creating an array to store 
  spom_homologues = []
  orthologues_list = []
  
  file = File.read(tair_homologues_file).split("\n") #reading the first file of homologues genes
  file.each do |pair|
    field = pair.split("\t")
    tair_homologues << [field[0],field[1]] #storing in an array the ID of each query and its best hit
  end
  
  file = File.read(spom_homologues_file).split("\n") #reading the second file of homologues genes
  file.each do |pair|
    field = pair.split("\t")
    spom_homologues << [field[0],field[1]] #storing in an array the ID of each query and its best hit
  end
  
  
  #comparing queries and hits of both homologues arrays to find reciprocal best hits
  tair_homologues.each do |tair_hom|
    spom_homologues.each do |spom_hom|
      if tair_hom[0] == spom_hom[1] and tair_hom[1][0..10] == spom_hom[0]
        orthologues_list << spom_hom
      end
    end
  end
  
  
  orthologues_list = orthologues_list.uniq #removing possible duplicated pairs of orthologues (just in case, but it should not be any)
  #puts orthologues_list
  orthologues_list.each do |line|
    File.open("orthologues_list.txt", "a") { |f| f.write("\t#{line[0]} & #{line[1]}\n") } #saving each pair of possible orthologues in a file
  end
  puts "Orthologues report is ready! Look at 'orthologues_list.txt' file"
end


find_orthologues("./best_tair.txt","./best_spom.txt") #creating a report file with pairs of possible orthologues between "A. thaliana" and "S. pombe"


