=begin
Assignment 3. GFF feature files and visualization
Author: Gema Castillo García
Date: 27/11/2021

FILES I HAVE USED:
      ArabidopsisSubNetwork_GeneList.txt

FILES I HAVE CREATED:
  - The main script that makes the GFF files and the report is:
      make_files.rb --> creates the following four files

  - The GFF files are:
      gene_coordinates.gff
      chr_coordinates.gff

  - The report file with the genes without the repeat is:
      genes_without_rep.txt

 - The file with the GFF track beside the AT2G46340 gene on the ENSEMBL website is:
      Arabidopsis_thaliana_219021079_19028603.png


### Final program to execute: $ ruby make_files.rb
=end



require 'bio'
require 'rest-client'  

#creating a function called "fetch" that we can re-use everywhere in our code ()
def fetch(url, headers = {accept: "*/*"}, user = "", pass="") #fetch definition was written by Mark Wilkinson
  response = RestClient::Request.execute({
    method: :get,
    url: url.to_s,
    user: user,
    password: pass,
    headers: headers})
  return response
  
  rescue RestClient::ExceptionWithResponse => e
    $stderr.puts e.inspect
    response = false
    return response  #now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue RestClient::Exception => e
    $stderr.puts e.inspect
    response = false
    return response  #now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue Exception => e
    $stderr.puts e.inspect
    response = false
    return response  #now we are returning 'False', and we will check that with an \"if\" statement in our main code
end


def read_file(file) #defining a function to read and store the content of "ArabidopsisSubNetwork_GeneList.txt"
  abort("Please, check the name of #{file}, I can't find it!\n") unless (File.exists?(file)) #aborting if the file does not exist
  puts "Reading the file...\n"
  @agi = File.read(file).split("\n") #reading the file
  return @agi
end


###EX 1. Using BioRuby, examine the sequences of the ~167 Arabidopsis genes from the last assignment by retrieving them from whatever database you wish.###  
def retrieve_seq(agi_list) #defining a function to obtain the sequences of the AGI codes and their positions for the target "CTTCTT"
  File.open("gene_coordinate.gff", "w") { |f| f.write("##gff-version 3\n") } #creating a GFF file for EX 4a
  File.open("chr_coordinate.gff", "w") { |f| f.write("##gff-version 3\n") } #creating a GFF file for EX 5

  fasta = {} #creating a hash for Bio::Sequence objects (key=AGI) - all the genes of the list -
  @repeats = {} #creating a hash for Bio::Sequence objects with Bio::Feature objects (key=AGI) - only the genes with repeats -
  if @agi.empty?
    puts "This file is empty!" #aborting if the file is empty
  else
    puts "Retrieving sequences from Ensembl with Dbfetch..."
    @agi.each do |geneid| #for each gene (geneid) of the file do...
      url = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&style=raw&id=#{geneid}" #retrieveing Ensembl information of 'geneid'
      res = fetch(url)
      if res  #res is either the response object, or False, so you can test it with 'if'
        body = res.body  #getting the "body" of the response
        if body.empty?
          puts "There is no entry for #{geneid} in Ensembl." #printing a friendly message if the gene has not an entry in Ensembl
        else
          embl = Bio::EMBL.new(res.body) #creating a Bio::EMBL object for 'geneid'
          #puts "\n\n\n" + embl.sv #showing chromosomic coordinates of 'geneid' for EX 5
          fasta[geneid] = Bio::Sequence.new(embl.to_biosequence) #creating a Bio::Sequence::NA object for 'geneid'
          #puts ">" + geneid + "\n" + fasta[geneid] #showing 'geneid' sequence in FASTA format
          
          
          #inserts will be targeted to the repeat "CTTCTT" on both the + and - strands of DNA
          rep_f = Bio::Sequence::NA.new("CTTCTT") #to search rep "CTTCTT" on + strand
          rep_f_re = /(?=(#{rep_f.to_s}))/ #this step allows to search for overlapping targets by tranforming 'rep_f' into a reg expression
          #as the EMBL's sequences are forward, I also have to search for the complementary sequence of that repeat, that is, "AAGAAG"
          rep_r = rep_f.reverse_complement #to search rep "CTTCTT" on - strand
          rep_r_re = /(?=(#{rep_r.to_s}))/ #this step allows to search for overlapping targets by tranforming 'rep_r' into a reg expression
          
          ###EX 2. Loop over every exon feature, and scan it for the "CTTCTT" sequence.###
          embl.features do |feature| #looking for the repetition "CTTCTT" in geneid's exons
            if feature.feature == "exon"
              if /:/.match(feature.position) #ignoring overlaped exons
                next
              elsif feature.position =~ /complement([^"]+)/ #searching for exons on - strand (complementary strands)
                #puts "\t\nEXON in reverse strand:" + $1[1..$1.length-2] #showing gene coordinates of the exon
                pos1 = $1[1..$1.length-2]
                pos1_1 = pos1.match(/(?<start>\d+)..(?<end>\d+)/)
                pos1_2 = (pos1_1[:start].to_i..pos1_1[:end].to_i) #getting coordinates of the exon in BIO-informatics style (first nucleotide is '1')
                exon = fasta[geneid][pos1_2.begin-1..pos1_2.end-1] #changing coordinates to informatics style (first nucleotide is '0') to retrieve the exon's sequence  
                #puts exon
                positions_f = exon.enum_for(:scan, rep_f_re).map { Regexp.last_match.begin(0) } #getting coordinates of the forward repetitions inside the exon
                puts "\texon positions with #{rep_f.to_s}: #{positions_f}" unless positions_f.empty? #showing the exon coordinates of the forward repeat
                positions_r = exon.enum_for(:scan, rep_r_re).map { Regexp.last_match.begin(0) } #getting coordinates of the reverse repetitions inside the exon
                #puts "\texon positions with #{rep_r.to_s}: #{positions_r}" unless positions_r.empty? #showing the exon coordinates of the reverse repeat
                
                
                
                if embl.sv =~ /chromosome:TAIR10:([^"]+):([^"]+):([^"]+):([^"]+)/ #searching for chromosomic coordinates for EX 5
                  chr = $1 #retrieving the chromosome
                  chr_pos = ($2.to_i..$3.to_i) #retrieving the chromosomic coordinates of the gene
                  
                  ###EX 3. Create a new Sequence Feature with the coordinates of every CTTCTT sequence and add it to the EnsEMBL Sequence object.### 
                  fasta[geneid].features = []  #creating an array to store the new feature for geneid's sequence (Bio::Sequence object)
                  positions_f.each do |f| ##creating a Bio::Feature object to store the information about the forward repeats
                    @feature = Bio::Feature.new('exon_rep',pos1_2,
                      [ Bio::Feature::Qualifier.new("cttctt", [pos1_2.select.with_index { |pos, idx| idx==f }[0], pos1_2.select.with_index { |pos, idx| idx==f+5 }[0]] ), 
                        Bio::Feature::Qualifier.new('strand', '+')
                        ])
                    fasta[geneid].features << @feature #adding the new Bio::Feature object to the Bio::Sequence object
                    #puts "\t\tgene position #{pos1_2.select.with_index { |pos, idx| idx==f }}: #{@feature}" #showing the gene coordinates of the Bio::Feature object

                    ###EX 4. Loop over each one of your CTTCTT features and create a GFF3-formatted file of these features###
                    File.open("gene_coordinate.gff", "a") { |f| f.write("#{seqid=geneid}\t#{source="."}\t#{type="insertion_site"}\t#{start=@feature.qualifiers[0].value[0]}\t#{final=@feature.qualifiers[0].value[1]}\t#{score="."}\t#{strand=@feature.qualifiers[1].value}\t#{phase="."}\t#{attributes="ID=#{geneid}"};repeat=#{type=@feature.qualifiers[0].qualifier}\n") }   
                    
                    ###EX 5. Re-execute your GFF file creation so that the CTTCTT regions are now in the full chromosome coordinates used by EnsEMBL###
                    File.open("chr_coordinate.gff", "a") { |f| f.write("#{seqid=chr}\t#{source="."}\t#{type="insertion_site"}\t#{start=chr_pos.select.with_index { |pos, idx| idx==@feature.qualifiers[0].value[0]-1 }}\t#{final=chr_pos.select.with_index { |pos, idx| idx==@feature.qualifiers[0].value[1]-1 }}\t#{score="."}\t#{strand=@feature.qualifiers[1].value}\t#{phase="."}\t#{attributes="ID=#{geneid}"};repeat=#{type=@feature.qualifiers[0].qualifier}\n") }  
                  end

                  positions_r.each do |r| ##creating a Bio::Feature object to store the information about the reverse repeats
                    @feature = Bio::Feature.new('exon_rep',pos1_2,
                      [ Bio::Feature::Qualifier.new("aagaag", [pos1_2.select.with_index { |pos, idx| idx==r }[0], pos1_2.select.with_index { |pos, idx| idx==r+5 }[0]] ),
                        Bio::Feature::Qualifier.new('strand', '-')
                        ])
                    fasta[geneid].features << @feature #adding the new Bio::Feature object to the Bio::Sequence object
                    #puts "\t\tgene position #{pos1_2.select.with_index { |pos, idx| idx==r }}: #{@feature}" #showing the gene coordinates of the Bio::Feature object

                    ###EX 4. Loop over each one of your CTTCTT features and create a GFF3-formatted file of these features###
                    File.open("gene_coordinate.gff", "a") { |f| f.write("#{seqid=geneid}\t#{source="."}\t#{type="insertion_site"}\t#{start=@feature.qualifiers[0].value[0]}\t#{final=@feature.qualifiers[0].value[1]}\t#{score="."}\t#{strand=@feature.qualifiers[1].value}\t#{phase="."}\t#{attributes="ID=#{geneid}"};repeat=#{type=@feature.qualifiers[0].qualifier}\n") }   
                    ###EX 5. Re-execute your GFF file creation so that the CTTCTT regions are now in the full chromosome coordinates used by EnsEMBL###
                    File.open("chr_coordinate.gff", "a") { |f| f.write("#{seqid=chr}\t#{source="."}\t#{type="insertion_site"}\t#{start=chr_pos.select.with_index { |pos, idx| idx==@feature.qualifiers[0].value[0]-1 }}\t#{final=chr_pos.select.with_index { |pos, idx| idx==@feature.qualifiers[0].value[1]-1 }}\t#{score="."}\t#{strand=@feature.qualifiers[1].value}\t#{phase="."}\t#{attributes="ID=#{geneid}"};repeat=#{type=@feature.qualifiers[0].qualifier}\n") }  
                  end
                end
                #puts fasta[geneid].features
                @repeats[geneid] = fasta[geneid]
                #puts @repeats[geneid]
                
              else #searching for exons on + strand
                #puts "\t\nEXON in forward strand:" + feature.position #showing gene coordinates of the exon
                pos1 = feature.position
                pos1_1 = pos1.match(/(?<start>\d+)..(?<end>\d+)/)
                pos1_2 = (pos1_1[:start].to_i..pos1_1[:end].to_i) #getting coordinates of the exon in BIO-informatics style (first nucleotide is '1')
                exon = fasta[geneid][pos1_2.begin-1..pos1_2.end-1] #changing coordinates to informatics style (first nucleotide is '0') to retrieve the exon's sequence  
                #puts exon #showing exon's sequence
                positions_f = exon.enum_for(:scan, rep_f_re).map { Regexp.last_match.begin(0) } #getting coordinates of "CTTCTT" repetitions inside the exon (not gene coordinates)     
                #puts "\texon positions with #{rep_f.to_s}: #{positions_f}" unless positions_f.empty? #showing the exon coordinates of the forward repeat
                positions_r = exon.enum_for(:scan, rep_r_re).map { Regexp.last_match.begin(0) } #getting coordinates of "AAGAAG" repetitions inside the exon (not gene coordinates)   
                #puts "\texon positions with #{rep_r.to_s}: #{positions_r}" unless positions_r.empty? #showing the exon coordinates of the reverse repeat
                
                
                
                if embl.sv =~ /chromosome:TAIR10:([^"]+):([^"]+):([^"]+):([^"]+)/ #searching for chromosomic coordinates for EX 5
                  chr = $1 #retrieving the chromosome
                  chr_pos = ($2.to_i..$3.to_i) #retrieving the chromosomic coordinates of the gene
                  
                  ###EX 3. Create a new Sequence Feature with the coordinates of every "CTTCTT" sequence and add it to the EnsEMBL Sequence object.### 
                  fasta[geneid].features = [] #creating an array to store the new feature for geneid's sequence (Bio::Sequence object)
                  positions_f.each do |f| ##creating a Bio::Feature object to store the information about the forward repeats
                    @feature = Bio::Feature.new('exon',pos1_2,
                      [ Bio::Feature::Qualifier.new("cttctt", [pos1_2.select.with_index { |pos, idx| idx==f }[0], pos1_2.select.with_index { |pos, idx| idx==f+5 }[0]] ), 
                        Bio::Feature::Qualifier.new('strand', '+')
                        ])
                    fasta[geneid].features << @feature #adding the new Bio::Feature object to the Bio::Sequence object
                    #puts "\t\tgene position #{pos1_2.select.with_index { |pos, idx| idx==f }}: #{@feature}" #showing the gene coordinates of the Bio::Feature object

                    ###EX 4a. Loop over each one of your CTTCTT features and create a GFF3-formatted file of these features###
                    File.open("gene_coordinate.gff", "a") { |f| f.write("#{seqid=geneid}\t#{source="."}\t#{type="insertion_site"}\t#{start=@feature.qualifiers[0].value[0]}\t#{final=@feature.qualifiers[0].value[1]}\t#{score="."}\t#{strand=@feature.qualifiers[1].value}\t#{phase="."}\t#{attributes="ID=#{geneid}"};repeat=#{type=@feature.qualifiers[0].qualifier}\n") }   
                    ###EX 5. Re-execute your GFF file creation so that the CTTCTT regions are now in the full chromosome coordinates used by EnsEMBL###
                    File.open("chr_coordinate.gff", "a") { |f| f.write("#{seqid=chr}\t#{source="."}\t#{type="insertion_site"}\t#{start=chr_pos.select.with_index { |pos, idx| idx==@feature.qualifiers[0].value[0]-1 }}\t#{final=chr_pos.select.with_index { |pos, idx| idx==@feature.qualifiers[0].value[1]-1 }}\t#{score="."}\t#{strand=@feature.qualifiers[1].value}\t#{phase="."}\t#{attributes="ID=#{geneid}"};repeat=#{type=@feature.qualifiers[0].qualifier}\n") }  
                  end
                  positions_r.each do |r| ##creating a Bio::Feature object to store the information about the reverse repeats
                    @feature = Bio::Feature.new('exon',pos1_2,
                      [ Bio::Feature::Qualifier.new("aagaag", [pos1_2.select.with_index { |pos, idx| idx==r }[0], pos1_2.select.with_index { |pos, idx| idx==r+5 }[0]] ),
                        Bio::Feature::Qualifier.new('strand', '-')
                        ])
                    fasta[geneid].features << @feature #adding the new Bio::Feature object to the Bio::Sequence object
                    #puts "\t\tgene position #{pos1_2.select.with_index { |pos, idx| idx==r }}: #{@feature}" #showing the gene coordinates of the Bio::Feature object

                    ###EX 4a. Loop over each one of your CTTCTT features and create a GFF3-formatted file of these features###
                    File.open("gene_coordinate.gff", "a") { |f| f.write("#{seqid=geneid}\t#{source="."}\t#{type="insertion_site"}\t#{start=@feature.qualifiers[0].value[0]}\t#{final=@feature.qualifiers[0].value[1]}\t#{score="."}\t#{strand=@feature.qualifiers[1].value}\t#{phase="."}\t#{attributes="ID=#{geneid}"};repeat=#{type=@feature.qualifiers[0].qualifier}\n") }   
                    ###EX 5. Re-execute your GFF file creation so that the CTTCTT regions are now in the full chromosome coordinates used by EnsEMBL###
                    File.open("chr_coordinate.gff", "a") { |f| f.write("#{seqid=chr}\t#{source="."}\t#{type="insertion_site"}\t#{start=chr_pos.select.with_index { |pos, idx| idx==@feature.qualifiers[0].value[0]-1 }}\t#{final=chr_pos.select.with_index { |pos, idx| idx==@feature.qualifiers[0].value[1]-1 }}\t#{score="."}\t#{strand=@feature.qualifiers[1].value}\t#{phase="."}\t#{attributes="ID=#{geneid}"};repeat=#{type=@feature.qualifiers[0].qualifier}\n") }  
                  end
                end
                #puts fasta[geneid].features
                @repeats[geneid] = fasta[geneid]
              end
            end
          end
        end
      end
    end
  end
  puts "\n\n\nDONE - now you can check your fasta sequences and exons with repeats above"
end


###EX 4a. Loop over each one of your CTTCTT features and create a GFF3-formatted file of these features###
def write_gene_gff_file(gfffile) #defining a function to delete the duplicated lines of the GFF file created with the previous function
  abort("Please, check the name of #{gfffile}, I can't find it!\n") unless (File.exists?(gfffile)) #aborting if the file does not exist
  gff = File.read(gfffile).split("\n")
  array = [] #creating an array to store the lines of the GFF file
  gff.each do |lines|
    array << lines
  end
  array = array.uniq #deleting duplicated lines

  array.each do |l|
    puts l
    File.open("gene_coordinates.gff", "a") { |f| f.write("#{l}\n") } #creating a new GFF file without the duplicated lines
  end

  File.delete(gfffile) #deleting the GFF file with duplicated lines
end


###EX 4b. Output a report showing which (if any) genes on your list do NOT have exons with the CTTCTT repeat###
def write_genes_wo_rep(agi_file,gff_file)
  abort("Please, check the name of #{agi_file}, I can't find it!\n") unless (File.exists?(agi_file)) #aborting if the agi_file does not exist
  abort("Please, check the name of #{gff_file}, I can't find it!\n") unless (File.exists?(gff_file)) #aborting if the gff_file does not exist
  gff_genes = [] #creating an array for the genes present in the GFF file
  gff = File.read(gff_file).split("\n") #reading the GFF file
  gff.each do |lines| #for each line of the GFF file do...
    gff_genes << lines[0..8] if lines =~ /A[Tt]\d[Gg]\d\d\d\d\d/ #retrieving the first 9 characters of each line if it contains an AGI code (i.e., ignoring the first line of the GFF file)
  end
  gff_genes = gff_genes.uniq #deleting duplicates of the same gene

  File.open("genes_without_rep.txt", "a") { |f| f.write("Genes without a CTTCTT repeat\n") } #creating a new file for the report with the genes without the "CTTCTT" repeat
  agi = File.read(agi_file).split #reading 'ArabidopsisSubNetwork_GeneList.txt'
  agi.each do |gene| #for each gene of the agi_file do...
    File.open("genes_without_rep.txt", "a") { |f| f.write("- #{gene}\n") } unless gff_genes.include? gene #adding the genes not present in the GFF file to the report with the genes without the "CTTCTT" repeat
  end
end


###EX 5. Re-execute your GFF file creation so that the CTTCTT regions are now in the full chromosome coordinates used by EnsEMBL###
def write_chr_gff_file(chrfile) #defining a function to delete the duplicated lines of the GFF file created with the previous function
  abort("Please, check the name of #{chrfile}, I can't find it!\n") unless (File.exists?(chrfile)) #aborting if the file does not exist
  gff = File.read(chrfile).split("\n")
  array = [] #creating an array to store the lines of the GFF file
  gff.each do |lines|
    array << lines.tr('[]', '')
  end
  array = array.uniq #deleting duplicated lines
  
  array.each do |l|
    puts l
    File.open("chr_coordinates.gff", "a") { |f| f.write("#{l}\n") } #creating a new GFF file without the duplicated lines
  end
  
  File.delete(chrfile) #deleting the GFF file with duplicated lines
end




## making files
read_file("ArabidopsisSubNetwork_GeneList.txt") #creating a variable with all the AGI codes of the list
retrieve_seq(@agi) #returns DNA sequences, exons, and more of "ArabidopsisSubNetwork_GeneList.txt" genes and creates a GFF file (use this GFF file as input for the following function)   
write_gene_gff_file("gene_coordinate.gff") #creates a GFF file without duplicated lines and with gene coordinates of the repeat "CTTCTT" 
write_genes_wo_rep("ArabidopsisSubNetwork_GeneList.txt","gene_coordinates.gff") #creates a txt file with the genes without "CTTCTT" repeat
write_chr_gff_file("chr_coordinate.gff") #creates a GFF file without duplicated lines and with chromosomic coordinates of the repeat "CTTCTT"


