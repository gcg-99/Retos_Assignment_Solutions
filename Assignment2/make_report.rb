=begin
Assignment 2. Intensive integration using Web APIs
Author: Gema Castillo García

FILES I HAVE USED:
      ArabidopsisSubNetwork_GeneList.txt

FILES I HAVE CREATED:
  - Scripts defining each Class:
      Interaction_Network.rb
      Uso_General_Annotation.rb

    **The 2 Classes are interconnected. The values of some Object Properties are other Objects: the networks are the key to link Interaction_Network > Uso_General_Annotation.**


  - The main script that uses the previous 2 scripts to make the report is:
      make_report.rb --> creates a new file called report.txt


### Final program to execute: $ ruby make_report.rb
=end



require 'rest-client'  

#creating a function called "fetch" that we can re-use everywhere in our code ()
def fetch(url, headers = {accept: "*/*"}, user = "", pass="") #the definition of fetch was written by Mark Wilkinson
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



#using classes
require './Interaction_Network.rb' #using Interaction_Network Class
require './Uso_General_Annotation.rb' #using Uso_General_Annotation Class



#creating networks
Interaction_Network.retrieve_networks #creating objects from 'ArabidopsisSubNetwork_GeneList.txt'
Interaction_Network.binding_networks #connecting the previous objects to make bigger networks and introducing them into new objects
Uso_General_Annotation.KEGG_annotation #retrieving KEGG pathways
Uso_General_Annotation.GO_annotation #retrieving GO terms
Uso_General_Annotation.object_annotation #creating new objects with the final networks and its annotations



##creating the final report with the networks and its annotations in report.txt

File.open("report.txt", "a") { |f| f.write("
### ========================================================================== REPORT ========================================================================== ###

Analysis of interactions between the 168 co-expressed genes of the list. It has been applied filters for species and MIscore cutoff to the results of BAR database:
    taxon ID for Arabidopsis thaliana = taxid:3702
    cutoff = 0.485 (value for optimal score predictions found by calcuating the maximal Matthews correlation coefficient (MCC) in https://europepmc.org/article/MED/25652942)

There are two kinds of interactions:
    direct: between two or more genes of the list
    indirect: two genes of the list has one or more interactors (not in the list) in common

Although some genes of the list are present in more than one of the following networks (probably because they form a protein complex with many proteins of the list), I have created NETWORKS THAT HAVE ONLY 2 GENES OF THE LIST that interact with each other directly or/and indirectly through other interactors.


### ========================================================================= NETWORKS ========================================================================= ###
") }



for i in 0..Interaction_Network.show_annotated_objects.length-1
  if Interaction_Network.show_annotated_objects[i].network_members.interactors == "direct"
    File.open("report.txt", "a") { |f| f.write("\n\nNETWORK ##{i+1}") }
    File.open("report.txt", "a") { |f| f.write("\nGene #{Interaction_Network.show_annotated_objects[i].network_members.agi_locus1} interacts directly with #{Interaction_Network.show_annotated_objects[i].network_members.agi_locus2}")}
    Interaction_Network.show_annotated_objects[i].kegg.each do |k|
      File.open("report.txt", "a") { |f| f.write("\n\tKEGG ID: #{k[0]}  Pathway Name: #{k[1]}") }
    end
    Interaction_Network.show_annotated_objects[i].go.each do |g|
      File.open("report.txt", "a") { |f| f.write("\n\tGO ID: #{g[0]} GO Term: #{g[1]}") }
    end
  else
    File.open("report.txt", "a") { |f| f.write("\n\nNETWORK ##{i+1}") }
    Interaction_Network.show_annotated_objects[i].network_members.interactors.each do |intn|
      if intn.class == Array
        intn.each do |intt|
          if Interaction_Network.show_annotated_objects[i].network_members.interactors[0].class == String
            File.open("report.txt", "a") { |f| f.write("\nGene #{Interaction_Network.show_annotated_objects[i].network_members.agi_locus1} interacts indirectly with #{Interaction_Network.show_annotated_objects[i].network_members.agi_locus2} by #{Interaction_Network.show_annotated_objects[i].network_members.interactors[0]}") }
            Interaction_Network.show_annotated_objects[i].kegg.each do |k|
              File.open("report.txt", "a") { |f| f.write("\n\tKEGG ID: #{k[0]}  Pathway Name: #{k[1]}") }
            end
            Interaction_Network.show_annotated_objects[i].go.each do |g|
              File.open("report.txt", "a") { |f| f.write("\n\tGO ID: #{g[0]} GO Term: #{g[1]}") }
            end
          else
            File.open("report.txt", "a") { |f| f.write("\nGene #{Interaction_Network.show_annotated_objects[i].network_members.agi_locus1} interacts indirectly with #{Interaction_Network.show_annotated_objects[i].network_members.agi_locus2} by #{Interaction_Network.show_annotated_objects[i].network_members.interactors[0].join(", ")}") }
          end
          break
        end
      elsif intn == "direct"
        File.open("report.txt", "a") { |f| f.write("\nGene #{Interaction_Network.show_annotated_objects[i].network_members.agi_locus1} interacts indirectly with #{Interaction_Network.show_annotated_objects[i].network_members.agi_locus2} by #{Interaction_Network.show_annotated_objects[i].network_members.interactors[0]}") }
        Interaction_Network.show_annotated_objects[i].kegg.each do |k|
          File.open("report.txt", "a") { |f| f.write("\n\tKEGG ID: #{k[0]}  Pathway Name: #{k[1]}") }
        end
        Interaction_Network.show_annotated_objects[i].go.each do |g|
          File.open("report.txt", "a") { |f| f.write("\n\tGO ID: #{g[0]} GO Term: #{g[1]}") }
        end
        break
      else
        File.open("report.txt", "a") { |f| f.write("\nGene #{Interaction_Network.show_annotated_objects[i].network_members.agi_locus1} interacts indirectly with #{Interaction_Network.show_annotated_objects[i].network_members.agi_locus2} by #{Interaction_Network.show_annotated_objects[i].network_members.interactors.join(", ")}") }
        Interaction_Network.show_annotated_objects[i].kegg.each do |k|
          File.open("report.txt", "a") { |f| f.write("\n\tKEGG ID: #{k[0]}  Pathway Name: #{k[1]}") }
        end
        Interaction_Network.show_annotated_objects[i].go.each do |g|
          File.open("report.txt", "a") { |f| f.write("\n\tGO ID: #{g[0]} GO Term: #{g[1]}") }
        end
        break
      end
    end
  end
end

puts "\nReport is ready! Look at 'report.txt' file\n"