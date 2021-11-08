##BONUS: this Class can hold any functional annotation as long as you create the corresponding attribute and function

require 'json'

class Uso_General_Annotation < Interaction_Network #Uso_General_Annotation is inheriting the properties of Interaction_Network
  #creating "attribute accessors" to read and write objects' properties
  attr_accessor :network_members 
  attr_accessor :kegg 
  attr_accessor :go 
  
  #initialization inside the class
  def initialize(network_members ,kegg ,go ) 
    @network_members = network_members #this attribute will contain the objects with the networks created in Interaction_Network Class
    @kegg = kegg #this attribute will contain the KEGG pathways of all members of each network
    @go = go #this attribute will contain the GO terms of all members of each network
  end
  
  
  def Uso_General_Annotation.GO_annotation
    puts "\nRetrieving GO terms...\n"
    goterms = [] #creating an array for the GO terms of each interactor
    @@allgoterms = [] #creating an array for the GO terms of each network
    Interaction_Network.show_all_networks.each do |net|
      res = fetch("https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id=#{net.agi_locus1}&style=raw") #accessing the URL of the first gene of the networks
      if res#res is either the response object, or False, so you can test it with 'if'
        body = res.body #getting the "body" of the response
        if body.empty?
          goterms = []
        else
          lines = res.body.split("\n")
          lines.each do |l|
            if l =~ /GO:([^"]+); P:([^"]+);/ #selecting the GO terms from the biological process part of the GO Ontology
              goterms << ["#{$1}", "#{$2}"] #storing the GO IDs and Terms, respectively
              goterms = goterms.uniq #deleting duplicated GO terms present in the same URL
            end
          end
        end
      end
      
      res = fetch("https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id=#{net.agi_locus2}&style=raw") #accessing the URL of the second gene of the networks
      if res #res is either the response object, or False, so you can test it with 'if'
        body = res.body #getting the "body" of the response
        if body.empty?
          goterms = []
        else
          lines = res.body.split("\n")
          lines.each do |l|
            if l =~ /GO:([^"]+); P:([^"]+);/ #selecting the GO terms from the biological process part of the GO Ontology
              goterms << ["#{$1}", "#{$2}"] #storing the GO IDs and Terms, respectively
              goterms = goterms.uniq #deleting duplicated GO terms present in the same URL
            end
          end
        end
      end
      
      #this block is very long because it has to access all the interactors (they are in a different format)
      if net.interactors.class == Array #accessing the URLs of the interactors in common of the networks
        net.interactors.each do |int|
          unless int.class == String
            int.each do |intt|
              res = fetch("https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id=#{intt}&style=raw") #accessing interactor's URLs
              if res #res is either the response object, or False, so you can test it with 'if'
                body = res.body #getting the "body" of the response
                if body.empty?
                  goterms = []
                else
                  lines = res.body.split("\n")
                  lines.each do |l|
                    if l =~ /GO:([^"]+); P:([^"]+);/ #selecting the GO terms from the biological process part of the GO Ontology
                      goterms << ["#{$1}", "#{$2}"] #storing the GO IDs and Terms, respectively
                      goterms = goterms.uniq #deleting duplicated GO terms present in the same URL
                    end
                  end
                end
              end
            end
          else
            unless int == "direct"
              res = fetch("https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id=#{int}&style=raw") #accessing interactor's URLs
              if res#res is either the response object, or False, so you can test it with 'if'
                body = res.body #getting the "body" of the response
                if body.empty?
                  goterms = []
                else
                  lines = res.body.split("\n")
                  lines.each do |l|
                    if l =~ /GO:([^"]+); P:([^"]+);/ #selecting the GO terms from the biological process part of the GO Ontology
                      goterms << ["#{$1}", "#{$2}"] #storing the GO IDs and Terms, respectively
                      goterms = goterms.uniq #deleting duplicated GO terms present in the same URL
                    end
                  end
                end
              end
            end
          end
        end
      end
      @@allgoterms << goterms #storing the GO terms of each network in an array
      goterms = []
    end
  end
  
  
  def Interaction_Network.show_go #defining a method to work with GO annotations
    return @@allgoterms
  end
  
  
  def Uso_General_Annotation.KEGG_annotation
    puts "\nRetrieving KEGG pathways...\n"
    kegg = [] #creating an array for the KEGG pathways of each interactor
    @@allkegg = [] #creating an array for the KEGG pathways of each network
    Interaction_Network.show_all_networks.each do |net|
      address = "http://togows.org/entry/kegg-genes/ath:#{net.agi_locus1}/pathways.json" #accessing the URL of the first gene of the networks
      response = RestClient::Request.execute(
        method: :get,
        url: address)
      kegg_data = JSON.parse(response.body)
      for elem in kegg_data[0].each
        if kegg_data[0].empty?
          kegg = []
        elsif kegg_data[0]
          kegg << [elem[0],elem[1]] #storing the KEGG IDs and Pathways, respectively
        end
      end
      
      address = "http://togows.org/entry/kegg-genes/ath:#{net.agi_locus2}/pathways.json" #accessing the URL of the second gene of the networks
      response = RestClient::Request.execute(
        method: :get,
        url: address)
      kegg_data = JSON.parse(response.body)
      for elem in kegg_data[0].each
        if kegg_data[0].empty?
          kegg = []
        elsif kegg_data[0]
          kegg << [elem[0],elem[1]] #storing the KEGG IDs and Pathways, respectively
        end
      end
      
      if net.interactors.class == Array
        net.interactors.each do |int|
          unless int.class == String
            int.each do |intt|
              address = "http://togows.org/entry/kegg-genes/ath:#{intt}/pathways.json" #accessing interactor's URLs
              response = RestClient::Request.execute(
                method: :get,
                url: address)
              kegg_data = JSON.parse(response.body)
              for elem in kegg_data[0].each
                if kegg_data[0].empty?
                  kegg = []
                elsif kegg_data[0]
                  kegg << [elem[0],elem[1]] #storing the KEGG IDs and Pathways, respectively
                end
              end
            end
          else
            unless int == "direct"
              address = "http://togows.org/entry/kegg-genes/ath:#{int}/pathways.json" #accessing interactor's URLs
              response = RestClient::Request.execute(
                method: :get,
                url: address)
              kegg_data = JSON.parse(response.body)
              for elem in kegg_data[0].each
                if kegg_data[0].empty?
                  kegg = []
                elsif kegg_data[0]
                  kegg << [elem[0],elem[1]] #storing the KEGG IDs and Pathways, respectively
                end
              end
            end
          end
        end
      end
      kegg = kegg.uniq #deleting duplicated KEGGs
      @@allkegg << kegg #storing the KEGG pathways of each network in an array
      kegg = []
    end
  end
  
  
  def Interaction_Network.show_kegg #defining a method to work with KEGG annotations
    return @@allkegg
  end
  
  
  def Uso_General_Annotation.object_annotation #defining a method to create objects with a network and its GO/KEGG annotations
    @@annotated_objects = [] #creating an array to store the Uso_General_Annotation Objects
    for i in 0..Interaction_Network.show_all_networks.length-1
      @@annotated_objects << Uso_General_Annotation.new(Interaction_Network.show_all_networks[i],Uso_General_Annotation.show_kegg[i],Uso_General_Annotation.show_go[i]) 
    end
  end
  
  
  def Interaction_Network.show_annotated_objects #defining a method to work with Uso_General_Annotation Class objects
    return @@annotated_objects
  end
  
  
end  
