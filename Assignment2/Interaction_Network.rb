#'ArabidopsisSubNetwork_GeneList.txt' contains 168 genes

#beginning a class definition to show the Interaction Network
class Interaction_Network
  #creating "attribute accessors" to read and write objects' properties
  attr_accessor :agi_locus1
  attr_accessor :agi_locus2
  attr_accessor :interactors
  
  #initialization inside the class
  def initialize(agi_locus1 ,agi_locus2 ,interactors ) 
    @agi_locus1 = agi_locus1 #this attribute is for the genes of the list that I have used to do the search in BAR database
    @agi_locus2 = agi_locus2 #this attribute is for the genes of the list that have some interactors in common with agi_locus1
    @interactors = interactors #this attribute is for the interactors of each gene that are not present in the list
  end
  
  
  def Interaction_Network.retrieve_networks #defining a method to store the interactors of each gene in an array called 'gene_interactors'
    puts "\nSearching for interactions in BAR database from UToronto...\n"
    agi = [] #creating a new array that will contain the genes' AGI of "ArabidopsisSubNetwork_GeneList.txt"
    @@networks = [] #creating an array to store the Interaction_Network Objects
    agi = File.read("ArabidopsisSubNetwork_GeneList.txt").split #retrieving the genes' AGI of "ArabidopsisSubNetwork_GeneList.txt"
    taxid = 'taxid:3702' #taxon ID for A.thaliana
    cutoff = 0.485 #MIscore cutoff value for optimal score predictions found by calcuating the maximal Matthews correlation coefficient (MCC) in https://europepmc.org/article/MED/25652942
    for $g in 0..agi.length-1 #for each gene of the list do...
      agi[$g]["T"]= "t" #changing 'T' for 't' to follow the BAR database nomenclature 
      gene_interactors = [] #creating an array for the AGI code of each interactor
      res = fetch("http://bar.utoronto.ca:9090/psicquic/webservices/current/search/query/#{agi[$g]}?format=tab25") #accessing genes URLs
      if res  #res is either the response object, or False, so you can test it with 'if'
        body = res.body  #getting the "body" of the response
        if body.empty?
          puts "There is no entry for #{agi[$g]} in BAR database from UToronto." #printing a friendly message if the gene has not any interactor
        else #formatting the URL content to retrieve the required information
          lines = res.body.split("\n")
          lines.each do |l|
            fields = l.split("\t")
            if fields[9].include? taxid and fields[10].include? taxid and fields[14][15..19].to_f > cutoff #selecting A.thaliana genes and interactions with optimal quality
              if fields[0] != fields[1] #not including the interactor if the protein interacts with itself
                if fields[2] =~ /tair:([^"]+)/ #searching for reg expressions in the first column of interactors
                  if $1 != agi[$g] #selecting only the interactors, NOT the query
                    gene_interactors << $1
                  end
                end
                if fields[3] =~ /tair:([^"]+)/ #searching for reg expressions in the second column of interactors
                  if $1 != agi[$g]  #selecting only the interactors, NOT the query
                    gene_interactors << $1
                  end
                end
              end
            end
          end
        end
      end
      gene_interactors = gene_interactors.uniq #removing repeated interactors
      #gene_interactors = agi & gene_interactors ##to remove AGI not present in "ArabidopsisSubNetwork_GeneList.txt"
                                                 ##(it is hidden because I am also considering the interactors that are not present in the list)  
      for $i in 0..gene_interactors.length-1
        #puts "#{agi[$g]} interacts with #{gene_interactors[$i]}"
        @@networks << Interaction_Network.new(agi[$g],"none",gene_interactors) #creating an array with the new Objects
      end
    end
  end
  
  
  def Interaction_Network.show_networks #defining a method to work with Interaction_Network Class objects
    return @@networks
  end
  
  
  def Interaction_Network.binding_networks #defining a method to connect different objects to make bigger networks
    puts "\nBinding networks...\n"
    @@all_networks = [] #creating an array to store the new Interaction_Network Objects
    common_interactors = [] #creating an array to store the genes of the list that interact with each other
    Interaction_Network.show_networks.each do |n1|
      Interaction_Network.show_networks.each do |n2|
        unless n1.agi_locus1 == n2.agi_locus1 #ignoring comparisons from a gene with itself
          unless n1.interactors & n2.interactors == [] #dismissing the genes of the list without interactors in common
            if n2.interactors.include? n1.agi_locus1 #searching direct interactions between the genes of the list
              common_interactors << [n1.agi_locus1, n2.agi_locus1, [n1.interactors & n2.interactors, "direct"]]
            else #searching indirect interactions between the genes of the list
              common_interactors << [n1.agi_locus1, n2.agi_locus1, n1.interactors & n2.interactors]
            end
          else #searching direct interactions between the genes of the list without interactors in common
            if n2.interactors.include? n1.agi_locus1
              common_interactors << [n1.agi_locus1, n2.agi_locus1, "direct"]
            end   
          end
        end 
  
      end
    end

    common_interactors = common_interactors.uniq #deleting the duplicated networks
    common_interactors2 = [] #creating an array to store definitive networks without duplicates
    common_interactors.each do |c|
      common_interactors2 << c.sort_by do |cc| #sorting by alphabet the members of each network to compare different networks and delete the duplicates
        cc.class == Array ? cc.to_s : cc
      end
    end
    
    common_interactors2 = common_interactors2.uniq
    common_interactors2.each do |ci|
      @@all_networks << Interaction_Network.new(ci[0],ci[1],ci[2]) #creating an array with the new Objects with bigger networks
    end
  end
  
  
  def Interaction_Network.show_all_networks #defining a method to work with Interaction_Network Class networks
    return @@all_networks
  end
  
  
end
