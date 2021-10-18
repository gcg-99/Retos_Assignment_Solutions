#beginning a class definition for gene_information.tsv data
class Gene
  #creating "attribute accessors" to read and write objects' properties
  attr_accessor :id
  attr_accessor :name
  attr_accessor :mutant_phenotype
  attr_accessor :linked_to
  
  #initialization inside the class
  def initialize(id ,name ,mutant_phenotype , linked_to) 
    @id = id
    @name = name
    @mutant_phenotype = mutant_phenotype
    @linked_to = linked_to
  end

  
  def Gene.load_from_file() ##BONUS: the object have a #load_from_file($gene_information.tsv)
    require 'csv'
    gene_info = [] #creating a new array that will contain the 'gene_information.tsv' data without the header by rows
    CSV.foreach('gene_information.tsv', headers: true, :quote_char => "|") do |row|
      gene_info << row
    end
    
    #creating objects for 'Gene Class' and passing parameters for initialization
    @@all_genes = [] 
    for $x in 0..4
      if $x==0
        #formatting the information contained in the array with the 'gene_information.tsv' data
        gene_infox = gene_info[$x][0].split("\t")
        gene_infox.map!{ |element| element.gsub(/"/, '') }
        if gene_infox[0] =~ /A[Tt]\d[Gg]\d\d\d\d\d/ ##BONUS: testing the format of the Gene ID 
          @@all_genes << Gene.new(gene_infox[0],gene_infox[1],gene_infox[2],"none") #creating an array with the new Objects
        else
          puts "The Gene ID is not correct. Please, check it!" ##BONUS: rejecting incorrect Gene ID without crashing
        end
      else
        #formatting the information contained in the array with the 'gene_information.tsv' data
        gene_infox = "#{gene_info[$x][0]}, #{gene_info[$x][1]}"
        gene_infox = gene_infox.split("\t")
        gene_infox.map!{ |element| element.gsub(/"/, '') }
        if gene_infox[0] =~ /A[Tt]\d[Gg]\d\d\d\d\d/ ##BONUS: testing the format of the Gene ID 
          @@all_genes << Gene.new(gene_infox[0],gene_infox[1],gene_infox[2],"none") #creating an array with the new Objects
        else
          puts "The Gene ID is not correct. Please, check it!" ##BONUS: rejecting incorrect Gene ID without crashing
        end
        $x +=1
      end
    end
    #puts @@all_genes #showing the content of the created objects for 'Gene Class'
    #puts
  end
  
  
  def Gene.show_genes #defining a method to work with Gene Class objects
    return @@all_genes
  end
  

  def Gene.linked_genes #defining a method to write a final report with linked genes
    puts "Final Report:"
    Gene.show_genes.each do |g|
      if g.linked_to != "none"
        puts "#{g.name} is linked to #{g.linked_to}"
      end
    end
  end
  
end

