#beginning a class definition for cross_data.tsv information
class Cross < Seed_stock #Cross is inheriting the properties of Seed_stock
  #creating "attribute accessors" to read and write objects' properties
  attr_accessor :parent_1
  attr_accessor :parent_2
  attr_accessor :f2_wild
  attr_accessor :f2_P1
  attr_accessor :f2_P2
  attr_accessor :f2_P1P2
  
  #initialization inside the class
  def initialize(parent_1 ,parent_2 ,f2_wild ,f2_P1 ,f2_P2 ,f2_P1P2 ) 
    @parent_1 = parent_1
    @parent_2 = parent_2
    @f2_wild = f2_wild
    @f2_P1 = f2_P1
    @f2_P2 = f2_P2
    @f2_P1P2 = f2_P1P2
  end
  
  
  def Cross.load_from_file() ##BONUS: the object have a #load_from_file($seed_stock_data.tsv)
    require 'csv'
    cross_data = [] #creating a new array that will contain the 'cross_data.tsv' data without the header by rows
    CSV.foreach('cross_data.tsv', headers: true, :quote_char => "|") do |row|
      cross_data << row
    end
    
    #creating objects for 'Cross Class' and passing parameters for initialization
    @@all_crosses = []
    for $c in 0..4
      cross_datax = cross_data[$c][0].split("\t") #formatting the information contained in the array with the 'cross_data.tsv' data
      @@all_crosses << Cross.new(cross_datax[0],cross_datax[1],cross_datax[2],cross_datax[3],cross_datax[4],cross_datax[5]) #creating an array with the new Objects
      $c +=1
    end
    
    #introducing a Seed_stock's Object in parent_1/2 properties values by seed stock ID (e.g., A348)
    Seed_stock.show_stocks.each do |yy|
      @@all_crosses.each do |zz|
        if yy.stock == zz.parent_1
          zz.parent_1 = yy
          #each Cross Object will have the corresponding Seed_stock's Object in its property @parent_1
        end
        if yy.stock == zz.parent_2
          zz.parent_2 = yy
          #each Cross Object will have the corresponding Seed_stock's Object in its property @parent_2
        end
      end
    end
    #puts @@all_crosses #showing the content of the created objects for 'Cross Class'
    #puts
  end

  
  def Cross.show_crosses #defining a method to work with Cross Class objects
    return @@all_crosses
  end
 
  
  def Cross.chi_square #defining a method to do the chi² test
    Cross.show_crosses.each do |c|
      total_F2 = c.f2_wild.to_i + c.f2_P1.to_i + c.f2_P2.to_i + c.f2_P1P2.to_i
      expected = [total_F2*(9.0/16), total_F2*(3.0/16), total_F2*(3.0/16), total_F2*(1.0/16)] #expected frequencies if not linked for dihybrids' (heterozygous for 2 characters) offspring
      chi² = ((((c.f2_wild.to_i - expected[0])**2)/expected[0]) + (((c.f2_P1.to_i - expected[1])**2)/expected[1]) + (((c.f2_P2.to_i - expected[2])**2)/expected[2]) + (((c.f2_P1P2.to_i - expected[3])**2)/expected[3]))
      if chi² > 7.8147 #results are statistically significant (p < 0.05) with df = 3 (4 phenotypes - 1)        
        c.parent_1.gene_ID.linked_to = c.parent_2.gene_ID.name
        c.parent_2.gene_ID.linked_to = c.parent_1.gene_ID.name
        puts "Recording: #{c.parent_1.gene_ID.name} is genetically linked to #{c.parent_2.gene_ID.name} with chisquare score #{chi²}"
      end
    end
  end
  
end


##BONUS: this Class represents the entire Seed Stock "database"
