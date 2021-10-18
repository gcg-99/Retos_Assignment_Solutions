#beginning a class definition for seed_stock_data.tsv information
class Seed_stock < Gene #Seed_stock is inheriting the properties of Gene
  #creating "attribute accessors" to read and write objects' properties
  attr_accessor :stock
  attr_accessor :gene_ID
  attr_accessor :last_planted
  attr_accessor :storage
  attr_accessor :grams_remaining
  
  
  #initialization inside the class
  def initialize(stock ,gene_ID ,last_planted ,storage ,grams_remaining ) 
    @stock = stock    
    @gene_ID = gene_ID    
    @last_planted = last_planted
    @storage = storage
    @grams_remaining = grams_remaining
  end


  def Seed_stock.load_from_file() ##BONUS: the object have a #load_from_file($seed_stock_data.tsv)
    require 'csv'
    seed_stock_data = [] #creating a new array that will contain the 'seed_stock_data.tsv' data without the header by rows
    CSV.foreach('seed_stock_data.tsv', headers: true, :quote_char => "|") do |row|
      seed_stock_data << row
    end
    
    @@all_stocks = []
    for $s in 0..4
      seed_stock_datax = seed_stock_data[$s][0].split("\t") #formatting the information contained in the array with the 'seed_stock_data.tsv' data
      @@all_stocks << Seed_stock.new(seed_stock_datax[0],seed_stock_datax[1],seed_stock_datax[2],seed_stock_datax[3],seed_stock_datax[4]) #creating an array with the new Objects
      $s +=1
    end
    
    #introducing a Gene's Object in gene_ID property value by Gene ID (ATxGgxxxxx)
    Gene.show_genes.each do |y|
      @@all_stocks.each do |z|
        if y.id == z.gene_ID
        z.gene_ID = y
        #each Seed_stock Object will have the corresponding Gene's Object in its property @gene_ID
        end
      end
    end
    #puts @@all_stocks #showing the content of the created objects for 'Seed_stock Class'
    #puts
  end
  
  
  def Seed_stock.show_stocks #defining a method to work with Seed_stock Class objects
    return @@all_stocks
  end
  
  
  def Seed_stock.updating_seed_stock(grams_to_plant, date) #defining a method to update the amount of remaining grams in a specific date
    ##BONUS: the object have a #write_database('new_stock_file.tsv')
    Seed_stock.show_stocks.each do |s|
    updated_grams = s.grams_remaining.to_i - grams_to_plant
    s.last_planted = date
      if updated_grams <= 0
        updated_grams = 0
        puts "WARNING: we have run out of Seed Stock #{s.stock}"
        File.open("new_stock_file.tsv", "a") { |f| f.write("#{s.stock}\t#{s.gene_ID.id}\t#{s.last_planted}\t#{s.storage}\t#{updated_grams}\n") }
      else
        File.open("new_stock_file.tsv", "a") { |f| f.write("#{s.stock}\t#{s.gene_ID.id}\t#{s.last_planted}\t#{s.storage}\t#{updated_grams}\n") }
      end
    end
  end
  
  
  def Seed_stock.get_seed_stock(somestock) ##BONUS: objects may access individual SeedStock objects based on their stock ID (e.g., 'A334')
    Seed_stock.show_stocks.each do |s|
      if s.stock.include?(somestock)
        puts "Stock ID: #{s.stock}, Gene ID: #{s.gene_ID.id}, last planted: #{s.last_planted}, storage: #{s.storage}, grams remaining: #{s.grams_remaining}"
        break
      else
        puts "The Stock ID #{somestock} is not present in the database."
        break
      end
    end
  end
  
end
