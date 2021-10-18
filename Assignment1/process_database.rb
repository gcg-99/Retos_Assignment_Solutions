=begin
Assignment 1. Creating Objects
Author: Gema Castillo García

FILES I HAVE USED:
      gene_information.tsv
      seed_stock_data.tsv
      cross_data.tsv

FILES I HAVE CREATED:
  - Scripts defining each Class:
      Gene.rb
      Seed_stock.rb
      Cross.rb

    **The 3 Classes are interconnected. The values of some Object Properties are other Objects: the Gene ID is the key to link Gene > Seed_stock and
      the seed stock ID is the key to link Seed_stock > Cross_data**. 


  - The main script that uses the previous 3 scripts to update the genebank information and to determine which genes are genetically-linked:
      process_database.rb --> creates a new file called new_stock_file.tsv


### Final program to execute: $ ruby process_database.rb  gene_information.tsv  seed_stock_data.tsv  cross_data.tsv new_stock_file.tsv
=end


#using classes
require './Gene.rb'
require './Seed_stock.rb'
require './Cross.rb'


#creating objects from tsv files
Gene.load_from_file()
Seed_stock.load_from_file() #Seed_stock has inherited Gene objects in the @gene_ID property
Cross.load_from_file() #Cross has inherited Seed_stock objects in @parent_1 and @parent_2 properties


### ========================================= TASK #1 ============================================ ###

##BONUS: the object have a #write_database('new_stock_file.tsv')
File.open("new_stock_file.tsv", "a") { |f| f.write("Seed_Stock\tMutant_Gene_ID\tLast_Planted\tStorage\tGrams_Remaining\n") }
Seed_stock.updating_seed_stock(7, "16/10/2021")

##see "BONUS: objects may access individual SeedStock objects based on their stock ID" in the jupyter notebook
##(I have not put this BONUS in this file because it is not required in the output specified in the instructions)


### ========================================= TASK #2 ============================================ ###

#Chi-square test on the F2 cross data
Cross.chi_square

#display of genetically-linked genes
Gene.linked_genes

