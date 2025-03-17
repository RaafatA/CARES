########################################
########################################
# Rafat Eissa
# Fash Diet Project 
# Date: 
########################################
########################################

# setwd("C:/Users/raafat.abdulmajeed/Documents/Fish_intestine")


# Either seqinr or Biostrings package should be installed for reading and writing fasta file
install.packages("seqinr", dependencies = TRUE)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install("Biostrings")

# or install Biostrings from bioconductor https://bioconductor.org/packages/release/bioc/html/Biostrings.html
# Now we show how to read the fasta file
# see https://github.com/ChiLiubio/file2meco to install file2meco
# rep_fasta_path <- system.file("extdata", "rep.fna", package="file2meco")

library(microeco)
setwd("C:/Users/raafat.abdulmajeed/Documents/Fish_intestine")
# Load the data 
samples_info = read.csv("FI.metadata-2.csv", row.names = 1)

taxonomy_table = read.csv("Results/denoise-dada2/taxonomy.biom/taxonomy.csv", row.names = 1)
env_data = read.csv("envs.csv", row.names = 1)
ASV_table = read.csv("Results/denoise-dada2/table-abundace-ASVs.biom/table-asv-abundance.csv", row.names = 1)
typeof(ASV_table)
rep_fasta <- seqinr::read.fasta(rep_fasta_path)
# or use Biostrings package
rep_fasta <- Biostrings::readDNAStringSet("Results/denoise-dada2/table-ASVs.biom/dna-sequences.fasta")
# try to create a microtable object with rep_fasta
data("otu_table_16S")

# In microtable class, all the taxa names should be necessarily included in rep_fasta
ASV_table <- ASV_table[rownames(ASV_table) %in% names(rep_fasta), ]
typeof(ASV_table)
test <- microtable$new(otu_table = ASV_table, rep_fasta = rep_fasta, 
                       tax_table = taxonomy_table, 
                       sample_table = samples_info)
test




library(microeco)

# load the example dataset from microeco package as there is a rep_fasta object in it
data(dataset)
dataset
# create a trans_func object
t1 <- trans_func$new(test)
# create a directory for result and log files
dir.create("test_prediction_17th March")
# https://chiliubio.github.io/microeco_tutorial/intro.html#tax4fun2 for installation
# ignore blast_tool_path parameter if blast tools have been in path
# the function can search whether blast tool directory is in the path, if not, automatically use provided blast_tool_path parameter
t1$cal_tax4fun2(blast_tool_path = "Ref/ncbi-blast-2.5.0+-x64-win64/ncbi-blast-2.5.0+/bin", path_to_reference_data = "Ref/Tax4Fun2_ReferenceData_v2",
                database_mode = "Ref99NR", path_to_temp_folder = "test_prediction_17th March")

# calculate functional redundancies
t1$cal_tax4fun2_FRI()


# prepare feature table and metadata
data(Tax4Fun2_KEGG)
# create a microtable object for pathways
func_mt <- microtable$new(otu_table = t1$res_tax4fun2_pathway, tax_table = Tax4Fun2_KEGG$ptw_desc, sample_table = test$sample_table)
func_mt$tidy_dataset()

# calculate relative abundances at three levels: Level.1, Level.2, Level.3
func_mt$cal_abund()

print(func_mt)

# bar plot at Level.1
tmp <- trans_abund$new(func_mt, taxrank = "Level.2", ntaxa = 15)
tmp$plot_bar(legend_text_italic = FALSE,facet = "conditon")

tmp$plot_bar(legend_text_italic = FALSE,facet = c("conditon", "group"))



# bar plot at Level.1
tmp <- trans_abund$new(func_mt, taxrank = "Level.2", ntaxa = 15, groupmean = "conditon")
tmp$plot_bar(legend_text_italic = FALSE)

# bar plot at Level.1
tmp <- trans_abund$new(func_mt, taxrank = "Level.3", ntaxa = 15)
tmp$plot_box(group='conditon', xtext_angle = 30)
tmp$plot_line(xtext_angle = 30)


tmp$plot_bar(legend_text_italic = FALSE,facet = c("conditon", "group"))


tmp <- trans_diff$new(dataset = func_mt, method = "rf", group = "group", taxa_level = "Level.3")
tmp$plot_line(plot_SE = TRUE)
tmp$plot_diff_bar(threshold = 3, width = 0.7)
