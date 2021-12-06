##Assignment 3 - BINF 6410
##By Shalvi Chirmade - December 13, 2021

##Searching for known motifs in the promoter region of coexpressed genes in Zea Mays. 


#Created a file that contains the names of the fasta files and their corresponding gff3 file. This file can be edited to represent the fasta and gff3 files the user would like to analyze. In this way, the names of the files are not hard-coded into the script.

#The lines of the file are structured to look like this: fasta file:gff3 file
#I will extract the names from this file to use in this script.

#User to enter the file path containing all the required files for this script.
while True:
    file_path = input("\nPlease enter the filepath to the directory you will using for this script. This directory should include the following files: Fasta-GFF3.txt, promoters.txt, zea_mays_genes.txt and all the corresponding FASTA and GFF3 files listed in Fasta-GFF3.txt\n\n")

    print("\nYour file path is", file_path)

    correct = input("\nPlease double-check your file path; it should NOT end with a backslash. Do you wish to continue? Enter y or n. ")

    if correct == "y":
        print("\nYou have entered the correct file path, we will continue.\n")
        break
    elif correct == "n":
        print("\nPlease re-enter your file path.")
    else:
        print("\nYou have not entered the appropriate letter. Please try again.")

#Extract the FASTA and GFF3 file names in Fasta-GFF3.txt
fg_file_name = file_path + "/Fasta-GFF3.txt"

fg_file = open(fg_file_name)
fasta_file_names = []
gff3_file_names = []

for line in fg_file.readlines():
    each_line = line.split(":")
    fasta_file_names.append(each_line[0])
    gff3_file_names.append(each_line[1].rstrip())

fg_file.close()

#Open each FASTA file and store its DNA sequence.
complete_dna = [] #A list where every element is the DNA sequence from each FASTA file
num = 1
print("This step takes a long time, please be patient.")

for name in fasta_file_names:
    fasta = open(file_path + "/" + name)  
    dna = ""
    for line in fasta.readlines():
        if line.startswith(">"):
            continue
        else:
            dna += line.rstrip()
    complete_dna.append(dna)
    print("\nChromosome", num, "FASTA done.")
    num += 1
    fasta.close()

#Open each GFF3 file and store lines for all genes.
complete_gff3 = [] #A list where every element is the gff3 file containing gene lines. Each element is a list, where each element is each gene line
num = 1
print("Inputting GFF3 files.")

for name in gff3_file_names:
    gff3 = open(file_path + "/" + name)
    genes = []

    for line in gff3.readlines():
        if line.startswith("#"):
            continue
        else:
            if (line.find("\tgene\t") > -1):
                genes.append(line)
    complete_gff3.append(genes)
    print("\nChromosome", num, "GFF3 done.")
    num += 1
    gff3.close()

#Read in the gene names and motif sequence files.
gene_name_file = open(file_path + "/zea_mays_genes.txt")
genes_list = [] #Make a list of the co-expressed genes
for line in gene_name_file:
    genes_list.append(line.rstrip())

gene_name_file.close()

motif_file = open(file_path + "/promoters.txt")
motif_list = [] #Make a list of the motif sequences
for line in motif_file:
    motif_list.append(line.rstrip())

motif_file.close()


#Create a GFF object to extract the features we require for further analysis. 
class GFF(object):    
    def __init__(self, data):  
        line_list = data.split("\t")
        self.chr = int(line_list[0]) #chromosome number
        self.start = int(line_list[3])
        self.end = int(line_list[4])
        self.strand = str(line_list[6])
        long_name = line_list[8].split(";")[0] #the first part of the name
        self.name =  long_name.split(":")[1] #extracting the exact gene name as written in zea_mays_genes.txt

#Use GFF class to extract this information for every gene in complete_gff list.
complete_class_gff3 = [] #Create a list for every gene in all the gff3 files

for file in complete_gff3: #Remember this was a list of lists
    for gene in file:
        complete_class_gff3.append(GFF(gene))


#Use the names of the co-expressed genes to create a list only containing the GFF3 lines for those genes.
selected_genes_class_gff3 = [] #Make a list of only the genes present in the co-expressed gene list
#Be aware that the gene order in this list is not the same as the gene order in genes_list.

for value in complete_class_gff3:
    if (value.name in genes_list):
        selected_genes_class_gff3.append(value)

#Check to see if there any duplicated genes in your list.
if len(genes_list) == len(selected_genes_class_gff3):
    print("\nThere are no multiple TSSs for your genes.")
elif len(genes_list) > len(selected_genes_class_gff3):
    print("\nThere are some genes from your list that have not been found. Please make sure you have included the correct gene names and restart this program.")
else:
    print("\nYou may have multiple TSSs for a few genes. Please choose the most upstream one.")




    


