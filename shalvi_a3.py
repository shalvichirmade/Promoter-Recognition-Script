##Assignment 3 - BINF 6410
##By Shalvi Chirmade - December 13, 2021

##Searching for known motifs in the promoter region of coexpressed genes in Zea mays. 


#Created a file that contains the names of the fasta files and their corresponding gff3 file. This file can be edited to represent the fasta and gff3 files the user would like to analyze. In this way, the names of the files are not hard-coded into the script.

#The lines of the file are structured to look like this: fasta file:gff3 file
#I will extract the names from this file to use in this script.

#User to enter the file path containing all the required files for this script.
while True:
    file_path = input("\nPlease enter the filepath to the directory you will using for this script. Make sure to add a slash at the end.\nThis directory should include the following files:\nFasta-GFF3.txt, promoters.txt, zea_mays_genes.txt and all the corresponding FASTA and GFF3 files listed in Fasta-GFF3.txt.\nThese files should not be in separate folders; all the files should be in one folder and this is the directory path you will be inputting.\n\n")

    print("\nYour file path is", file_path)

    correct = input("\nPlease double-check your file path; it should end with a slash. Do you wish to continue? Enter y or n. ")

    if correct == "y":
        print("\nYou have entered the correct file path.\n")
        
        gene_file = input("\nPlease enter the file name contianing the list of coexpressed genes you wish to analyze.\nFor example, if your file name is zea_mays_genes.txt then please enter: zea_mays_genes.txt\n\n")
        motif_file = input("\nPlease enter the file name contianing the motif list you wish to analyze.\nFor example, if your file name is promoters.txt then please enter: promoters.txt\n\n")

        print("\nThe file names you have chosen are:", gene_file, "for your list of coexpressed genes\nand", motif_file, "for your motif list.")

        cont = input("\nPlease double check these names. Do you wish to continue? Enter y or n ")

        if cont == "y":
            print("\nYou have entered the correct file names. Let us continue.\n")
            break
        elif cont == "n":
            print("\nPlease re-enter your file names.")
        else:
            print("\nYou have not entered the approproate letter. Please try again.")

    elif correct == "n":
        print("\nPlease re-enter your file path.")
    else:
        print("\nYou have not entered the appropriate letter. Please try again.")

#Extract the FASTA and GFF3 file names in Fasta-GFF3.txt
fg_file_name = file_path + "Fasta-GFF3.txt"

fg_file = open(fg_file_name)
fasta_file_names = []
gff3_file_names = []

for line in fg_file.readlines():
    each_line = line.split(":")
    fasta_file_names.append(each_line[0])
    gff3_file_names.append(each_line[1].rstrip())

fg_file.close()

#Open each FASTA file and store its DNA sequence.
#Note: I use a Mac and I found an error on the ways some functions are run on Windows. When I first wrote my code, I did not have the variables, strList or lines, I had to add those because of the way Windows was utilizing the append function. On a Mac, my previous code read each FASTA file within two seconds but when my code was run on a Windows, it tooks over 40 minutes per FASTA file.. I found a solution; I had to add strList and create complete_dna using "".join instead of just using .append. Now, this code reads each FASTA file within two seconds on both Mac and Windows.

complete_dna = [] #A list where every element is the DNA sequence from each FASTA file
num = 1
print("This step takes a long time, please be patient.")

for name in fasta_file_names:
    fasta = open(file_path + name)  
    dna = ""
    lines = fasta.readlines() #Added variable 1 to work on Windows
    strList = [] #Added variable 2 to work on Windows
    for line in lines: #Before, this read: for line in fasta.readlines():
        if line.startswith(">"):
            continue
        else:
            strList.append(line.rstrip()) #Line 1 that had to be changed
            #dna += line.rstrip() #Before, this was the only line used and then complete_dna.append(dna)
    complete_dna.append("".join(strList)) #Line 2 that had to be changed
    print("\nChromosome", num, "FASTA done.")
    num += 1
    fasta.close()

#Open each GFF3 file and store lines for all genes.
complete_gff3 = [] #A list where every element is the gff3 file containing gene lines. Each element is a list, where each element is each gene line
num = 1
print("\nInputting GFF3 files.")

for name in gff3_file_names:
    gff3 = open(file_path + name)
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


#A function to read in files per line.
def read_file(file_name):
    file = open(file_name)
    file_list = []
    for line in file:
        file_list.append(line.rstrip())
    
    file.close()
    return (file_list)


#Read in the gene names and motif sequence files. 
genes_list = read_file(file_path + gene_file)
motif_list = read_file(file_path + motif_file)


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


#Check to see if there are any duplicated gene entries. If yes then there may be more than one entry for a particular gene. If not then there are no multiple TSSs for any of the Zea mays genes.
complete_gene_set = set() #Sets cannot have duplicates
for gene in complete_class_gff3:
    complete_gene_set.add(gene.name)

if (len(complete_gene_set) == len(complete_class_gff3)): #If the length of the set matches the number of complete genes, then there are no duplicates
    print("\nThere are no duplicates in our whole collection of genes from Zea mays; there are no multiple TSSs for any of our genes.\n")
else:
    print("\nThere are gene duplicates.")


#Use the names of the co-expressed genes to create a list only containing the GFF3 lines for those genes.
selected_genes_class_gff3 = [] #Make a list of only the genes present in the co-expressed gene list
#Be aware that the gene order in this list is not the same as the gene order in genes_list.

for object in complete_class_gff3:
    if (object.name in genes_list):
        selected_genes_class_gff3.append(object)



#Now we need to extract the promoter region for each gene of interest; this corresponds to 500 nucleotides upstream from the start of the gene. If the gene is on the + strand, 500 nucleotides is 500-start position. If the gene is on the - strand, 500 nucleotides is 500+end position.

#A function to reverse complement the - strand.
def reverse_complement(dna):
    complement_dict = {"A" : "T",
                        "T" : "A", 
                        "C" : "G", 
                        "G" : "C",
                        "N" : "N"}

    return("".join(complement_dict[base] for base in reversed(dna)))

#A function to see if there any Ns in the promoter sequence and only extract the downstream bases from the Ns.
def downstream_seq(dna):
    dna_length = len(dna)
    last_N = dna.rfind("N")
    new_dna = dna[last_N+1:dna_length]
    
    return (new_dna)

#A function to create the list of promoter sequences. The function takes into account if the gene is on the positive or negative strand. As the gene does not need to correspond to these sequences, we will not be addressing the gene names anymore.
def promoters(gff_list):
    promoter_list = []

    for object in gff_list:
        chromosome_number = object.chr

        if object.strand == "+":
            start_seq = object.start
            seq = complete_dna[chromosome_number-1][start_seq-501:start_seq-1]
            promoter_list.append(seq)
        
        elif object.strand == "-":
            start_seq = object.end
            seq = complete_dna[chromosome_number-1][start_seq:start_seq+500]
            rc_seq = reverse_complement(seq)
            promoter_list.append(rc_seq)
    
    return (promoter_list)


promoter_sequence = promoters(selected_genes_class_gff3)


#If there are any N's in the promoter sequences, the instructions ask us to disregard any bases before the N's and only carry on with the downstream bases.
iter = 0
for seq in promoter_sequence:
    if "N" in seq:
        promoter_sequence[iter] = downstream_seq(seq)
        iter += 1

#There are no N's in our selected promoters.


#Using the list of known promoter motifs, find the number of times each motif is seen in the list of promoter sequences. We were told to count every occurence; even if the motif is displayed multiple times in one sequence, count all of them, even overlapping sequences.

import re #Regular expression library

#A function to count the number of times the motif sequence is present and write the output to a new file.
def motif_output_file(m_list, pr_seq, file_name):
    
    motif_dict = {} #Make a dictionary counting the number of times the motif sequence is present

    for motif in m_list:
        count = 0
        for seq in pr_seq:
            match = str("(?=" + motif + ")") #To count all occurences of the motif
            count += len(re.findall(match, seq, re.I)) #To ignore cases
            
        motif_dict[motif] = count

    #Write the motif_dict into a file for submission.
    genes_output = open(file_name, "w")

    genes_output.write("Motifs\tCounts\n\n") #Title for output file
    for motif in motif_dict:
        genes_output.write(motif + "\t" + str(motif_dict[motif]) + "\n")

    genes_output.close()

    return (print("The file, " + file_name + " has been created."))


#Write the selected genes output to a new file for submission.
motif_output_file(motif_list, promoter_sequence, file_path + "Coexpressed_Genes_Output.txt")



##Randomly select 594 genes from complete_class_gff3 and conduct the same analysis. This number was chosen as the zea_mays_genes.txt files contains a list of 594 genes.

import random #The only way I found to be able to randomly select entries
random_genes_class_gff3 = random.sample(complete_class_gff3, 594) 

#Extract promoter sequences 
promoter_sequence = promoters(random_genes_class_gff3)

#Remove any Ns; choose downstream
iter = 0
for seq in promoter_sequence:
    if "N" in seq:
        promoter_sequence[iter] = downstream_seq(seq)
        iter += 1

#Write the random genes output to a new file for submission.
motif_output_file(motif_list, promoter_sequence, file_path + "Random_Genes_Output.txt")



##Output a file containing a brief comparison of the coexpressed and random genes.
comparison_expression_file = open(file_path + "Comparison_Expression.txt", "w")

comparison_expression_file.write("In the mutiple tries I have run, I have found that every time I chose a different selection of random genes, the number of motifs found are consistent with the values of the selected genes." + "\n\n" + "There are no motifs, that I have seen, that are notably over or under represented amongst the coexpressed genes in comparison to the random genes!" + "\n\n" + "For both the coexpressed and random gene analysis, motif AAAG had over 1000 hits while the others were either  0, single digits or less than 500.")



