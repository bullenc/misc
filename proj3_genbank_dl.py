################################################################################
# proj3_genbank_dl.py
# script to get rbcL/matK accessions for species from GenBank with e-utilities
# Need to have e-utilities installed; see link below
# https://www.ncbi.nlm.nih.gov/books/NBK25501/
#
# Run script twice for each gene; matK and rbcL
#
# useage: python3 proj3_genbank_dl.py [rbcL/matK] [list_of_species_names_file.csv]
#
# MAKE SEPARATE FOLDER FOR EACH GENE AND COPY SCRIPT AND SPECIES LIST THERE
#
# NOTE: if you find for one gene that there are species with no seqs at all in 
# genbank, you can remove them from the species list for the second gene to 
# speed up script runtime
#
# specify matK or rbcL, as well as provide list of species names as a file
# input file should have column headers removed and consist of one column
# and NOT be saved as a .txt file (i.e. csv or tsv is fine)
#
# NOTE: you cannot run a matK run and an rbcL run at the same time;
# NCBI will kick you off for too many requests at once
#
#
# Catie Ausland
# NIU
# Dec 2020
################################################################################

####import packages
import subprocess
import sys
import os
from collections import Counter
import re

####initialize directories
if not os.path.isdir('LISTS'):
	os.mkdir('LISTS')
if not os.path.isdir('DOWNLOADS'):
	os.mkdir('DOWNLOADS')

#read in command line arguments
gene = sys.argv[1]

if gene == 'rbcL':
	full_gene = "ribulose-1,5-bisphosphate carboxylase/oxygenase"
else:
	full_gene = "maturase K"

file = [i.strip() for i in open(sys.argv[2].strip()).readlines()] #input list

# first round of queries against GenBank
'''for i in file:
	#removing variety and subspecies to increase chance of returning hits
	#and then querying against GenBank
	
	if 'var.' in i or 'subsp.' in i: 
		
		#take binomial species name only
		j = " ".join(i.split(" ")[:2])
		
		#draft command
		command = 'esearch -db nuccore -query "' + gene + ' [PROT] OR ' + \
		full_gene + ' [PROT] AND ' + j + \
		' [ORGN]" | efetch -format docsum | xtract -pattern DocumentSummary \
		-element Title AccessionVersion Slen > ' + "_".join(i.split())+".txt"
		
		#call command
		subprocess.call(command, shell=True)
	
	else: #no var. or subsp. in name
		
		#draft command
		command = 'esearch -db nuccore -query "' + gene + ' [PROT] OR ' + \
		full_gene + ' [PROT] AND ' + i + \
		' [ORGN]" | efetch -format docsum | xtract -pattern DocumentSummary \
		-element Title AccessionVersion Slen > ' + "_".join(i.split())+".txt"
		
		#call command
		subprocess.call(command, shell=True)'''
				
# second round to redo those that were missed/file is empty from first round query
# grabbing ALL seqs for species in case gene was missed from initial query

'''file2 = [" ".join(i.replace(".txt", "").split('_')) for i in os.listdir(".") if \
os.path.getsize(i) == 0 and not os.path.isdir(i)]

for i in file2:
	#removing variety and subspecies to increase chance of returning hits
	#and then querying against GenBank
	if 'var.' in i or 'subsp.' in i:
	
		#take binomial species name only
		j = " ".join(i.split(" ")[:2])
		
		#draft command
		command = 'esearch -db nuccore -query " ' + j + \
		' [ORGN]" | efetch -format docsum | xtract -pattern \
		DocumentSummary -element Title AccessionVersion Slen > ' + \
		"_".join(j.split())+"_2.txt"
		
		#run command
		subprocess.call(command, shell=True)
	
	else:
	
		#draft command
		command = 'esearch -db nuccore -query " ' + i + \
		' [ORGN]" | efetch -format docsum | xtract -pattern \
		DocumentSummary -element Title AccessionVersion Slen > ' + \
		"_".join(i.split())+"_2.txt"
		
		#run command
		subprocess.call(command, shell=True)'''

###### record species with no sequences in genbank #######

# find empty files
empties = Counter([" ".join(i.replace("_2", "").split('_')) \
for i in os.listdir(".") if os.path.getsize(i) == 0])

# find empty entries after 2 rounds (i.e. no hits in GenBank at all)
no_seq_gb = sorted([i for i in empties.keys() if empties[i] == 2])

with open("LISTS/no_seq_genbank.txt", "w+") as w:
        w.write("\n".join(no_seq_gb))
w.close()

####### determine which are partial seqs, seqs from complete genomes #######
####### and species with no gene hit in Genbank #######

#initialize lists of accessions for download

#seqs from complete genome
if os.path.isfile("DOWNLOADS/"+gene+"_complete_genomes_download.txt"):
	os.remove("DOWNLOADS/"+gene+"_complete_genomes_download.txt")
comp_dl = open("DOWNLOADS/"+gene+"_complete_genomes_download.txt", "a+")

#parital seqs
if os.path.isfile("DOWNLOADS/"+gene+"_partials_download.txt"):
	os.remove("DOWNLOADS/"+gene+"_partials_download.txt")
partial_dl = open("DOWNLOADS/"+gene+"_partials_download.txt", "a+")

#no complete genomes or gene hits in genbank
if os.path.isfile("LISTS/"+gene+"_no_hits_genbank.txt"):
	os.remove("LISTS/"+gene+"_no_hits_genbank.txt")
no_hits_genbank = open("LISTS/"+gene+"_no_hits_genbank.txt", "a+") 

# find complete genomes or partials with rbcL in non-subsp./var. plants
# record those with complete genomes
# if no complete genome, determine if rbcL actually in sequences; 
# if so, record the longest one in partials file


non_sb_var = [i for i in os.listdir(".") if os.path.getsize(i) > 0 and \
'.txt' in i and 'subsp.' not in i and 'var.' not in i and not os.path.isdir(i)]


for i in non_sb_var:
	completes = [] # list to store complete genomes
	partials = [] # list to store partials/complete cds
	par_genome = [] # list to store 'partial genomes'
	
	species = " ".join(i.replace(".txt", "").split("_"))
	
	r = [j.strip().split("\t") for j in open(i).readlines()] # read in file
	
	for j in r:
		if 'complete genome' in j[0] and 'chloroplast' in j[0]:
			completes.append(j)
		elif gene in j[0] or full_gene in j[0]:
			partials.append(j)
		elif 'partial genome' in j[0]:
			par_genome.append(j)
			
	if len(completes) > 0: #check for any complete genomes
		
		if len(completes) > 1: #check if multiple complete genomes
			found = False
			for c in completes:
				if 'NC_' in c[1]: #check for reference genome
					store = c
					found = True
			if found == True: #if reference genome found, record accession
				comp_dl.write("\t".join([species]+store)+"\n")
			else: #no reference genome, but still multiple complete genomes;
			#determine which complete genome is the longest then
				longest = 0
				for c in completes:
					if int(c[-1]) > longest:
						longest = int(c[-1])
						store = c
				comp_dl.write("\t".join([species]+store)+"\n")
		
		#if only one complete genome found
		else:
			comp_dl.write("\t".join([species]+completes[0])+"\n")

	elif len(par_genome) > 0:
		if len(par_genome) == 1: #if one gene found
			comp_dl.write("\t".join([species]+par_genome[0])+"\n")
		else:
			longest = 0
			for c in par_genome:
				if int(c[-1]) > longest:
					longest = int(c[-1])
					store = c
			comp_dl.write("\t".join([species]+store)+"\n")

	elif len(partials) > 0: #find longest partials
		if len(partials) == 1: #if one gene found
			partial_dl.write("\t".join([species]+partials[0])+"\n")
		else: #if multiple partials, find longest one
			longest = 0
			for p in partials:
				if int(p[-1]) > longest:
					longest = int(p[-1])
					store = p
			partial_dl.write("\t".join([species]+store)+"\n")
	
	else: 
		#so at this point, these are essentially genomes that had no hits in first round,
		#and then had no hits in second round. 
		#Those that had no hits in '_2' round means no hits in first either 
		#(hence need for second round)
		#will record in 'no_hits_to_genbank.txt'
		no_hits_genbank.write(i.replace("_2", "").replace(".txt", "").replace("_", " ")+"\n")


# find complete genomes or partials with rbcL in subsp./var. plants

sb_var = [i for i in [i for i in os.listdir('.') if 'subsp.' in i or 'var.' in i]\
 if os.path.getsize(i) > 0 and '.txt' in i and not os.path.isdir(i)]


#list of plants with sequences found but no match to var./subsp.
if os.path.isfile("LISTS/"+gene+"_no_varsubsp_genbank.txt"):
	os.remove("LISTS/"+gene+"_no_varsubsp_genbank.txt")
no_varsubsp_genbank = open("LISTS/"+gene+"_no_varsubsp_genbank.txt", "a+") 

for i in sb_var:
	
	species = " ".join(i.replace(".txt", "").split("_"))
	sv = i.split("_")[3].replace(".txt", "") #record var./subsp.
	completes = [] # list to store complete genomes
	partials = [] # list to store partials/complete cds
	
	r = [j.strip().split("\t") for j in open(i).readlines()] # read in file
	
	#### check if subsp./var. found in seqs returned for those species ####
	
	#first check if subsp/var is the same as species name
	if not re.search(sv, i.split("_")[1], re.IGNORECASE):
			
		sv_hits = [] #hits to specific var/subsp
		
		for line in r:
			#if sv in line: ####REWRITE using regex instead!!!!!
			if re.search(sv, line[0], re.IGNORECASE):
				sv_hits.append(line)
	else:
		sv_hits = []
			
	if len(sv_hits) == 0:
		no_varsubsp_genbank.write(i.replace("_2", "")+"\n")
		
		
		#check for completes, genes and partial genomes of non-matches
		#record still under var/subsp from list though for recordkeeping
		completes = [] # list to store complete genomes
		partials = [] # list to store partials/complete cds
		par_genome = [] # list to store 'partial genomes'
		
		for j in r:
			if 'complete genome' in j[0] and 'chloroplast' in j[0]:
				completes.append(j)
			elif gene in j[0] or full_gene in j[0]:
				partials.append(j)
			elif 'partial genome' in j[0]:
				par_genome.append(j)
				
		if len(completes) > 0: #check for any complete genomes
		
			if len(completes) > 1: #check if multiple complete genomes
				found = False
				for c in completes:
					if 'NC_' in c[1]: #check for reference genome
						store = c
						found = True
				if found == True: #if reference genome found, record accession
					comp_dl.write("\t".join([species]+store)+"\n")
				else: #no reference genome, but still multiple complete genomes;
				#determine which complete genome is the longest then
					longest = 0
					for c in completes:
						if int(c[-1]) > longest:
							longest = int(c[-1])
							store = c
					comp_dl.write("\t".join([species]+store)+"\n")

			#if only one complete genome found
			else:
				comp_dl.write("\t".join([species]+completes[0])+"\n")
				
		elif len(par_genome) > 0:
			if len(par_genome) == 1: #if one gene found
				comp_dl.write("\t".join([species]+par_genome[0])+"\n")
			else:
				longest = 0
				for c in par_genome:
					if int(c[-1]) > longest:
						longest = int(c[-1])
						store = c
				comp_dl.write("\t".join([species]+store)+"\n")
				
		elif len(partials) > 0: #find longest partials
			if len(partials) == 1: #if one gene found
				partial_dl.write("\t".join([species]+partials[0])+"\n")
			else: #if multiple partials, find longest one
				longest = 0
				for p in partials:
					if int(p[-1]) > longest:
						longest = int(p[-1])
						store = p
				partial_dl.write("\t".join([species]+store)+"\n")

		else: 
			#so at this point, these are essentially genomes that had no hits in first round,
			#and then had no hits in second round. 
			#Those that had no hits in '_2' round means no hits in first either 
			#(hence need for second round)
			#will record in 'no_hits_to_genbank.txt'
			no_hits_genbank.write(i.replace("_2", "").replace(".txt", "").replace("_", " ")+"\n")

	else: # len(sv_hits) > 0:
		#go thru all of the sv_hits; determine if there are any hits to gene or complete genome
		completes = [] # list to store complete genomes
		partials = [] # list to store partials/complete cds
		par_genome = [] # list to store 'partial genomes'

		for j in sv_hits:
			if 'complete genome' in j[0] and 'chloroplast' in j[0]:
				completes.append(j)
			elif gene in j[0] or full_gene in j[0]:
				partials.append(j)
			elif 'partial genome' in j[0]:
				par_genome.append(j)
				
		if len(completes) > 0: #check for any complete genomes
		
			if len(completes) > 1: #check if multiple complete genomes
				found = False
				for c in completes:
					if 'NC_' in c[1]: #check for reference genome
						store = c
						found = True
				if found == True: #if reference genome found, record accession
					comp_dl.write("\t".join([species]+store)+"\n")
				else: #no reference genome, but still multiple complete genomes;
				#determine which complete genome is the longest then
					longest = 0
					for c in completes:
						if int(c[-1]) > longest:
							longest = int(c[-1])
							store = c
					comp_dl.write("\t".join([species]+store)+"\n")

			#if only one complete genome found
			else:
				comp_dl.write("\t".join([species]+completes[0])+"\n")
		
		elif len(par_genome) > 0:
			if len(par_genome) == 1: #if one gene found
				comp_dl.write("\t".join([species]+par_genome[0])+"\n")
			else:
				longest = 0
				for c in par_genome:
					if int(c[-1]) > longest:
						longest = int(c[-1])
						store = c
				comp_dl.write("\t".join([species]+store)+"\n")
				
		elif len(partials) > 0: #find longest partials
			if len(partials) == 1: #if one gene found
				partial_dl.write("\t".join([species]+partials[0])+"\n")
			else: #if multiple partials, find longest one
				longest = 0
				for p in partials:
					if int(p[-1]) > longest:
						longest = int(p[-1])
						store = p
				partial_dl.write("\t".join([species]+store)+"\n")
		
		else: 
			#so at this point, these are essentially genomes that had no hits in first round,
			#and then had no hits in second round. 
			#Those that had no hits in '_2' round means no hits in first either 
			#(hence need for second round)
			#will record in 'no_hits_to_genbank.txt'
			no_hits_genbank.write(i.replace("_2", "").replace(".txt", "").replace("_", " ")+"\n")
			#if not, then go thru whole file again and look for generalized hits


				

				
				
				
				
				