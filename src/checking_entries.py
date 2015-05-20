#! /usr/bin/env python
# -*- coding: utf8 -*-

""" This script is a module for the main script that launch Trimmomatic (an 
	adapteur and quality trimming tool for high throughput sequencing data). 
	
	This module contains all functions that checks XML file entries by the user 
"""

__author__ = "Anita Annamalé"
__version__  = "0.0.1"
__copyright__ = "copyleft"
__date__ = "2015/05"

#------------------------ Importation des modules -----------------------------#

import xml.etree.ElementTree as ET
import sys
import os.path

#------------------------ Définition des fonctions ----------------------------#



def no_blank(text) :
	"""
	Function that remove blank from a text.
	
	Takes one arguments : text [string]
	
	Returns one arguments : clean_text [string] : text without blank
	"""
	
	# remove blank (' ') in start and end of the text
	clean_text = text.strip()
	
	return clean_text



def not_empty(text) :
	"""
	Function that check if the text is empty or not.
	
	Takes one argument : text [string]
	
	Returns one argument :
		- clean_text [string] : a clean text without blank
		- 0 [integer] : if it's empty
	"""
	# checking that the text is not none
	if(text==None):
		return 0
	
	# removing blank from the text
	clean_text=no_blank(text)
	
	# checking if there is a text
	res = bool(clean_text)
	
	# convert the result into integer
	res = int(res)
	
	if(res==1) :
		return clean_text

	# return 0 if there was only blanks
	return res



def check_layout(text) :
	"""
	Function that check if the layout isn't empty, and if it's 'PE' or 'SE'.
	
	Takes one argument : text [string] : the layout
	
	Returns one argument :
		- layout [string] : if the layout is conform
		- or quit : if the layout isn't conform
	"""
	
	# check that the layout isn't empty

	layout = not_empty(text)

	if(layout):
		
		# upper the text
		layout = layout.upper()
		
		if(layout == 'SE' or layout == 'PE') :
			return layout
		
		else :
			sys.exit("/!\ The layout can only be 'PE' or 'SE'.")	

	else :
		sys.exit("/!\ You must enter a type of layout for your reads.")



def check_fastq_extension(text):
	""" 
	Booleen that checks if the extension of a input file is '.fastq' or '.fq'.

	Takes one argument : text [string]

	Returns one argument:
		- 1 [integer] : if it is
		- or quit : if the extension is neither '.fastq' nor '.fq'
	"""

	# getting the file extension
	ext = os.path.splitext(text)[1]

	if(ext=='.fastq' or ext=='.fq'):		
			return 1

	else:
		return 0



def check_extension(text):
	"""
	Booleen that check if a file is a fastq file and compressed (bz2 or gz) or not

	Takes one argument: text [string]

	Returns one argument :
		- 1 [integer] : if the extension is correct : it's a fastq file [.fastq or .fq] 
		- 2 [integer] : if the extension is correct : it's a compressed fastq file [fastq/fq.bz2 or fastq/fq.gz]
		- or quit : if the extension isn't known
	"""

	# check if the file have fastq extension, if yes return 1
	if(check_fastq_extension(text) == 1) : 
		return 1

	# split the file
	ext1 = os.path.splitext(text)
	
	# check the last extension (if compressed)
	if(ext1[1]=='.gz' or ext1[1]=='.bz2'):
		
		# check that the first extension is .fastq or .fq (check if it's a fastq file)
		if(check_fastq_extension(ext1[0]) == 1) : 
			return 2 
	
	else :
		sys.exit('Input files have not the right extension [fastq/fq.bzip2 , \
fastq/fq.gzip or fastq/fq]\n')



def check_input(text) :
	"""
	Function that check that the entry of the input file is not empty and that 
	it has the right extension

	Takes one argument: text [string]

	Returns one argument :
		- input_file [string] : if the extension is correct : it's a fastq file [.fastq or .fq] 
		- 0 [integer] : if the extension is correct : it's a compressed fastq file [fastq/fq.bz2 or fastq/fq.gz]
	"""	

	# check that it's not empty

	input_file = not_empty(text)

	if(input_file) :
		
		# check extensions of the file
		if(check_extension(input_file)):
			return input_file
		
	else :
		return 0



def check_child_number(tree, number):
	"""
	Boolen that checks that the tree have the expected number of child in the XML
	file
	
	Takes two arguments : - tree [ElementTree] : the tree
						  - number [integer] : the expected number of child
	
	Return one argument:
		- 1 [integer] : if the number find is equal to the expected one
		- 0 [integer] : else
	"""
	
	# getting the number of child
	nb_child = len(tree.getchildren())
	
	# checking if it's equal to the expected one
	if(nb_child == number): 
		return 1
	
	return 0
	


def get_pe_input(Paired_end, param):
	"""
	Function that gets the input and output file for paired-end (PE) data.
	
	Takes two arguments : - Paired_end [ElementTree] : the subtree Paired-end
						  - param [dict] : dictionnary containning all parameters
	
	Returns param[dict] where have been added inputs parameters for PE
	"""

	# find input files
	for filename in Paired_end.findall('input') :
		
		if(filename.get('name') == 'read 1') :
			
			# check the input file for read 1
			Read1 = check_input(filename.text)
				
			if not Read1: 
				sys.exit("/!\ You must give the fastq file containning all \
reads 1 for paired-end data")
		
		
		elif(filename.get('name') == 'read 2') :
			
			# check the input file for read 2
			Read2 = check_input(filename.text)

			if not Read2:	
				sys.exit("/!\ You must give the fastq file containning all \
reads 2 for paired-end data")
		
		else : 
			sys.exit("/!\ The value of 'name' in paired-end section have been \
modified")

	# add the checked input files 
	param['input']= Read1, Read2

	return param



def get_se_input(Single_end, param):
	"""
	Function that gets the input and output file for single-end (SE) data.
	
	Takes two arguments : - Single_end [ElementTree] : the subtree Single-end
						  - param [dict] : dictionnary containning all parameters
	
	Returns param[dict] where have been added inputs parameters for SE
	"""

	# find input file
	filename=Single_end.find('input').text
	
	filename = check_input(filename)

	# check the filename
	if(filename):
		# add to dic
		param['input'] = filename
		
	else :
		sys.exit("/!\ You must give a fastq file containing all reads \
for single-end data\n")
		
	return param



def check_skip(text,location):
	"""
	Function that check if the value for skip is not empty and if it is 'yes' or 'no'.
	
	Takes two arguments : - text [string] : the value for skip option
						  - location [string] : location of the skip option
	
	Returns :
		- clean_text : if the value is yes or no
		- or quit : if it's empty or not yes or no
	"""

	# check is the text is not empty
	checked_text = not_empty(text)

	if(checked_text):
		
		# lower the text
		checked_text = checked_text.lower()
		
		if(checked_text =='yes' or checked_text == 'no'):
			return checked_text
		
		else :
			sys.exit("/!\ Value for skip can only be 'yes' or 'no'.\n\
Please correct your skip value for %s." %location)
			
	else :
		sys.exit("/!\ Please enter a skip value for %s." %location)



def check_clip_file(text):
	"""
	Function that check if the skip file is not empty
	
	Takes one arguments : - text [string] : the file
	
	Returns :
		- clip_file : if a file have been entered
		- or quit : if no file have been given
	"""	

	# check that a file have been given
	clip_file = not_empty(text)

	if(clip_file) :
		return clip_file
				
	else : 
		sys.exit("/!\ You must enter a fasta file containing adapter sequences \
if you want to do adapter trimminig\n")



def check_integer(text, location):
	"""
	Function that check if the text isn't empty and is an integer
	
	Takes two arguments : - text [string] 
						  - location [string] : location of integer
	
	Returns :
		- checked_text : if it's an integer
		- or quit : if not
	"""	

	checked_text=not_empty(text)

	if(checked_text):

		# convert into integer
		checked_text=int(checked_text)
		
		return checked_text

	else :
		sys.exit("/!\ You haven't enter an integer for %s." %location )



def check_float(text, location):
	"""
	Function that check if the text isn't empty and is a float
	
	Takes two arguments : - text [string] 	
						  - location [string] : location of float
	
	Returns :
		- checked_text : if it's a float
		- or quit : if not
	"""		

	checked_text=not_empty(text)

	if(checked_text):
		# convert into float
		checked_text=float(checked_text)

		return checked_text

	else:
		sys.exit("/!\ You haven't enter a float for %s." %location)



def check_true_false(text, location):
	"""
	Function that check if the text isn't empty and is either 'true' either 'false'
	
	Takes two arguments : - text [string] 	
						  - location [string] : location of text
	
	Returns :
		- checked_text : if it's 'true' or 'false'
		- or quit : if not
	"""	

	checked_text = not_empty(text)
	
	if(checked_text):
		
		checked_text = checked_text.lower()
		
		if(checked_text == 'true' or checked_text == 'false'):
			return checked_text
		
		else :
			sys.exit("/!\ Value for %s can only be 'true' or 'false'" %location)
		
	else :
		sys.exit("/!\ You haven't enter a text for %s." %location)



def check_yes_no(text, location):
	"""
	Function that check if the text isn't empty and is either 'yes' either 'no'
	
	Takes two arguments : - text [string] 	
						  - location [string] : location of text
	
	Returns :
		- checked_text : if it's 'yes' or 'no'
		- or quit : if not
	"""	

	checked_text = not_empty(text)

	if(checked_text):

		checked_text = checked_text.lower()

		if(checked_text == 'yes' or checked_text == 'no'):
			return checked_text

		else :
			sys.exit("/!\ Value for show in %s can only be 'yes' or 'no'." %location) 

	else :
		sys.exit("/!\ You haven't enter a text for %s." %location)	



def get_min_adapt_length(skip,min_adapt_len):
	"""
	Function that get the argument for minimum adapter length in adapter trimming.
	
	Takes two arguments : - skip [string] 	
						  - min_adapt_len [integer]
	
	Returns :
		- 8 [integer] : argument by default
		- or minlen [integer] : value choosen by user
	"""		

	checked_skip = check_skip(skip, 'min-adapter-length in adapter trimming')

	if(checked_skip == 'yes') :
		return 8
	
	else :
		minlen = check_integer(min_adapt_len.find('value').text, 
			      'min-adapter-length in adapter trimming')

		return minlen



def get_keep(skip, keep_both_reads):
	"""
	Function that get the argument for keep both reads in adapter trimming.
	
	Takes two arguments : - skip [string] 	
						  - keep_both_reads [string]
	
	Returns :
		- 'true' [string] : argument by default
		- or keep [string] : value choosen by user
	"""		

	checked_skip = check_skip(skip, 'keep-both-reads in adapter trimming')
	
	if(checked_skip == 'yes'):
		return 'true'
	
	else :
		keep = check_true_false(keep_both_reads.find('value').text,
								'keep-both-reads in adapter trimming')
		return keep



def check_format(text, location):
	"""
	Function that check if the format wanted for compression is accepted
	
	Takes two arguments : - text [string] 	
						  - location [string] = location of the text
	
	Returns :
		- checked_text [string] : if the compression format is acceptable
		- or quit : if not
	"""	

	checked_text = not_empty(text)

	if(checked_text):

		checked_text = checked_text.lower()

		if(checked_text == 'bz2' or checked_text == 'gz'):
			return checked_text

		else :
			sys.exit("Value for format in %s can only be 'bz2' or 'gz'." %location)

	else :
		sys.exit("You haven't enter a text for %s." %location)



def get_file_prefix(text):
	"""
	Function that gets the prefix of a file without all extensions.

	Takes one argument : text [string]

	Returns one argument : prefix [string]
	"""

	tmp = check_fastq_extension(text)

	prefix = os.path.splitext(text)[0]

	if(tmp==1) : 
		return prefix
	
	else :
		prefix = os.path.splitext(prefix)[0]
		return prefix

	