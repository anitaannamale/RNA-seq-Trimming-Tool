#! /usr/bin/env python
# -*- coding: utf8 -*-

""" This script is a module for the main script that launch Trimmomatic (an 
	adapteur and quality trimming tool for high throughput sequencing data). 
	
	This module contains all functions to generate Trimmomatic commandline from 
	a argument given directly on console. 
	But it depend on 2 other modules : checking_entries(check entered arguments)
	and commandline (generate commandline) """

__author__ = "Anita Annamalé"
__version__  = "0.0.1"
__copyright__ = "copyleft"
__date__ = "2015/05"


#-------------------------- Modules Importation -------------------------------#


import argparse
import sys
import os.path
from argparse import RawTextHelpFormatter

from checking_entries import *
from commandline import *


#------------------------- Definition Of Functions ----------------------------#



def Trimmomatic_parser():
	"""
	Menu of the programme Trimmomatic (An adapter and Quality Trimming Tool).
	
	Don't takes any argument
	
	Return parser
	"""
	
	parser = argparse.ArgumentParser( 
		formatter_class=RawTextHelpFormatter,
		description='A python program that launch Trimmomatic Tool : A Adapter and Quality Trimming Tools for high throughput sequencing data',
		epilog="\nSome examples:\n\n\t\
'python Filtrage.py -h' \t for help\n\t\
'python Filtrage.py --XML' \t to launch trimmomatic using parameters from the XML file\n\t\
'python Filtrage.py PE read_1.fastq read_2.fastq -illuminaclip fasta-file.fa:2:10:30' \t for adapter trimming\n\t\
'python Filtrage.py SE read_1.fq -slidingwindow 10:30 -leading 30 -minlen 36' \t for quality trimming\n\t\
'python Filtrage.py PE read_1.fq.bz2 read_2.fq.bz2 -illuminaclip fasta-file.fa:2:10:30 -slidingwindow 10:30 -minlen 36' \t for adapter and quality trimming.")

	parser.add_argument("--XML", 
						action='store_const', 
						const='XML', 
						help='use the XML file or not"')
	
	parser.add_argument("layout", 
						type=str, 
						nargs='?',
						action='store', 
						help = "layout of reads.\n  Usage: 'SE' (single-ends) or 'PE' (paired-ends)")
	
	parser.add_argument("input", 
						type=str, 
						nargs= '*',
						action='store',
						help = "input read file.\n  Usage:\n  For SE 'readfile.fastq' || 'readfile.fq.gz'\n  For PE 'read_1.fastq read_2.fastq' || 'read_1.fastq.gz/bz2 read_2.fastq.gz/bz2'")
	
	parser.add_argument("-threads",
						type=int,
						action='store',
						default= 1,
						help = 'number of threads to use')
	
	parser.add_argument("-phred",
						type=int,
						action='store',
						help = 'phred quality of the data, 33 or 64')

	parser.add_argument("-illuminaclip",
						type=str,
						action='store',
						help= 'Adapter trimming.\n  Usage: -illuminaclip <fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>\n  Recommended: fasta-file:2:30:10')
	
	parser.add_argument("-slidingwindow",
						type=str,
						action='store',
						help= 'Quality trimming.\n  Usage: -slidingwindow <window-size>:<required-quality>\n  Recommended: 4:30 or 10:30')
	
	parser.add_argument("-maxinfo",
						type=str,
						action='store',
						help='Adaptative quality trimming.\n  Usage: -maxinfo <target-length>:<strictness>')
	
	parser.add_argument("-leading",
						type=int,
						action='store',
						help="minimum quality required to keep a base in 5' end.\n  Usage: '-leading <required-quality>'")
	
	parser.add_argument("-trailing",
						type=int,
						action='store',
						help="minimum quality required to keep a base in 3' end.\n  Usage: '-trailing <required-quality>'")
	
	parser.add_argument("-crop",
						type=int,
						action='store',
						help="number of base to remove from the end of the read.\n  Usage: '-crop <number>'")
	
	parser.add_argument("-headcrop",
						type=int,
						action='store',
						help="number of base to remove from the start of the read.\n  Usage: '-head <number>'")
	
	parser.add_argument("-minlen",
						type=int,
						action='store',
						help="minimum required length of a read to be kept.\n  Usage: '-minlen <length>'")
		
	parser.add_argument("-tophred33",
						action='store_const',
						const='TOPHRED33',
						help="(re)encodes the quality part of the FASTQ file to base 33.\n  Usage: '-tophred33'")
	
	parser.add_argument("-tophred64",
						action='store_const',
						const='TOPHRED64',
						help="(re)encodes the quality part of the FASTQ file to base 64.\n  Usage: '-tophred64'")
	
	return parser



def check_layout(arg):
	"""
	Function that check that a layout have been entered and it's either 'SE'(Single-Ends) 
	or 'PE'(Paired-Ends).
	
	Takes one argument :
		arg [dict] : dictionnary containning all the command argument entries.
	
	Returns:
		- layout [string] : clean layout ('PE' or 'SE')
		or quit the program
	"""
	
	# get de layout text
	layout = arg['layout']

	if(layout != None):
		# upper the layout
		layout = layout.upper()
		
		# if it's 'PE' or 'SE' return it or else quit
		if(layout == 'PE' or layout == 'SE'):
			return layout
		else :
			sys.exit("Error : Layout can only be 'SE' (for single-end) or 'PE' (for paired-ends).")



def nb_inputs(arg):
	"""
	Function that check the length of input entries depending on the layout
	
	Takes one argument :
		arg [dict] : dictionnary containning all the command argument entries.
	
	Returns:
		- 1 [integer] : if length of inputs are the expected one
		or
		quit
	"""
	
	# For Single Ends, one input file is expected
	if(arg['layout']=='SE'):
		if not len(arg['input']) == 1 :
			sys.exit("Error : For layout 'SE' you must give one read file.")

	# For Paired Ends, two input files are expected
	else:
		if not len(arg['input']) == 2 :
			sys.exit("Error : For layout 'PE' you must give two read files.")
	
	return 1



def check_input(arg):
	"""
	Function that check the extension of input files.
	
	Takes one argument :
		arg [dict] : dictionnary containning all the command argument entries.
	
	Returns:
		- 1 [integer] : if input files extension is '.fastq' or '.fq'
		- 2 [integer] : if input files extension is '.fastq/fq.bz2' or '.fastq/fq.gz'
		or
		quit
	"""
	
	# check each file given in input
	for files in arg['input']:
		check_extension(files) # from checking_entries module



def check_phred(arg):
	"""
	Function that check phred quality is 33 or 64.
	
	Takes one argument :
		args [dict] : dictionnary containning all the command argument entries.
	
	Returns:
		- 1 [integer] : if phred quality is 33 or 64
		or
		quit
	"""
	
	# get phred text
	phred = arg['phred']

	
	if phred != None :
		# if it's 33 or 64 return 1 or quit
		if not (phred == 33 or phred == 64):
			sys.exit("Error : Value for the option 'phred' can only be '33' or '64'.")
		else :
			return 1



def check_illuminaclip(arg):
	"""
	Function that check illuminaclip (adapter trimming) arguments.
	
	Takes one argument :
		arg [dict] : dictionnary containning all the command argument entries.
	
	Returns:
		- 1 [integer] : if all arguments are confrom
		or
		quit
	"""
	
	# get illumina clip arguments
	illum = arg['illuminaclip']

	if illum != None :
		
		# separates arguments
		illum = illum.split(':')
		
		# verify the number of arguments between ':'.
		if not len(illum) == 4:
			sys.exit("Error : Option illuminaclip must avec 4 elements between ':'")
		
		# get fasta file
		fasta_file=illum[0]
		
		# check fasta file extension
		ext = os.path.splitext(fasta_file)

		if not (ext[1] == '.fa' or ext[1] == '.fasta'):
			sys.exit("Error: Adapter file must be a fasta file.")
		
		
		# check that an integer is entered for seed mismatches
		if not illum[1].isdigit():
			sys.exit("Error: Value for seed mismatches must be an integer")

		# check that an integer is entered for palindrome clip threshold
		if not illum[2].isdigit():
			sys.exit("Error: Value for palindrome clip threshold must be an integer")

		# check that an integer is entered for single clip threshold
		if not illum[3].isdigit():
			sys.exit("Error: Value for simple clip threshold must be an integer")
		
		return 1
		
	else:
		return 1



def check_slidingwindow(arg):
	"""
	Function that check slidingwindow (quality trimming) arguments.
	
	Takes one argument :
		arg [dict] : dictionnary containning all the command argument entries.
	
	Returns:
		- 1 [integer] : if all arguments are confrom
		or
		quit
	"""
	
	# get slidingwindow parameter
	slidw = arg['slidingwindow']

	if slidw != None:
		
		# separates arguments
		slidw = slidw.split(':')
		
		# verify the number of arguments between ':'.
		if not len(slidw) == 2 :
			sys.exit("Error : Option silidingwindow must have 2 elements between ':'")

		# check that an integer is entered for window size
		if not slidw[0].isdigit():
			sys.exit("Error : Value for window-size (sliding-window) must be an integer")
	
		# check that an integer is entered for required quality
		if not slidw[1].isdigit():
			sys.exit("Error : Value for required-quality (sliding-window) must be an integer")
		
		return 1
	
	else :
		return 1
		

def check_maxinfo(arg):
	"""
	Function that check maxinfo (quality trimming) arguments.
	
	Takes one argument :
		arg [dict] : dictionnary containning all the command argument entries.
	
	Returns:
		- 1 [integer] : if all arguments are confrom
		or
		quit
	"""
	
	# get maxinfo parameters
	maxinfo = arg['maxinfo']

	if maxinfo != None : 
		
		# separates arguments
		maxinfo = maxinfo.split(':')
		
		# verify the number of arguments between ':'.
		if not len(maxinfo) == 2:
			sys.exit("Error : Option maxinfo must have two elements between ':'")

		# check that an integer is entered for targer length
		if not maxinfo[0].isdigit():
			sys.exit("Error : Value for target-length (maxinfo) must be an integer")


		# convert strictness value into float
		nb = float(maxinfo[1])
		# check if it's between 0 and 1
		if (nb < 0 or nb > 1):
			sys.exit("Error : Value for strictness must be a float between 0 and 1")
		
		return 1
	
	else :
		return 1



def check_args(arg):
	"""
	Function that check maxinfo (quality trimming) arguments.
	
	Takes one argument :
		arg [dict] : dictionnary containning all the command argument entries.
	
	Returns:
		arg [dict] or quit
	"""

	# check layout and add to args dictionnary	
	layout = check_layout(arg)
	arg['layout'] = layout

	# check number of input according to layout
	nb_inputs(arg)
	
	# check input files extension
	check_input(arg)

	# check phred quality
	check_phred(arg)

	# check illuminaclip arguments
	check_illuminaclip(arg)

	# check slidingwindow arguments
	check_slidingwindow(arg)

	# check maxinfo arguments
	check_maxinfo(arg)
		
	# save the place of working directory
	arg['output'] = '.'
	
	
	return arg


def argparse_commandline_step_1(arg,nb, inout):
	"""
	Function that generate step1 (adapter trimming) command line for Trimmomatic.
	
	Takes 3 argument :
		arg [dict] : dictionnary containning all the command argument entries.
		nb [integer] : number of executed command
		inout [dict] : dictionnary containning all the input and output files
	
	Returns:
		- cmd [string] : step1 command line
		or
		None : if no adapter trimming is done
	"""
	
	# base command line
	cmd = 'java -jar trimmomatic-0.33.jar'

	# adding layout and input and output files
	cmd, inout = commandline_input_output(arg, cmd, nb, inout) # from commandline module

	# if illuminaclip step, add argument to command line or else skip the step
	if(arg['illuminaclip'] != None):
		cmd += ' ILLUMINACLIP:{0}'.format(arg['illuminaclip'])

	else :
		cmd = None

	return cmd



def argparse_commandline_step_2(arg, nb, inout):
	"""
	Function that generate step1 2 (quality trimming) command line for Trimmomatic.
	
	Takes 3 argument :
		arg [dict] : dictionnary containning all the command argument entries.
		nb [integer] : number of executed command
		inout [dict] : dictionnary containning all the input and output files
	
	Returns:
		- cmd [string] : step 2 command line
		or
		None : if no quality trimming is done
	"""
	
	# base command line
	cmd = 'java -jar trimmomatic-0.33.jar'
	
	# adding layout and input and output files
	cmd, inout = commandline_input_output(arg, cmd, nb, inout)
	
	cmd_1 = cmd

	if(arg['slidingwindow'] != None):
		cmd += ' SLIDINGWINDOW:{0}'.format(arg['slidingwindow'])

	if(arg['maxinfo'] != None):
		cmd += ' MAXINFO:{0}'.format(arg['maxinfo'])

	if(arg['leading'] != None):
		cmd += ' LEADING:{0}'.format(arg['leading'])

	if(arg['trailing'] != None):
		cmd += ' TRAILING:{0}'.format(arg['trailing'])

	if(arg['crop'] != None):
		cmd += ' CROP:{0}'.format(arg['crop'])

	if(arg['headcrop'] != None):
		cmd += ' HEADCROP:{0}'.format(arg['headcrop'])

	if(arg['minlen'] != None):
		cmd += ' MINLEN:{0}'.format(arg['minlen'])

	if(arg['tophred33']):
		cmd += ' {0}'.format(arg['tophred33'])

	if(arg['tophred64']):
		cmd += ' {0}'.format(arg['tophred64'])

	if(cmd == cmd_1):
		cmd = None

	return cmd
	
	
	