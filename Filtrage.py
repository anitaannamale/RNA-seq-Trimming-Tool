#! /usr/bin/env python
# -*- coding: utf8 -*-

""" This script launch Trimmomatic, an adapteur and quality trimming tool for 
	high throughput sequencing data. 

	This script need three modules to function : parseXML, commandline and 
	argparse_commandline """

__author__ = "Anita Annamalé"
__version__  = "0.0.1"
__copyright__ = "copyleft"
__date__ = "2015/05"


#-------------------------- Modules Importation -------------------------------#

import argparse
import shlex, subprocess
import os
import sys

# Give the location of RNA-seq-Trimming-Tool directory
loc="./src"
sys.path.append(loc)

# personal modules
from parseXML import *
from commandline import *
from argparse_commandline import *

#------------------------- Definition Of Functions ----------------------------#


if __name__ == '__main__' :

	# parsing the commandline arguments.	
	parser=Trimmomatic_parser()	
	
	# getting commandline entries
	result = parser.parse_args()
	# change into dictionnary
	arguments = dict(result._get_kwargs())
	
	
	# checking :
	if len(sys.argv) < 2 :
		sys.exit("Usage : python Filtrage.py --XML or python Filtrage.py\
 layout read_files\nFor more informations use the option '-h'.")

	# if use the commandline argument then:
	if arguments['layout'] != None :
		if len(sys.argv) == 2 :
			sys.exit("Error : A read file must be specified if you\
 use the layout 'SE' or 2 read files if the layout 'PE' is specified")


	# Use XML file to launch Trimmomatic
	if(arguments['XML'] != None):
		
		# parse XML file
		config=parse_xml_file("configuration.xml")
		
		# separate sub trees
		in_out, trimmo = separate_steps(config)
		
		# separate trimming categories of Trimmomatic
		adapter, quality, useful=separate_categories_Trimmo(trimmo)

		# create an empty dict()
		param={} 
		
		# fill param
		param = get_input_output_parameters(in_out, param)
		param = get_adapter_parameters(adapter, param)
		param = get_quality_parameters(quality, param)
		param = get_useful_parameters(useful, param)

		# initializing nb (nb of exécuted commandline)
		nb = 0 
		
		# creating io [dict()] which will contain created files.
		io = dict()
		
		# generating step1 (adapter trimming) commandline
		cmd_step1, inout = commandline_step_1(loc,param,nb,io)
		
		# If Adapter Trimming
		if(cmd_step1 != None):
			
			# split the commandline
			args_1 = shlex.split(cmd_step1)
			
			# launch step1 
			with open("output_file_step1.out","w") as file_out1:
				prog_1 = subprocess.check_call(args_1, stderr=file_out1)
			
			# nb become 1 (first step done)
			nb = 1
	
		
		# If Quality Trimming Step
		tmp_cmd="test"
		if(commandline_quality(param,tmp_cmd)!=None):
			
			if(nb==1):
				# change step1 output files to step2 input files
				io = change_output_as_input(io, param)
		
			# generating step2 (quality trimming) commandline
			cmd_step2 = commandline_step_2(loc,param,nb,io)	
			
			# split commandline
			args_2 = shlex.split(cmd_step2)
		
			# launch step2
			with open("output_file_step2.out","w") as file_out2:
				prog_2 = subprocess.check_call(args_2, stderr=file_out2)
			
			if(nb==1):
				# delete temporary files
				os.system('rm {0} {1}'.format(io['trimmed'][0], io['trimmed'][1]))
	
	
	# Use Arguments line to launch Trimmomatic
	else:

		# check given arguments
		arguments=check_args(arguments)
		
		# initializing nb to 0
		nb = 0
		
		# creating io [dict()] which will contain created files.
		io= dict()
		
		# generating step1 (adapter trimming) commandline
		cmd_step1 = argparse_commandline_step_1(loc,arguments, nb,io)
	
		# If Adapter Trimming
		if(cmd_step1 != None):
			
			# split the commandline
			args_1 = shlex.split(cmd_step1)
			
			# launch step1 
			with open("output_file_step1.out","w") as out:
				prog_1 = subprocess.check_call(args_1, stderr=out)
			
			# nb become 1 (first step done)
			nb=1

		
		# If Quality Trimming step
		tmp_cmd="test"
		if(argparsecmd_quality(arguments,tmp_cmd)!=None):
			
			if(nb==1):	
				# change step1 output files to step2 input files
				io = change_output_as_input(io, arguments)
			
			# generating step2 (quality trimming) commandline
			cmd_step2 = argparse_commandline_step_2(loc,arguments, nb, io)
				
			# split commandline
			args_2 = shlex.split(cmd_step2)

			# launch step2
			with open("output_file_step2.out","w") as out:
				prog_2 = subprocess.check_call(args_2, stderr=out)
			
			if(nb==1):	
				# delete temporary files
				os.system('rm {0} {1}'.format(io['trimmed'][0], io['trimmed'][1]))
