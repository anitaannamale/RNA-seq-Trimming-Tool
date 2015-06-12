#! /usr/bin/env python
# -*- coding: utf8 -*-

""" This script is a module for the main script that launch Trimmomatic (an 
	adapteur and quality trimming tool for high throughput sequencing data). 
	
	This module contains all functions to generate Trimmomatic commandline from 
	a XML file but is depend on another module which parse the XML file : 
	parseXML."""

__author__ = "Anita Annamalé"
__version__  = "0.0.1"
__copyright__ = "copyleft"
__date__ = "2015/05"


#-------------------------- Modules Importation -------------------------------#


from parseXML import *

import os
import os.path


#------------------------- Definition Of Functions ----------------------------#


def commandline_input_output(param, cmd, nb, inout):
	"""
	Function that add to 'cmd' the input and output files commandline depending
	on if the user want to do one or two trimming (infomation given by 'nb').
	
	Takes 4 arguments :
		- param [dict] : dictionnary containning all parameters
		- cmd [string] : base command line of the programme (trimmomatic)
		- nb [integer] : one or two step of trimming
		- inout [dict] : dictionnary containning all generated files on the 
						user's working directory
	
	Returns:
		- cmd [string] : the command line with the input and output files
		- inout [dict] : containing the new files names
	"""

	# For Single Ends data
	if(param['layout'] == 'SE'):
		cmd += ' SE'
		cmd += ' -threads {0}'.format(param['threads'])

		# get the prefix of the input file and creating the new files name
		prefix = get_file_prefix(param['input'])
		trimmed = "{0}/trimmed_{1}.fastq".format(param['output'],prefix)
		
		# add the compression format if choosen
		if 'compress' in param :
			trimmed += ".{0}".format(param['compress'])

		# if it's a the first trimming
		if (nb == 0):
			cmd += ' {0} {1}'.format(param['input'], trimmed)
			
			# adding the input and output files to inout
			inout['input'] = param['input']
			inout['trimmed'] = trimmed	
		
		# if nb == 1, the first step (adapter trimming) have been done
		# SO input of step 2 are output of step1
		elif (nb == 1):
			cmd += ' {0} {1}'.format(inout['trimmed'], trimmed)
		
	
	# For Paired Ends data	
	elif(param['layout'] == 'PE') :
		cmd += ' PE'
		cmd += ' -threads {0} '.format(param['threads'])

		# get the prefix of the input file and creating the new files name
		prefix_1 = get_file_prefix(param['input'][0])
		prefix_2 = get_file_prefix(param['input'][1])

		trimmed_1 = "{0}/trimmed_{1}.fastq".format(param['output'],prefix_1)
		trimmed_2 = "{0}/trimmed_{1}.fastq".format(param['output'],prefix_2)
		
		single_1 = "{0}/single_{1}.fastq".format(param['output'],prefix_1)
		single_2 = "{0}/single_{1}.fastq".format(param['output'],prefix_2)


		# add the compression format if choosen
		if 'compress' in param :

			trimmed_1 += ".{0}".format(param['compress'])
			trimmed_2 += ".{0}".format(param['compress'])
			single_1 += ".{0}".format(param['compress'])
			single_2 += ".{0}".format(param['compress'])


		# if it's a the first trimming
		if(nb==0):
			cmd += ' {0} {1} {2} {3} {4} {5}'.format(param['input'][0], 
				param['input'][1], trimmed_1, single_1, trimmed_2, single_2)
			
			# adding the input and output files to inout
			inout['input'] = param['input']
			inout['trimmed'] = trimmed_1, trimmed_2
			inout['single'] = single_1, single_2
		
		
		# if nb == 1, the first step (adapter trimming) have been done
		# SO input of step 2 are output of step1	
		elif(nb==1):
			cmd += '{0} {1} {2} {3} {4} {5}'.format(inout['trimmed'][0],
													inout['trimmed'][1],
													trimmed_1, single_1,
													trimmed_2, single_2)

	return cmd,inout



def change_output_as_input(inout, param):
	"""
	Function that change step1 output files to step2 input files.
	
	Takes 2 arguments :
		- inout [dict] : dictionnary containning all generated files on the 
						user's working directory
		- param [dict] : dictionnary containning all parameters	
	
	Returns:
		inout [dict] : containing the new files names
	"""
		
	# For Single Ends data
	if(param['layout'] == 'SE'):
		
		# getting the prefix of the file
		filename = get_file_prefix(param['input'])
		
		# new temporary filename
		new_name = '{0}/tmp{1}.fastq'.format(param['output'], filename)
	
		if 'compress' in param:
			new_name += '.{0}'.format(param['compress'])

		# renaming step1 output file
		os.system("mv {0} {1}".format(inout['trimmed'], new_name))
		# replacing the filenames in the dict
		inout['trimmed'] = new_name


	# For Paired Ends data
	else :
		
		# getting the prefix of files to create new tmp files
		filename_1 = get_file_prefix(param['input'][0])
		filename_2 = get_file_prefix(param['input'][1])
		
		# creating the tmp filename
		new_name_1 = './tmp{0}.fastq'.format(filename_1)
		new_name_2 = './tmp{0}.fastq'.format(filename_2)
		
		if 'compress' in param:
			new_name_1 += '.{0}'.format(param['compress'])
			new_name_2 += '.{0}'.format(param['compress'])
		
		# renaming step1 output files.
		os.system('mv {0} {1}'.format(inout['trimmed'][0], new_name_1))
		os.system('mv {0} {1}'.format(inout['trimmed'][1], new_name_2))
		
		# replacing the filenames in the dict
		inout['trimmed'] = new_name_1, new_name_2
		
		# deleting step1 singleton read files 
		os.system('rm {0} {1}'.format(inout['single'][0], inout['single'][1]))


	# return inout with replaced step1 outputfiles
	return inout	



def commandline_adapter(param, cmd):
	"""
	Function that add to commandline 'cmd' adapter trimming parameters.
	
	Takes 2 arguments :
		- param [dict] : dictionnary containning all parameters
		- cmd [string] : base commandline for trimmomatic
	
	Returns:
		- cmd [string] : 'cmd' with illuminaclip arguments 
		or
		- None if there is no adapter trimming step
	"""
	
	if 'clip' in param :

		cmd += ' ILLUMINACLIP:{0}'.format(param['clip'])
	
	# if no adapter trimming step
	else :
		cmd = None
		
	return cmd



def commandline_quality(param,cmd):
	"""
	Function that add to commandline 'cmd' quality trimming parameters.
	
	Takes 2 arguments :
		- param [dict] : dictionnary containning all parameters
		- cmd [string] : base commandline for trimmomatic
	
	Returns:
		- cmd [string] : 'cmd' with illuminaclip arguments 
		or
		- None if there is no adapter trimming step
	"""
	
	# if atleast one quality trimming parameter is used, flag will be switched 
	# to 1
	
	flag = 0
	
	if 'Crop' in param :
		cmd += ' CROP:{0}'.format(param['Crop'])
		flag = 1

	if 'Headcrop' in param :
		cmd += ' HEADCROP:{0}'.format(param['Headcrop'])
		flag = 1

	if 'Lead' in param :
		cmd += ' LEADING:{0}'.format(param['Lead'])
		flag = 1

	if 'Tail' in param :
		cmd += ' TRAILING:{0}'.format(param['Tail'])
		flag = 1

	if 'SW' in param :
		cmd += ' SLIDINGWINDOW:{0}'.format(param['SW'])
		flag = 1

	if 'MI' in param :
		cmd += ' MAXINFO:{0}'.format(param['MI'])
		flag = 1

	if 'Minlen' in param :
		cmd += ' MINLEN:{0}'.format(param['Minlen'])
		flag = 1

	if 'Avgqual' in param :
		cmd += ' AVGQUAL:{0}'.format(param['Avgqual'])
		flag = 1
	
	
	if(flag == 0):
		cmd = None

	return cmd



def commandline_step_1(loc,param, nb,inout):
	"""
	Function that generate step 1 command line for Trimmomatic.
	
	Takes 3 arguments :
		- param [dict] : dictionnary containning all parameters
		- nb [integer] : number of executed step
		- cmd [string] : base commandline for trimmomatic
	
	Returns:
		- cmd [string] : the commandline for step1
		- inout [dict] : the new filenames generated by step1
	"""
	
	# base command line
	cmd = 'java -jar {0}trimmomatic-0.33.jar'.format(loc[:-3])
	
	# adding input and output filename
	cmd, inout = commandline_input_output(param, cmd, nb, inout)
	
	# adding adapter trimming parameters
	cmd = commandline_adapter(param, cmd)
	
	# add compression format if choosen
	if(cmd != None):
		if 'convert' in param :
			cmd += ' {0}'.format(param['convert'])
	
	
	return cmd, inout



def commandline_step_2(loc,param, nb, inout):
	"""
	Function that generate step 2 command line for Trimmomatic.
	
	Takes 3 arguments :
		- param [dict] : dictionnary containning all parameters
		- nb [integer] : number of executed step
		- cmd [string] : base commandline for trimmomatic
	
	Returns:
		- cmd [string] : the commandline for step2
	"""
	
	# base command line
	cmd = 'java -jar {0}trimmomatic-0.33.jar'.format(loc[:-3])

	# adding input and output filename
	cmd, inout = commandline_input_output(param,cmd,nb,inout)
	
	# adding quality trimming parameters
	cmd = commandline_quality(param, cmd)

	# add compression format if choosen
	if(cmd != None):
		if 'convert' in param :
			cmd += ' {0}'.format(param['convert'])

	return cmd


