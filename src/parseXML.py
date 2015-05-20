#! /usr/bin/env python
# -*- coding: utf8 -*-

""" This script is a module for the main script that launch Trimmomatic (an 
	adapteur and quality trimming tool for high throughput sequencing data). 
	
	This module contains all functions to parse a XML file but is depend on 
	another module which check the entries of user : checking_entries. """

__author__ = "Anita Annamalé"
__version__  = "0.0.1"
__copyright__ = "copyleft"
__date__ = "2015/05"

#-------------------------- Modules Importation -------------------------------#


import xml.etree.ElementTree as ET
import sys
from checking_entries import *


#------------------------- Definition Of Functions ----------------------------#


def parse_xml_file(filename_xml) :
	"""
	Function that parse an xml file and gets the root of the tree as an 
	ElementTree object.
	
	Takes one argument :
		filename_xml [string] : the xml file
	
	Returns:
		root [ElementTree] = the root of the tres
	"""
	
	# parse the xml file
	tree = ET.parse(filename_xml)
	
	# getting the root of the tree
	root = tree.getroot()
	
	return root



def separate_steps(root):
	"""
	Function that separate different steps (Inputs/Outputs and different Programs).
	
	Takes one argument : root [ElementTree] : the root of the tree
	
	Returns two arguments :
		- Puts [ElementTree] : subtree which contains inputs and outputs 
		  information's
		- Trimmomatic [ElementTree] : subtree which contains Trimmomatic 
		  program parameters
	"""
	
	# checking that the root have 2 child
	if not check_child_number(root,2) :
		sys.exit("/!\ Warning : The XML file must contain exactly one Input and \
Output section and one program (Trimmomatic)")
	
	# separting the different steps	
	for element in root :
		
		if(element.tag == 'input-output') :
			Puts = element
		
		elif(element.get('name')=='trimmomatic') :
			Trimmomatic=element
		
		else :
			sys.exit("/!\ Oops! Atleast the name of the section 'input-output' \
or the name of program (trimmomatic) have been modified")
	
	return Puts, Trimmomatic



def separate_categories_Trimmo(Trimmomatic):
	"""
	Function that separates the subtree (Trimmomatic program) according to categories.
	
	Takes one argument : Trimmomatic [ElementTree] : tree 
	
	Returns four arguments :
		- Adapter [ElementTree] : subtree which contains adapter trimming parameters
		- Quality [ElementTree] : subtree which contains quality trimming parameters
		- Useful [ElementTree] : subtree which contains uselful parameters
	"""
	
	# checking that the number of categories is 3   
	if not check_child_number(Trimmomatic,3) : 
		sys.exit("/!\ Warning : The XML file must contain exactly 3 categories \
for Trimmomatic!")
	
	# separates categories of Trimmomatic
	for category in Trimmomatic :
		
		# getting parameters for 'Adapter-Trimming' in a subtree
		if(category.get('name') == 'adapter-trimming') :
			Adapter = category
			continue
		
		# getting parameters for 'Quality-Trimming' in a subtree
		elif(category.get('name') == 'quality-trimming') : 
			Quality = category
			continue
			
		# getting parameters for 'Usefull-Parameters' in a subtree
		elif(category.get('name') == 'useful-parameters') :
			Useful = category
			continue
		
		else :
			sys.exit("/!\ Oops! Atleast one category haven't been recognized.\n\
Please, have a look at the name of Trimmomatic categories, atleast one of them \
have been modified.")
	
	return Adapter, Quality, Useful



def get_input_output_parameters(Puts, param):
	"""
	Function that gets inputs and outputs parameters
	
	Takes two arguments : - Puts [ElementTree] : subtree which contains inputs and
							outputs information's
						  - param [dict] : dictionnary containning all parameters
	
	Returns param[dict] where have been added inputs and outputs parameters
	"""
	
	# checking if the section 'input-output' contains 4 parameters
	if not check_child_number(Puts,4):
		sys.exit("/!\ Warning : The XML file must contain exactly 4 parameters \
in the section 'input-output'!")
	
	# get the layout text
	layout = Puts.find('layout').text
	
	# checking that the layout is not empty and is 'SE' or 'PE'
	checked_layout = check_layout(layout)
	
	# adding the layout to the dictionnary
	param['layout'] = checked_layout
	
	
	# getting input(s)

	if(checked_layout == 'SE') :
		# getting the subtree which contain SE input
		SE = Puts.find('single-ends')

		# getting the input file after checking
		param = get_se_input(SE, param)

	else:
		# getting the subtree which contain PE input
		PE = Puts.find('paired-ends')
		
		# getting the input files after checking
		param = get_pe_input(PE, param)


	# getting output directory
		# getting the working-directory where the new files must be generated
	directory = Puts.find('working-directory').text
	
	# checking if the working-directory is not empty
	if(not_empty(directory)):
		checked_directory = not_empty(directory)
		
		# adding into the dictionnary
		param['output']=checked_directory
	
	else :
		sys.exit("/!\ You must enter a 'working directory' where the new output\
 will be written")
	
	return param



def get_adapter_parameters(Adapter, param):
	"""
	Function that gets adapter trimming parameters
	
	Takes 2 arguments:
		- Adapter [ElementTree] : subtree which contains adapter trimming parameters
		- param [dict] : dictionnary containning all parameters

	Returns param [dict] with added adapter trimming parameters
	"""

	# checking that Adapter section have 2 child (skip and parameters)
	if not check_child_number(Adapter,2):
			sys.exit("/!\ Warning : The XML file must contain exactly 1 skip \
option and 1 parameter for adapter trimming")

	# getting the skip option text
	skip = Adapter.find('skip').text
	
	# checking if it's not empty and either yes either no
	checked_skip = check_skip(skip, 'adapter trimming')
	
	# if skip = 'no' getting the parameters
	if(checked_skip == 'no') :
		
		# getting the subtree parameter
		parameter = Adapter.find('parameter')
		
		# getting the parameter which name is 'illuminaclip'
		if(parameter.get('name')=='illuminaclip') :
			
			Clip = parameter
			
			# get obligatory parameters
			clip_file = Clip.find('adapters-fasta-file').text
			checked_clip_file = check_clip_file(clip_file)
			
			mismatches = check_integer(Clip.find('seed-mismatches').text, 
												'mismatches in illuminaclip')
			P_thres = check_integer(Clip.find('palindrome-clip-threshold').text,
									'palindrome clip threshold in illuminaclip')
			S_thres = check_integer(Clip.find('simple-clip-threshold').text,
									'simple clip threshold in illuminaclip')

			# optional parameters
				# get min-adapter-length value
			minlen = Clip.find('min-adapter-length')
			minlen_skip = minlen.find('skip').text
			min_length = get_min_adapt_length(minlen_skip,minlen)
			
				# get keep-both-reads value 
			keepreads = Clip.find('keep-both-reads')
			keepreads_skip = keepreads.find('skip').text
			keep = get_keep(keepreads_skip, keepreads)
			
			
			# adding all the parameters to the dictionnary
			param['clip']="{0}:{1}:{2}:{3}:{4}:{5}".format(checked_clip_file,mismatches,P_thres,S_thres,min_length,keep)
			

			return param

		# if parameter name is not 'illuminaclip'
		else :
			sys.exit("/!\ Warning : The name of parameter 'illuminaclip'\
have been modified or replaced by something else. Please rename it 'illuminaclip'")
	
	# if skip = 'no'	
	else :
		# return param without the argument for illuminaclip
		return param 



def get_quality_parameters(Quality, param) :
	"""
	Function that gets quality trimming parameters.

	Takes 2 arguments:
		- Quality [ElementTree] : subtree which contains quality trimming parameters
		- param [dict] : dictionnary containning all parameters

	Returns param [dict] with added quality trimming parameters
	"""

	# checking that quality trimming subtree have 9 child (skip and 8 parameters)
	if not check_child_number(Quality,9) :
			sys.exit("/!\ Warning : The XML file must contain exactly one skip \
option skip and 8 parameters for quality trimming.")
			
	# getting the skip option text
	skip = Quality.find('skip').text
	# checking that it is not empty and either 'yes' either 'no'
	checked_skip = check_skip(skip, 'quality trimming')
	
	# getting the parameters if skip = no
	if(checked_skip == 'no'):
		
		for parameter in Quality.findall('parameter'):
			
			# getting the skip option for each parameter
			param_skip = parameter.find('skip').text
			# checking the skip option text
			checked_param_skip = check_skip(param_skip, 
								'%s in quality trimming'%parameter.get('name'))
			
			# getting the arguments for the parameter and adding to dict if skip = no
			if(checked_param_skip == 'no'):
				
				if(parameter.get('name') == 'sliding-window'):
					SW_size = check_integer(parameter.find('window-length').text,
						      'window length in sliding window trimming')
					SW_quality = check_integer(parameter.find('required-quality').text,
								'required-quality for sliding window trimming')
					
					param['SW'] = '{0}:{1}'.format(SW_size,SW_quality)
					continue
				
				elif(parameter.get('name')== 'maxinfo'):
					MI_length = check_integer(parameter.find('target-length').text,
								'target-length in Maxinfo quality trimming')
					MI_strictness = check_float(parameter.find('strictness').text,
								'strictness in Maxinfo quality trimming')
					
					param['MI'] = '{0}:{1}'.format(MI_length, MI_strictness)
					continue
					
				elif(parameter.get('name') == 'leading'):
					lead_quality = check_integer(parameter.find('required-quality').text,
								   "required-quality in 'leading' quality trimming")
					
					param['Lead'] = lead_quality
					continue
		
				elif(parameter.get('name') == 'trailing'):
					tail_quality = check_integer(parameter.find('required-quality').text,
								   "required-quality in 'trailing' quality trimming")
					
					param['Tail'] = tail_quality
					continue
				
				elif(parameter.get('name') == 'crop'):
					crop_length = check_integer(parameter.find('length').text,
						"length in 'crop' quality trimming")
					
					param['Crop'] = crop_length
					continue
				
				elif(parameter.get('name') == 'headcrop'):
					headcrop_length = check_integer(parameter.find('length').text,
									  "length in 'headcrop' quality trimming")
					
					param['Headcrop'] = headcrop_length
					continue	
					
				elif(parameter.get('name') == 'minlen'):
					min_len = check_integer(parameter.find('length').text,
							 "length in 'minlen' quality trimming")
					
					param['Minlen'] = min_len
					continue
				
				elif(parameter.get('name') == 'average-quality'):
					avg_qual = check_integer(parameter.find('required-quality').text,
							   "required-quality in 'average quality' trimming")
					
					param['Avgqual'] = avg_qual
					continue
				
				else :
					sys.exit("You have modified a quality trimming parameter \
name or enter a new one which have not been recognized")
		
		return param

	# if skip = 'yes' return param(dict) without quality trimming parameters
	else :
		return param	



def get_useful_parameters(Useful, param) :
	"""
	Function that useful parameters.

	Takes 2 arguments:
		- Useful [ElementTree] : subtree which contains useful parameters
		- param [dict] : dictionnary containning all parameters

	Returns param [dict] with added quality trimming parameters
	"""
	
	# checking that useful parameter have 4 child
	if not check_child_number(Useful,4):
		sys.exit("/!\ Warning : The XML file must contain exactly 4 useful parameters\n")


	for parameter in Useful.findall('parameter') :
		
		# getting the name of the parameter
		if(parameter.get('name') == 'singleton-reads'):
			
			# checking if text for show is not empty and either yes or no
			show = check_yes_no(parameter.find('show').text, 'singleton-reads in useful-parameters')

			# adding to the dict if show = yes
			if(show == 'yes'):
				param['Show'] = 'yes'
			continue


		elif(parameter.get('name') == 'convert-to-phred') :
			
			# getting the skip text
			skip = parameter.find('skip').text

			# checking if the skip text is not empty and if it's yes or no
			checked_skip = check_skip(skip, 'convert-to-phred in useful parameters')

			if(checked_skip == 'no'): 
				format_num = check_integer(parameter.find('format').text,
						 'format in convert-to-phred in useful parameters.')
				
				# if format is phred33 or phred64, add to the dict
				if(format_num == 33 or format_num == 64):
					param['convert'] = "TOPHRED{0}".format(format_num)
					continue
				
				else :
					sys.exit("/!\ Warning : Quality score can only be converted \
to phred33 or phred64.")


		elif (parameter.get('name') == 'threads') :
			number = check_integer(parameter.find('number').text, 
					 'number in threads in useful parameters.')

			param['threads'] = number
			continue


		elif(parameter.get('name') == 'compressed-output'):
			
			skip = parameter.find('skip').text
			checked_skip = check_skip(skip, 'compressed-output in useful parameters.')
			
			if(skip == 'no') :
				comp = check_format(parameter.find('format').text, 
					   'format in compressed-output in useful parameters.')

				param['compress'] = comp
				continue

		else :
			sys.exit("You have modified a useful parameter name or enter a new \
one which have not been recognized\n")

	return param


