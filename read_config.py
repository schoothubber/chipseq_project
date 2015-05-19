#!/usr/bin/env




def get_file_names(config_file):	
	"""
	Parse a chip-seq configuration text file.
	
	The text file contains 2 columns:
	column 1: BAM File
	column 2: Control File
	"""


	with open(config_file, 'r') as fo:
		data = fo.readlines()
		data = [line.split() for line in data]
	
	#print data
		
	return data
