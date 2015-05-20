#!/usr/bin/env

from subprocess import Popen, PIPE
import os

#test values
#bam_folder = "/home/wouter/chipseq_project/BAM_files/"
#aln_folder = "/home/wouter/chipseq_project/ALN_files/"


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


def remove_duplicates(bam_test_file, bam_folder):
	"""
	Remove or mark the duplicates in the bam files
	using sambamba
	"""
	#make the filenames including the full path
	bam_in = "%s%s"%(bam_folder, bam_test_file)
	dedup_bam_out = "%s%s%s"%(bam_folder, "dedupped_", bam_test_file)
	
	cmd = ["sambamba", "markdup", bam_in, dedup_bam_out]
	
	print " Removing duplicates from %s"%bam_in
	
	pipe = Popen(cmd, stdout=PIPE, stderr=PIPE)
	stdout, stderr = pipe.communicate()	


def converter(bam_test_file, bam_control_file, bam_folder, aln_folder):
	"""
	Pass a bam file to the actual process of converting
	Check for each bam file to be processed if a corresponding
	aln file is already present in the output file.
	"""
	
	
	#make the filenames including the full path
	bam_test = "%s%s%s"%(bam_folder, "dedupped_", bam_test_file)
	bam_ctrl = "%s%s"%(bam_folder, bam_control_file)
	
	#take everything except the .bam extension
	base_name_test = os.path.splitext(bam_test_file)[0]
	base_name_ctrl = os.path.splitext(bam_control_file)[0]
	
	#set aln output name
	aln_out_test = "%s%s.aln"%(aln_folder, base_name_test)
	aln_out_ctrl = "%s%s.aln"%(aln_folder, base_name_ctrl)
	
	#make a filelist.txt for seqpeak
	aln_files = [aln_out_test, aln_out_ctrl]
	make_seqpeaklist(base_name_test, aln_files, aln_folder)
	
	#check if this aln TEST file is already present in the aln folder
	if not os.path.isfile(aln_out_test):
		print "start BAM to ALN conversion for %s"%aln_out_test
		bam_to_aln_maker(bam_folder, bam_test, aln_folder, aln_out_test)
	else:
		print "%s is already present, skipping conversion..."%aln_out_test
		pass
		
	#check if this aln CONTROLfile is already present in the aln folder
	if not os.path.isfile(aln_out_ctrl):
		print "start BAM to ALN conversion for %s"%aln_out_ctrl
		bam_to_aln_maker(bam_folder, bam_ctrl, aln_folder, aln_out_ctrl)
	else:
		print "%s is already present, skipping conversion..."%aln_out_ctrl
		pass
		
		

def bam_to_aln_maker(bam_folder, bam_in, aln_folder, aln_out):
	"""
	Take all bam files
	Take a bam folder 
	Convert the bam files into aln files
	Store the ALN files in the ALN_folder
	
	NOTE: The BAM files and their corresponding control BAM files will
	be given by a text file containing 2 columns:
	column 1: BAM File
	column 2: Control File
	"""
	
	#set the Bedtools command line
	pipe_cmd = ["bedtools", "bamtobed", "-i", bam_in]
	

	#create a separate aln file for each process 
	fo = open(aln_out, 'w')
	print "created %s"%aln_out

	#run bedtools on each bamfile
	#and extract data for aln files
	proc = Popen(pipe_cmd, stdout=PIPE)
	perform = True
	while perform:
		line = proc.stdout.readline()
		if line != '':
			data = line.split()
			chrn = data[0]
			start = int(data[1])
			end = int(data[2])
			if isinstance(start, int) and isinstance(end, int):
				mid = (start+end)/2
			strand = str(data[5])
			if mid and strand:
				fo.write("chr%s\t%s\t%s\n"%(chrn, mid, strand))
		else:
			fo.close()
			perform = False


def make_seqpeaklist(base_name, aln_files, aln_folder):
	"""
	Prepare configuration files for seqpeak.
	Base this configuration on the original chip-seq configuration
	
	The seqpeak configuration file contains 2 rows and 2 columns:
	
	"full path of output folder"/ALN_folder/filename_control.aln	1
	"full path of output folder"/ALN_folder/filename_test.aln	0
	
	"""
	
	fn = "%s_filelist.txt"%base_name
	fn_path = "%s%s"%(aln_folder, fn)

	with open(fn_path, 'w') as fo:
		fo.write("%s\t%s\n"%(aln_files[0], 1))
		fo.write("%s\t%s\n"%(aln_files[1], 0))
		






