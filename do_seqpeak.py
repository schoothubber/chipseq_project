#!/usr/bin/env python2.7

######################################
###Author: W. van der Schoot
######################################

from subprocess import Popen, PIPE
import os

from prepare_data import get_base_name


def call_peaks(fileargs, peakargs):
	"""
	#file to be checked for presence
	#basename_peak.cod	
	
	The peak calling is started by the following command:
	
	"seqpeak", "-i", seqpeaklist, "-d", seqpeak_folder, 
	"-o", peak_filename, "-e", read_ext_len, "-b", bin_size, 
	"-w", half_win_size, "ts", stand_win_stat, 
	"-maxgap", max_peak_gap, "-minlen", min_reg_len, 
	"-br", refine_boundary
	
	Some values are actually a yes or a no but the input is a 1 or a 0,
	respectively
	"""
	#Tool Path
	seqpeak = "%s/bin/seqpeak"%fileargs['cisgenome_path']
	#folder parameters
	bam_test_file = fileargs['bam_test_file']
	aln_folder = fileargs['aln_folder']
	seqpeak_folder = fileargs['seqpeak_folder']

	#tool parameters
	read_ext_len = peakargs['seqpeak_readextension']
	bin_size = peakargs['seqpeak_binsize']
	half_win_size = peakargs['seqpeak_halfwinsize']
	stand_win_stat = peakargs['seqpeak_standardize']
	max_peak_gap = peakargs['seqpeak_maxpeakgap']
	min_reg_len = peakargs['seqpeak_minreglen']
	refine_boundary = peakargs['seqpeak_boundary']
	
	#take the base name of the original bam file
	#use it to name the corresponding output files
	base_name_test = get_base_name(bam_test_file)
	
	#check if the actual output file is already present in seqpeak_folder
	check_file = "%s%s_peak.cod"%(seqpeak_folder, base_name_test)
	if not os.path.isfile(check_file):
	
		peak_filename = "%s"%base_name_test
		seqpeaklist = "%s/%s_filelist.txt"%(aln_folder, base_name_test)
		
		cmd = [
			seqpeak, "-i", seqpeaklist, "-d", seqpeak_folder, 
			"-o", peak_filename, "-e", read_ext_len, "-b", bin_size, 
			"-w", half_win_size, "-ts", stand_win_stat, 
			"-maxgap", max_peak_gap, "-minlen", min_reg_len, 
			"-br", refine_boundary
			]
		
#		cmd = "seqpeak -i %s -d %s -o %s -e %s -b %s -w %s -ts %s -maxgap %s -minlen %s -br %s"%
#				(
#				seqpeaklist,seqpeak_folder,peak_filename, 
#				read_ext_len, bin_size, half_win_size, 
#				stand_win_stat, max_peak_gap, min_reg_len, 
#				refine_boundary
#				)
		
		print "Running seqpeak..."	
#		pipe = check_output(cmd, shell=True)
		pipe = Popen(cmd, stdout=PIPE, stderr=PIPE)
		stdout, stderr = pipe.communicate()
	
	else:
		print "The peaks are already called for this ALN file"
		pass


def delete_files(fileargs):
	"""	
	#files to be deleted
	#basename.cgw
	#basename_log2fc.bar
	#basename_t.bar
	
	These files are produced by cisgenomes seqpeak tool
	However they are not used and are therefor declared obsolete
	"""
	#folder parameters
	bam_test_file = fileargs['bam_test_file']
	seqpeak_folder = fileargs['seqpeak_folder']
	
	#get the basename
	base_name_test = get_base_name(bam_test_file)
	
	#make filenames to-be-deleted
	file_cgw = "%s%s.cgw"%(seqpeak_folder, base_name_test)
	file_2fcbar = "%s%s_log2fc.bar"%(seqpeak_folder, base_name_test)
	file_tbar = "%s%s_t.bar"%(seqpeak_folder, base_name_test)
	
	obsolete_files = [file_cgw, file_2fcbar, file_tbar]
	
	#delete the obsolete files
	print "Removing obsolete files..."
	for obsolete_file in obsolete_files:
		if os.path.isfile(obsolete_file):
			os.remove(obsolete_file)
			print "removed %s"%obsolete_file
		else:
			pass














