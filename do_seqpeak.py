#!/usr/bin/env

from subprocess import Popen, PIPE
import os



def call_peaks(bam_test_file, aln_folder, seqpeak_folder):
	"""
	#file to be checked for presence
	#basename_peak.cod	
	"""
	
	#take the base name of the original bam file
	#use it to name the corresponding output files
	base_name_test = os.path.splitext(bam_test_file)[0]
	
	#check if the actual output file is already present in seqpeak_folder
	check_file = "%s%s_peak.cod"%(seqpeak_folder, base_name_test)
	if not os.path.isfile(check_file):
	
		peak_filename = "%s"%base_name_test
		seqpeaklist = "%s/%s_filelist.txt"%(aln_folder, base_name_test)
		
		cmd = [
				"seqpeak", "-i", seqpeaklist, "-d", seqpeak_folder, 
				"-o", peak_filename, "-e", "150", "-maxgap", "200", 
				"-minlen", "200"
				]
		
		print "Running seqpeak..."	
		#pipe = check_output(cmd)
		pipe = Popen(cmd, stdout=PIPE, stderr=PIPE)
		stdout, stderr = pipe.communicate()
	
	else:
		print "The peaks are already called for this ALN file"
		pass


def delete_files(bam_test_file, seqpeak_folder):
	"""	
	#files to be deleted
	#basename.cgw
	#basename_log2fc.bar
	#basename_t.bar
	
	These files are produced by cisgenomes seqpeak tool
	However they are not used and are therefor declared obsolete
	"""
	#get the basename
	base_name_test = os.path.splitext(bam_test_file)[0]
	
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














