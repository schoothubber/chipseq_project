from subprocess import Popen, PIPE
import os
import fnmatch



def call_peaks(aln_folder, seqpeak_folder):
	"""
	
	"""
	"/home/wouter/chipseq_project/ALN_files/seqpeak_filelist.txt"
	"/home/wouter/chipseq_project/cisgenome2.0/output_data/"



	
	peak_filename = "test_peak"
	seqpeaklist = "%s/seqpeak_filelist.txt"%aln_folder
	
	cmd1 = [
			"seqpeak", "-i", seqpeaklist, "-d", seqpeak_folder, "-o", 
			peak_filename, "-e", "150", "-maxgap", "200", "-minlen", "200"
			]
	print "Running seqpeak..."	
	#pipe1 = check_output(cmd1)
	pipe1 = Popen(cmd1, stdout=PIPE, stderr=PIPE)
	stdout, stderr = pipe1.communicate()
