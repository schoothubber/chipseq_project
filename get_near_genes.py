from subprocess import Popen, PIPE, check_output
import os
import fnmatch


#test values
#annotation_file = "/home/wouter/chipseq_project/REF_files/hg19/annotation/refLocus_sorted.txt"
#peak_folder = "/home/wouter/chipseq_project/cisgenome2.0/output_data/"
#output_folder = "/home/wouter/chipseq_project/cisgenome2.0/output_data/neighborgenes/"

#neighborgenes folder might have to be created first

def neighbours(annotation_file, seqpeak_folder, neighbor_folder):
	"""
	By utilizing cisgenome's reflocus_getneighborgenes app, a masterfile
	is created that contains information on the identified peaks in relation 
	to the genes on the reference genome.
	This master file is filtered and split up into 4 separate files.
	Each file contains the same format as the master file, but only those
	rows with a certain vicinity of a gene to the TSS.
	"""
	pattern = "*.cod"
	peak_files = [file for file in os.listdir(seqpeak_folder) if fnmatch.fnmatch(file, pattern)]

	neighbour_files = []
	
	for peak_file in peak_files:
		print "peak_File: %s"%peak_file
		peak_in = "%s%s"%(seqpeak_folder, peak_file)
		peak_out = "%s%s%s"%(neighbor_folder, "neighbor_", peak_file)
		neighbour_files.append(peak_out)
		
		print "Running reflocus_getneighborgenes on %s"%peak_file

		#run this command for every peak_file.cod in the peak_folder
		cmd1 = [
				"reflocus_getneighborgenes", "-d" , annotation_file, 
				"-s", "human", "-i", peak_in, "-o", peak_out, "-g", "1000000", 
				"-up", "1000", "-down", "1000"
				]
		#pipe1 = check_output(cmd1)
		pipe1 = Popen(cmd1, stdout=PIPE, stderr=PIPE)
		stdout, stderr = pipe1.communicate()



	
		#delete master file
		#os.remove(neighbour_file)

#neighbours(annotation_file, peak_folder, output_folder)

	#delete the bar and cgw files

	
	
	
	
