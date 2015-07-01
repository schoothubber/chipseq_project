
######################################
###Author: W. van der Schoot
######################################

from subprocess import Popen, PIPE
import os

from prepare_data import get_base_name


def neighbours(fileargs, annoargs):
	"""
	By utilizing cisgenome's reflocus_getneighborgenes app, a masterfile
	is created that contains information on the identified peaks in 
	relation to the genes on the reference genome.
	This master file is filtered and split up into 4 separate files.
	Each file contains the same format as the master file, but only 
	those rows with a certain vicinity of a gene to the TSS.
	
	
	"reflocus_getneighborgenes", "-d" , annotation_file, 
	"-s", species, "-i", peak_in, "-o", peak_out, "-g", 
	g_distance, "-up", upstream_genes, "-down", downstream_genes
	
	"""
	#folder parameters
	bam_test_file = fileargs['bam_test_file']
	seqpeak_folder = fileargs['seqpeak_folder']
	neighbor_folder = fileargs['neighbor_folder']
	
	#tool parameters
	species = annoargs["nspec"]
	annotation_file = annoargs["nanno"]
	g_distance = annoargs["ngdist"]
	upstream_genes = annoargs["nusgenes"]
	downstream_genes = annoargs["ndsgenes"]
	
	#process variables
	base_name_test = get_base_name(bam_test_file)
	peak_file = "%s_peak.cod"%base_name_test
	peak_in = "%s%s"%(seqpeak_folder, peak_file)
	peak_out = "%s%s%s"%(neighbor_folder, "neighbour_", peak_file)
	
	print "Running reflocus_getneighborgenes on %s"%peak_file
	
	cmd = [
			"reflocus_getneighborgenes", "-d" , annotation_file, 
			"-s", species, "-i", peak_in, "-o", peak_out, "-g", 
			g_distance, "-up", upstream_genes, "-down", downstream_genes
			]
			
	#pipe = check_output(cmd)
	pipe = Popen(cmd, stdout=PIPE, stderr=PIPE)
	stdout, stderr = pipe.communicate()
	
	
	
	
