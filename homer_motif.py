
######################################
###Author: W. van der Schoot
######################################

from subprocess import Popen, PIPE, check_output
import os

from prepare_data import get_base_name




def HOMER_make_motifs(fileargs, motifargs):
	"""
	findMotifsGenome.pl 	/home/wouter/chipseq_project/homer/findpeak_output/peaks2.txt 
							hg19 
							/home/wouter/chipseq_project/homer/motif_output2/
							-size given
							-S
							-mis
							-len
	"""
	#folder parameters
	bam_test_file = fileargs['bam_test_file']
	seqpeak_folder = fileargs['seqpeak_folder']
	motif_folder = fileargs['motif_folder']
	
	#tool parameters
	genome = motifargs["mgenom"]
	motif_size = motifargs["msize"]
	motif_len = motifargs["mlen"]
	motif_num = motifargs["mnum"]
	mismatches = motifargs["mism"]
	
	#Prepare the BED file names for input
	#take everything except the .bam extension
	base_name_test = get_base_name(bam_test_file)

	#Make a name for the rewritten peak file
	peak_file = "%s%s_peak.cod"%(seqpeak_folder, base_name_test)
	
	#Make a name for the motif folder
	motif_output = "%s%s_motifoutput/"%(motif_folder, base_name_test)
	
	#cmd = ['findMotifsGenome.pl', peak_file, 'hg38', motif_output, '-size', 'given']
	#cmd = 'findMotifsGenome.pl %s hg38 %s -size given'%(peak_file, motif_output)
	
	cmd = 'findMotifsGenome.pl %s %s %s -size %s -S %s -mis %s -len %s'%(
							peak_file, genome, motif_output, motif_size,
							motif_num, mismatches, motif_len
							)
	
	print "Start finding motifs for %s"%base_name_test
	check_output(cmd, shell=True)
	#pipe = Popen(cmd, stdout=PIPE, stderr=PIPE)
	#stdout, stderr = pipe.communicate()








