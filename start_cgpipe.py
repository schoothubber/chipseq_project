
import time
import argparse, os

from prepare_data import get_file_names, remove_duplicates, converter
from do_seqpeak import call_peaks, delete_files
from get_near_genes import neighbours
from classify_peaks import (
				peak_classifier, peak_distributor,
				peak_localizer, peaks_per_chromosome
							)

parser = argparse.ArgumentParser()

parser.add_argument(
		"-p", "--cgpipeline", help="Start the cisgenome pipeline",
		action="store_true")
parser.add_argument(
		"cf", help="full pathname of the Configuration file folder"
		)
parser.add_argument(
		"bf", help="full pathname of the BAM files folder"
		)
parser.add_argument(
		"of", help="full pathname of the output folder"
		)
parser.add_argument(
		"af", help="full pathname of the annotation file"
		)
	

#command:
#python start.py -p <configfile.ini> <bam file folder> <output folder> <refLocus_sorted.txt>

#python start.py -p /home/wouter/chipseq_project/cisgenome2.0/Config/chipseq_configuration.ini /home/wouter/chipseq_project/BAM_files/ /home/wouter/chipseq_project/TESToutput/ /home/wouter/chipseq_project/REF_files/hg19/annotation/refLocus_sorted.txt


args = parser.parse_args()

if args.cgpipeline:
	
	print "Initiate the cisgenome pipeline"
	
	#create main output folder
	if not os.path.exists(args.of):
		os.mkdir(args.of)

	#create output sub-folders	
	aln_folder = "%sALN_files/"%args.of
	if not os.path.exists(aln_folder):
		os.mkdir(aln_folder)
		
	seqpeak_folder = "%sseqpeak_output/"%args.of
	if not os.path.exists(seqpeak_folder):
		os.mkdir(seqpeak_folder)
		
	neighbor_folder = "%sneighborgenes/"%args.of
	if not os.path.exists(neighbor_folder):
		os.mkdir(neighbor_folder)
	
	#get the filenames of the bam files which are written in a config file
	bam_files_data = get_file_names(args.cf)
	
	
	for bam_test_file, bam_control_file in bam_files_data:
		
		tic = time.clock()
		
		remove_duplicates(bam_test_file, args.bf)
		converter(bam_test_file, bam_control_file, args.bf, aln_folder)
		call_peaks(bam_test_file, aln_folder, seqpeak_folder)
		delete_files(bam_test_file, seqpeak_folder)
		neighbours(args.af, seqpeak_folder, neighbor_folder, bam_test_file)
		peak_classifier(neighbor_folder, bam_test_file)
		peak_distributor(neighbor_folder, bam_test_file)
		peak_localizer(neighbor_folder, bam_test_file)
		peaks_per_chromosome(seqpeak_folder, neighbor_folder, bam_test_file)
		
		toc = time.clock()
		print toc - tic

else:
	print "lets not do anything thanks to you"
