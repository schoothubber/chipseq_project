
import time
import argparse, os


from prepare_data import get_file_names, remove_duplicates
from homer_motif import (
				BED_make_beds, HOMER_make_tags, 
				HOMER_make_peaks, PYTHON_rewrite_peakfile, 
				HOMER_make_motifs
						)

def check_the_time(tic):
	toc = time.clock()
	print toc - tic


parser = argparse.ArgumentParser()

parser.add_argument(
		"-m", "--homermotif", help="Start HOMERs motif analysis",
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

args = parser.parse_args()
		
if args.homermotif:
	
	print "Initiate the Homer de novo Motif Analysis pipeline"

	#create main output folder
	if not os.path.exists(args.of):
		os.mkdir(args.of)
		
	#create output sub-folders	
	bed_folder = "%sBED_files/"%args.of
	if not os.path.exists(bed_folder):
		os.mkdir(bed_folder)
		
	tag_folder = "%sTAG_files/"%args.of
	if not os.path.exists(tag_folder):
		os.mkdir(tag_folder)
		
	peak_folder = "%sPEAK_files/"%args.of
	if not os.path.exists(peak_folder):
		os.mkdir(peak_folder)

	motif_folder = "%sMOTIF_files/"%args.of
	if not os.path.exists(motif_folder):
		os.mkdir(motif_folder)

	#get the filenames of the bam files which are written in a config file
	bam_files_data = get_file_names(args.cf)
	
	for bam_test_file, bam_control_file in bam_files_data:
		
		tic = time.clock()
		
		remove_duplicates(bam_test_file, args.bf)
		check_the_time(tic)
		
		BED_make_beds(bam_test_file, bam_control_file, args.bf, bed_folder)
		check_the_time(tic)
		
		HOMER_make_tags(bam_test_file, bam_control_file, tag_folder, bed_folder)
		check_the_time(tic)
		
		HOMER_make_peaks(bam_test_file, bam_control_file, peak_folder, tag_folder)
		check_the_time(tic)
		
		PYTHON_rewrite_peakfile(bam_test_file, peak_folder)
		check_the_time(tic)
		
		HOMER_make_motifs(bam_test_file, peak_folder, motif_folder)
		check_the_time(tic)
		

else:
	print "lets not do anything thanks to you"





"""
MOTIF COMMAND LINE USAGE:

python start_motifpipe.py -m /home/wouter/chipseq_project/cisgenome2.0/Config/chipseq_configuration.ini /home/wouter/chipseq_project/BAM_files/ /home/wouter/chipseq_project/denovomotif_output/

"""




