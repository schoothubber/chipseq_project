import argparse, os

from make_ALN import bam_to_aln_maker, make_seqpeaklist
from do_seqpeak import call_peaks
from get_near_genes import neighbours
from classify_peaks import peak_classifier, peak_distributor

parser = argparse.ArgumentParser()


parser.add_argument(
		"-p", "--pipeline", help="Start the cisgenome pipeline",
		action="store_true")
parser.add_argument(
		"bf", help="full pathname of the BAM files folder"
		)
parser.add_argument(
		"of", help="full pathname of the output folder"
		)
parser.add_argument(
		"af", help="full pathname of the annotation file"
		)


args = parser.parse_args()

if args.pipeline:
	print "Initiate pipeline"
	
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
	
	
	#bam_to_aln_maker(args.bf, aln_folder)
	#make_seqpeaklist(aln_folder)
	#call_peaks(aln_folder, seqpeak_folder)
	#neighbours(args.af, seqpeak_folder, neighbor_folder)
	#peak_classifier(neighbor_folder)
	peak_distributor(neighbor_folder)

else:
	print "lets not do anything thanks to you"
