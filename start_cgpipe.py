
######################################
###Author: W. van der Schoot
######################################


import time
import argparse, os

from prepare_data import converter
from do_seqpeak import call_peaks, delete_files
from get_near_genes import neighbours
from classify_peaks import (
				peak_classifier, peak_distributor,
				peak_localizer, peaks_per_chromosome
						)
from homer_motif import HOMER_make_motifs

parser = argparse.ArgumentParser()

parser.add_argument(
		"-p", "--cgpipeline", help="Start the cisgenome pipeline",
		action="store_true")

#files and output folder
parser.add_argument("bt", help="Full pathname BAM test file")
parser.add_argument("bc", help="Full pathname BAM control file")
parser.add_argument("of", help="Full pathname of the output folder")
#seqpeak parameters
parser.add_argument("srext", help="")
parser.add_argument("sbsize", help="")
parser.add_argument("shwsize", help="")
parser.add_argument("sstan", help="")
parser.add_argument("smpgap", help="")
parser.add_argument("smrlen", help="")
parser.add_argument("sbound", help="")
#neighbor parameters
parser.add_argument("nspec", help="")
parser.add_argument("nanno", help="")
parser.add_argument("ngdist", help="")
parser.add_argument("nusgenes", help="")
parser.add_argument("ndsgenes", help="")
#motif paramters
parser.add_argument("mgenom", help="")
parser.add_argument("msize", help="")
parser.add_argument("mlen", help="")
parser.add_argument("mnum", help="")
parser.add_argument("mism", help="")


args = parser.parse_args()


if args.cgpipeline:
	
	print "Initiate the chip-seq analysis pipeline"
	
	main_output = args.of
	
	#create main output folder
	if not os.path.exists(main_output):
		os.mkdir(main_output)

	#create output sub-folders	
	aln_folder = "%s/ALN_files/"%main_output
	if not os.path.exists(aln_folder):
		os.mkdir(aln_folder)

	seqpeak_folder = "%s/seqpeak_output/"%main_output
	if not os.path.exists(seqpeak_folder):
		os.mkdir(seqpeak_folder)
		
	neighbor_folder = "%s/neighborgenes/"%main_output
	if not os.path.exists(neighbor_folder):
		os.mkdir(neighbor_folder)
		
	motif_folder = "%s/motif/"%main_output
	if not os.path.exists(motif_folder):
		os.mkdir(motif_folder)
		
	random_folder = "%s/random/"%main_output
	if not os.path.exists(random_folder):
		os.mkdir(random_folder)
	
	#some hash:
	
	#folder parameters
	fileargs = {
		'main_output' : main_output,
		'bam_test_file' : args.bt,
		'bam_ctrl_file' : args.bc,
		'aln_folder' : aln_folder,
		'seqpeak_folder' : seqpeak_folder,
		'neighbor_folder' : neighbor_folder,
		'motif_folder' : motif_folder,
		'random_folder' : random_folder,
		'nanno' : args.nanno
		}
	
	#seqpeak parameters
	peakargs = {
		"srext" : args.srext,
		"sbsize" : args.sbsize,
		"shwsize" : args.shwsize,
		"sstan" : args.sstan,
		"smpgap" : args.smpgap,
		"smrlen" : args.smrlen,
		"sbound" : args.sbound
		}
	
	#neighbor parameters
	annoargs = {
		"nspec" : args.nspec,
		"nanno" : args.nanno,
		"ngdist" : args.ngdist,
		"nusgenes" : args.nusgenes,
		"ndsgenes" : args.ndsgenes
		}
	
	#motif paramters
	motifargs = {
		"mgenom" : args.mgenom,
		"msize" : args.msize,
		"mlen" : args.mlen,
		"mnum" : args.mnum,
		"mism" : args.mism
		}
	
	converter(fileargs)
	call_peaks(fileargs, peakargs)
	delete_files(fileargs)
	neighbours(fileargs, annoargs)
	peak_classifier(fileargs)
	peak_distributor(fileargs, annoargs)
	peak_localizer(fileargs)
	peaks_per_chromosome(fileargs)
	HOMER_make_motifs(fileargs, motifargs)
	

else:
	print "Something went wrong, better check the input parameters!"
