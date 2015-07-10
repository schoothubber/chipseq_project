#!/usr/bin/env python2.7

######################################
#######Author: W. van der Schoot######
######################################

import argparse, os, sys

from logging import parameter_logger
from prepare_data import converter
from do_seqpeak import call_peaks, delete_files
from get_near_genes import neighbours
from classify_peaks import (peak_classifier, peak_distributor,
				peak_localizer, peaks_per_chromosome)
from homer_motif import make_motifs

parser = argparse.ArgumentParser()

parser.add_argument(
		"-p", "--cgpipeline", help="Start the cisgenome pipeline",
		action="store_true")

#tool path names
parser.add_argument("sam", help="Full pathname to sambamba")
#files and output folder
parser.add_argument("bt", help="Full pathname BAM test file")
parser.add_argument("bc", help="Full pathname BAM control file")
parser.add_argument("of", help="Full pathname of the output folder")
#seqpeak parameters
parser.add_argument("srext", help="reads with be extended E base pairs to 3 ends. The extended DNA fragments will be used to determine genome coverage.")
parser.add_argument("sbsize", help="the genome will be divided into non-overlapping bins, each B bp long. For each bin, the program will count how many DNA fragments cover the bin. Later, the counts will be normalized across samples, taking into account of sequencing depths, and the normalized counts will be used to call peaks.")
parser.add_argument("shwsize", help="To determine enrichment signal, each bin and its 2*W flanking bins (W on the left and W on the right) will be combined into a window. The normalized DNA fragment counts in these 2*W+1 bins will be compared between IP and control samples. A t-statistic similar to the TileMap t-statistic (formula 6 in the paper) will be computed. The t-statistic will be used as the window statistic for the center bin. For transcription factors (narrow peaks), the value of W is usually chosen such that (2W+1)*B = 150 ~ 500. For histone modifications (broad peaks), the value of W is usually chosen such that (2W+1)*B = 500 ~1000.")
parser.add_argument("sstan", help="If checked, the window statistics across all bins will be standardized by subtracting genome-wide mean and dividing by genome-wide standard deviation.")
parser.add_argument("smpgap", help="If the distance between two peaks are smaller than Max_Gap bp, then they will be merged into one peak.")
parser.add_argument("smrlen", help="If the peak length is smaller than Min_Peak_Length bp, then the peak will not be reported.")
parser.add_argument("sbound", help="After peaks are detected, we use a finer bin (5 bp by default or Y bp specified by Boundary Refinement Resolution) to analyze forward and reverse strand reads separately. The offset between the forward and reverse strand peaks will be used to determine the refined peak boundaries. This is very useful for narrowing down the transcription factor binding sites to a resolution of <50 bp.")
#neighbor parameters
parser.add_argument("nspec", help="current support: human and mouse.")
parser.add_argument("nanno", help="The full path of he location to the annotation reference file.")
parser.add_argument("ngdist", help="distance upper limit")
parser.add_argument("nusgenes", help="no. of upstream genes")
parser.add_argument("ndsgenes", help="no. of downstream genes")
#motif parameters
parser.add_argument("mgenom", help="current support: hg38. mm10")
parser.add_argument("msize", help="The size of the region used for motif finding is important.  If analyzing ChIP-Seq peaks from a transcription factor, Chuck would recommend 50 bp for establishing the primary motif bound by a given transcription factor and 200 bp for finding both primary and co-enriched motifs for a transcription factor.  When looking at histone marked regions, 500-1000 bp is probably a good idea (i.e. H3K4me or H3/H4 acetylated regions).  In theory, HOMER can work with very large regions (i.e. 10kb), but with the larger the regions comes more sequence and longer execution time.  These regions will be based off the center of the peaks.  If you prefer an offset, you can specify -size -300,100 to search a region of size 400 that is centered 100 bp upstream of the peak center (useful if doing motif finding on putative TSS regions).  If you have variable length regions, use the option -size given and HOMER will use the exact regions that were used as input.")
parser.add_argument("mlen", help="motif length, default=8,10,12) [NOTE: values greater 12 may cause the programto run out of memory - in these cases decrease the number of sequences analyzed (-N),or try analyzing shorter sequence regions (i.e. -size 100)]")
parser.add_argument("mnum", help="Number of motifs to optimize, default: 25")
parser.add_argument("mism", help="global optimization: searches for strings with # mismatches, default: 2")


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

	log_folder = "%s/logs/"%main_output
	if not os.path.exists(log_folder):
		os.mkdir(log_folder)
	
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
		'nanno' : args.nanno,
		'log_folder' : log_folder,
		'sambamba_path' : args.sam
		}
	
	#seqpeak parameters
	peakargs = {
		"seqpeak_readextension" : args.srext,
		"seqpeak_binsize" : args.sbsize,
		"seqpeak_halfwinsize" : args.shwsize,
		"seqpeak_standardize" : args.sstan,
		"seqpeak_maxpeakgap" : args.smpgap,
		"seqpeak_minreglen" : args.smrlen,
		"seqpeak_boundary" : args.sbound
		}
	
	#neighbor parameters
	annoargs = {
		"neighbor_species" : args.nspec,
		"neighbor_annotationfile" : args.nanno,
		"neighbor_gdistance" : args.ngdist,
		"neighbor_upstreamgenes" : args.nusgenes,
		"neighbor_downstreamgenes" : args.ndsgenes
		}
	
	#motif parameters
	motifargs = {
		"motif_genome" : args.mgenom,
		"motif_size" : args.msize,
		"motif_length" : args.mlen,
		"motif_number" : args.mnum,
		"motif_mismatch" : args.mism
		}
	
	parameter_logger(fileargs, peakargs, annoargs, motifargs)
	converter(fileargs)
	call_peaks(fileargs, peakargs)
	delete_files(fileargs)
	neighbours(fileargs, annoargs)
	peak_classifier(fileargs)
	peak_distributor(fileargs, peakargs, annoargs)
	peak_localizer(fileargs)
	peaks_per_chromosome(fileargs)
	make_motifs(fileargs, motifargs)
	

else:
	print "Something went wrong, better check the input parameters!"
