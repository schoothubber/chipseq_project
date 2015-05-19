from subprocess import Popen, PIPE
import os
import fnmatch

#test values
bam_folder = "/home/wouter/chipseq_project/BAM_files/"
bed_folder = "/home/wouter/chipseq_project/BED_files/"


def bam_to_aln_maker(bam_folder, bed_folder):
	"""
	Take all bam files from a bam folder and convert them to bed files

	"""

	#Take all files from the bamfolder that have the .bam extension
	ext = "*.bam"
	bam_files = [file for file in os.listdir(bam_folder) if fnmatch.fnmatch(file, ext)]
	#for each bam file perform bedtools and extract aln info
	for bam_file in bam_files:
		print bam_file
		
		#set the Bedtools command line
		bam_in = "%s%s"%(bam_folder, bam_file)
		
		base_name = os.path.splitext(bam_file)[0]
		bed_out = "%s%s%s"%(bed_folder, base_name, '.bed')
		
		pipe_cmd = "bedtools bamtobed -i %s > %s"%(bam_in,bed_out)
		pipe1 = Popen(pipe_cmd, stdout=PIPE, stderr=PIPE, shell=True)
		stdout, stderr = pipe1.communicate()





bam_to_aln_maker(bam_folder, bed_folder)
