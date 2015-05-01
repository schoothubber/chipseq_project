from subprocess import Popen, PIPE
import os
import fnmatch

#test values
#bam_folder = "/home/wouter/chipseq_project/BAM_files/"
#aln_folder = "/home/wouter/chipseq_project/ALN_files/"


def bam_to_aln_maker(bam_folder, aln_folder):
	"""
	Take all bam files from a bam folder and convert them to aln files

	"""

	#Take all files from the bamfolder that have the .bam extension
	pattern = "*.bam"
	bam_files = [file for file in os.listdir(bam_folder) if fnmatch.fnmatch(file, pattern)]
	#for each bam file perform bedtools and extract aln info
	for bam_file in bam_files:
		print bam_file
		
		#set the Bedtools command line
		bam_in = "%s%s"%(bam_folder, bam_file)
		pipe_cmd = ["bedtools", "bamtobed", "-i", bam_in]
		
		#take everything except the .bam extension
		base_name = os.path.splitext(bam_file)[0]
		#set aln output name
		aln_out = "%s%s.aln"%(aln_folder, base_name)
		
		#create a separate aln file for each process 
		fo = open(aln_out, 'w')
		print "created %s"%aln_out

		#run bedtools on each bamfile
		#and extract data for aln files
		proc = Popen(pipe_cmd, stdout=PIPE)
		perform = True
		while perform:
			line = proc.stdout.readline()
			if line != '':
				data = line.split()
				chrn = data[0]
				start = int(data[1])
				end = int(data[2])
				if isinstance(start, int) and isinstance(end, int):
					mid = (start+end)/2
				strand = str(data[5])
				if mid and strand:
					fo.write("chr%s\t%s\t%s\n"%(chrn, mid, strand))
			else:
				fo.close()
				perform = False


def make_seqpeaklist(aln_folder):
	"""
	Take all ALN files and prepare a configuration text file
	containing the names of the ALN files and
	"""
	#Take all files from the alnfolder that have the .aln extension
	pattern = "*.aln"
	aln_files = [aln_folder + file for file in os.listdir(aln_folder) if fnmatch.fnmatch(file, pattern)]

	fn = "seqpeak_filelist.txt"
	fn_path = "%s%s"%(aln_folder, fn)
	with open(fn_path, 'w') as fo:
		for aln_file in aln_files:
			if "Input" in aln_file:
				sample = 0
			else:
				sample = 1
			fo.write("%s\t%s\n"%(aln_file, sample))






