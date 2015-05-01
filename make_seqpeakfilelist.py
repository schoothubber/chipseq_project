import os
import fnmatch

#create a tab-delimited text file that contains all DNA samples to be analyzed.
#The first column gives the full file path. 
#The second column indicates the DNA sample type (1 = IP, 0 = control).

aln_folder = "/home/wouter/chipseq_project/ALN_files/"

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
