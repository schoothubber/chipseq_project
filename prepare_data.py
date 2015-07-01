#!/usr/bin/env

from subprocess import Popen, PIPE
import os
import random as rd



def converter(fileargs):
	"""
	Pass a bam file to the actual process of converting
	Check for each bam file to be processed if a corresponding
	aln file is already present in the output file.
	"""
	
	bam_test_file = fileargs['bam_test_file']
	bam_ctrl_file = fileargs['bam_ctrl_file']
	aln_folder = fileargs['aln_folder']
	
	tofile = [bam_test_file, bam_ctrl_file]
	aln_files = []
	
	for fbam in tofile:
	    base_name = get_base_name(fbam)
	    aln_out =  "%s%s%s"%(aln_folder, base_name, ".aln")
	    aln_files.append(aln_out)
	    
	    if not os.path.isfile(aln_out):
		print "start BAM to ALN conversion for %s"%aln_out
		make_aln_from_bam(fbam, aln_out)
	    else:
		stat_size = os.stat(aln_out)[6]
		if stat_size == 0:
		    make_aln_from_bam(fbam, aln_out)
		else:
		    print "%s is already present, skipping conversion..."%aln_out
		    pass
	
	base_name = get_base_name(bam_test_file)
	make_seqpeaklist(tofile, aln_files, aln_folder)
	

def get_base_name(fullpath):
	"""
	Take a full path and extract the filename  without the extension
	"""
	filename = os.path.split(fullpath)[1]
	base_name = os.path.splitext(filename)[0]
	
	return base_name


		
		
			
def make_aln_from_bam(bam_in, aln_out):
	"""
	The name says it all
	"""

	cmd = ["sambamba", "view", "-t", "4", "-F", "not duplicate", bam_in]
	fo = open(aln_out, 'w')
	proc = Popen(cmd, stdout=PIPE)
	perform = True
	while perform:
		line = proc.stdout.readline()
		if line != '':
			data = line.split()
			
			strand_flag = int(data[1])
			if strand_flag == 16:
				strand = "-"
			else:
				strand = "+"
			
			chrn = data[2]
			
			start = int(data[3])
			seq = data[9]
			end = start + len(seq)
			if isinstance(start, int) and isinstance(end, int):
				mid = (start+end)/2
			
			if mid and strand:
				fo.write("chr%s\t%s\t%s\n"%(chrn, mid, strand))
		else:
			fo.close()
			perform = False
	




def make_seqpeaklist(tofile, aln_files, aln_folder):
	"""
	Prepare configuration files for seqpeak.
	Base this configuration on the original chip-seq configuration
	
	The seqpeak configuration file contains 2 rows and 2 columns:
	
	"full path of output folder"/ALN_folder/filename_control.aln	1
	"full path of output folder"/ALN_folder/filename_test.aln	0
	
	"""
	base_name = get_base_name(tofile[0])
	fn = "%s_filelist.txt"%base_name
	fn_path = "%s%s"%(aln_folder, fn)

	with open(fn_path, 'w') as fo:
		fo.write("%s\t%s\n"%(aln_files[0], 1))
		fo.write("%s\t%s\n"%(aln_files[1], 0))
		




def make_random_peak_list(random_folder, aln_folder, base_name_ctrl, N):
	"""
	N = number of peaks in the test file
	the amount of random peaks should be 10N < 100000
	
	A seed is produced as a random integer of 8 numbers
	This seed is stored in a text file so that the random
	peaks list is reproducible.
	"""
	
	
	
	aln_path = "%s%s%s"%(aln_folder, base_name_ctrl, '.aln')
	with open(aln_path, 'r') as aln_fo:
	    aln_data = []
	    for line in aln_fo:
		info = line.split()
		
		try:
			aln_data.append([info[0], int(info[1]), info[2]])
		
		except ValueError:
			pass
		
	nr_rand_peaks = N * 10
	if nr_rand_peaks > len(aln_data):
		nr_rand_peaks = len(aln_data)
	elif nr_rand_peaks > 100000:
		nr_rand_peaks = 100000
	else:
		pass
	
	seedling = rd.randint(11111111, 99999999)
	fn_seed = "%s%s%s"%('seed_' ,str(seedling), '.txt')
	seed_path = "%s%s"%(random_folder, fn_seed)
	with open(seed_path, 'w') as sd_fo:
	    sd_fo.write("seed: %s"%seedling)
	
	rd.seed(seedling)
	rd.shuffle(aln_data)
	random_peaks = aln_data[0:nr_rand_peaks]
	
	#write the random peaks to a file
	random_path = "%s%s"%(random_folder, "random_data_peak.cod")
	rank = 1
	with open(random_path, 'w') as rnd_fo:
		rnd_fo.write("%s\t%s\t%s\t%s\t%s\n"%("rank", "chromosome", "start", "end", "strand"))
		for c, p, s in random_peaks:
			rnd_fo.write("%s\t%s\t%s\t%s\t%s\n"%(rank, c, int(p)-500, int(p)+500, s))
			rank += 1





















