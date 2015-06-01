#!/usr/bin/env

from subprocess import Popen, PIPE
import os



def BED_make_beds(bam_test_file, bam_control_file, bam_folder, bed_folder):
	"""
	"""
	
	
	#make the filenames including the full path
	bam_test = "%s%s%s"%(bam_folder, "dedupped_", bam_test_file)
	bam_ctrl = "%s%s"%(bam_folder, bam_control_file)
	
	#take everything except the .bam extension
	base_name_test = os.path.splitext(bam_test_file)[0]
	base_name_ctrl = os.path.splitext(bam_control_file)[0]
	
	#set bed output name
	bed_out_test = "%s%s.bed"%(bed_folder, base_name_test)
	bed_out_ctrl = "%s%s.bed"%(bed_folder, base_name_ctrl)
	
	
	#check if this bed TEST file is already present in the bed folder
	if not os.path.isfile(bed_out_test):
		print "start BAM to BED conversion for %s"%bed_out_test
		cmd_test = "bedtools bamtobed -i %s > %s"%(bam_test,bed_out_test)
		pipe_t = Popen(cmd_test, stdout=PIPE, stderr=PIPE, shell=True)
		stdout, stderr = pipe_t.communicate()
		
	else:
		print "%s is already present, skipping conversion..."%bed_out_test
		pass
		
	#check if this bed CONTROLfile is already present in the bed folder
	if not os.path.isfile(bed_out_ctrl):
		print "start BAM to BED conversion for %s"%bed_out_ctrl
		cmd_ctrl = "bedtools bamtobed -i %s > %s"%(bam_ctrl,bed_out_ctrl)
		pipe_c = Popen(cmd_ctrl, stdout=PIPE, stderr=PIPE, shell=True)
		stdout, stderr = pipe_c.communicate()
		
	else:
		print "%s is already present, skipping conversion..."%bed_out_ctrl
		pass
		



def HOMER_make_tags(bam_test_file, bam_control_file, tag_folder, bed_folder):
	"""
	makeTagDirectory 	/home/wouter/chipseq_project/homer/tagfolder/control/ 
						-format bed 
						/home/wouter/chipseq_project/BED_files/control/L01_Input1_F3_20141107.bed 


	"""
	
	
	#Prepare the BED file names for input
	#take everything except the .bam extension
	base_name_test = os.path.splitext(bam_test_file)[0]
	base_name_ctrl = os.path.splitext(bam_control_file)[0]
	
	#set bed input name
	bed_in_test = "%s%s.bed"%(bed_folder, base_name_test)
	bed_in_ctrl = "%s%s.bed"%(bed_folder, base_name_ctrl)

	#make the tag folders
	#tag test folder
	tag_test_folder = "%s%s/"%(tag_folder,base_name_test)
	#tag control folder
	tag_ctrl_folder = "%s%s/"%(tag_folder,base_name_ctrl)
	
	#actually make the folders if they dont exist yet
	if not os.path.exists(tag_test_folder):
		os.mkdir(tag_test_folder)
	if not os.path.exists(tag_ctrl_folder):
		os.mkdir(tag_ctrl_folder)


	#check if this bed TEST file is present in the bed folder
	if os.path.isfile(bed_in_test):
		print "Start tagging %s"%tag_test_folder
		cmd_test = ['makeTagDirectory', tag_test_folder, '-format', 'bed', bed_in_test]
		pipe_t = Popen(cmd_test, stdout=PIPE, stderr=PIPE)
		stdout, stderr = pipe_t.communicate()
		
	else:
		print "%s is already present, skipping..."%tag_test_folder
		pass
		
	#check if this bed CONTROL file is present in the bed folder
	if os.path.isfile(bed_in_ctrl):
		print "Start tagging %s"%tag_ctrl_folder
		cmd_ctrl = ['makeTagDirectory', tag_ctrl_folder, '-format', 'bed', bed_in_ctrl]
		pipe_c = Popen(cmd_ctrl, stdout=PIPE, stderr=PIPE)
		stdout, stderr = pipe_c.communicate()
		
	else:
		pass
	


def HOMER_make_peaks(bam_test_file, bam_control_file, peak_folder, tag_folder):
	"""
	findPeaks 			/home/wouter/chipseq_project/homer/tagfolder/ 
						-style factor 
						-o /home/wouter/chipseq_project/homer/findpeak_output/peaks.txt 
						-i /home/wouter/chipseq_project/homer/tagfolder/control/
	"""

	#Prepare the BED file names for input
	#take everything except the .bam extension
	base_name_test = os.path.splitext(bam_test_file)[0]
	base_name_ctrl = os.path.splitext(bam_control_file)[0]

	#make the tag folders
	#tag test folder
	tag_test_folder = "%s%s/"%(tag_folder,base_name_test)
	#tag control folder
	tag_ctrl_folder = "%s%s/"%(tag_folder,base_name_ctrl)
	
	#Make a name for the peak file
	output_peak = "%s%s_peaks.txt"%(peak_folder, base_name_test)
	
	cmd = [
			'findPeaks', tag_test_folder, '-style', 'factor', 
			'-o', output_peak, '-i', tag_ctrl_folder
			]

	if not os.path.isfile(output_peak):
		print "Start peaking %s"%output_peak
		#pipe = check_output(cmd1)
		pipe = Popen(cmd, stdout=PIPE, stderr=PIPE)
		stdout, stderr = pipe.communicate()
	
	else:
		print "%s is already present, skipping..."%output_peak



def PYTHON_rewrite_peakfile(bam_test_file, peak_folder):
	"""
	The plan is to rewrite HOMERs peakfile.txt so that the chromosome
	numbers will have the format 'chr#' where # is an integer or X or Y
	
	"""
	Homer_Peak_File1 = "/home/wouter/chipseq_project/homer/findpeak_output/peaks.txt"
	Homer_Peak_File2 = "/home/wouter/chipseq_project/homer/findpeak_output/peaks3.txt"
	
	#Prepare the BED file names for input
	#take everything except the .bam extension
	base_name_test = os.path.splitext(bam_test_file)[0]
	
	#Make a name for the peak file to be rewritten
	output_peak = "%s%s_peaks.txt"%(peak_folder, base_name_test)
	#Make a name for the rewritten peak file
	rewrite_peak = "%s%s_peaks_rewrite.txt"%(peak_folder, base_name_test)
	
	with open(output_peak, 'r') as fo:
		data = [line.rstrip() for line in fo]
		data = [line.split() for line in data if line]
	
	with open(rewrite_peak, 'w') as fo:
		for line in data:
			#Homers peak.txt file starts with lines that report on the
			#used parameters. All these lines start with a #
			#So just copy those directly in the new file
			if line[0][0] == "#":
				for item in line:
					fo.write("%s "%item)
				fo.write("\n")
				
			else:
				#Next follows the actual peak data given in 15 columns
				#The chromosome number is given in the second column
				for i in range(14):
					if i != 1:
						fo.write("%s\t"%line[i])
					elif line[i][0:3] != 'chr':
						fo.write("chr%s\t"%line[i])
					else:
						fo.write("%s\t"%line[i])
				fo.write("\n")
				
				




def HOMER_make_motifs(bam_test_file, peak_folder, motif_folder):
	"""
	findMotifsGenome.pl 	/home/wouter/chipseq_project/homer/findpeak_output/peaks2.txt 
							hg19 
							/home/wouter/chipseq_project/homer/motif_output2/ 
							-size given
	"""
	
	#Prepare the BED file names for input
	#take everything except the .bam extension
	base_name_test = os.path.splitext(bam_test_file)[0]

	#Make a name for the rewritten peak file
	rewrite_peak = "%s%s_peaks_rewrite.txt"%(peak_folder, base_name_test)
	
	#Make a name for the motif folder
	motif_output = "%s%s_motifoutput/"%(motif_folder, base_name_test)
	
	cmd = ['findMotifsGenome.pl', rewrite_peak, 'hg19', motif_output, '-size', 'given']
	
	print "Start finding motifs for %s"%base_name_test
	#pipe = check_output(cmd1)
	pipe = Popen(cmd, stdout=PIPE, stderr=PIPE)
	stdout, stderr = pipe.communicate()








