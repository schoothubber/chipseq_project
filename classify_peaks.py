#!/usr/bin/env python2.7

######################################
###Author: W. van der Schoot
######################################

import os
import matplotlib
#use a non-interactive backend
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import arange as ar
from collections import Counter

from prepare_data import get_base_name, make_random_peak_list
from get_near_genes import neighbours


def open_master_file(neighbor_folder, bam_test_file):
	"""
	Use this function to read and parse the Master file
	"""
	
	base_name_test = get_base_name(bam_test_file)
	neighbour_file = "%sneighbour_%s_peak.cod"%(
								neighbor_folder, 
								base_name_test
								)

	with open(neighbour_file, 'r') as fo:
		
		#skip the first line which containes the headers
		#next(fo)
		#remove all empty lines
		data = [line.rstrip() for line in fo]
		data = [line.split() for line in data if line]
		
		header = data[0]
		data = data[1:]
		
	return header, data, base_name_test


def peak_classifier(fileargs):
	"""
	classify all peaks based on distance from TSS
	but only for specific distances
	"""
	#folder parameters
	bam_test_file = fileargs['bam_test_file']
	neighbor_folder = fileargs['neighbor_folder']
	
	#classify genes based on TSS vicinity :2kb, 5kb, 20kb, 100kb
	vics = [2000, 5000, 20000, 100000]
	
	#read data
	header, data, base_name = open_master_file(
								neighbor_folder, 
								bam_test_file
								)
	
	for vic in vics:
		
		fn = "%s%s_%skb.txt"%(neighbor_folder, base_name, str(vic/1000))
		with open(fn, 'w') as fo:
			for head in header:
				fo.write("%s\t"%head)
			fo.write("\n")
		
		for line in data:
			try:
				if abs(int(line[8])) <= vic:
					with open(fn, 'a') as fo:
						for item in line:
							fo.write("%s\t"%item)
						fo.write("\n")
					
			except ValueError:
				pass


def peak_distributor(fileargs, peakargs, annoargs):
	"""
	classify all peaks based on distance from TSS
	do the same with random regions
	
	display data in a bar graph
	
	"""
	#folder parameters
	bam_test_file = fileargs['bam_test_file']
	bam_ctrl_file = fileargs['bam_ctrl_file']
	aln_folder = fileargs['aln_folder']
	neighbor_folder = fileargs['neighbor_folder']
	random_folder = fileargs['random_folder']
	annotation_file = fileargs['nanno']
	random_folder = fileargs['random_folder']
	
	#read data
	header, data, base_name = open_master_file(
								neighbor_folder, 
								bam_test_file
								)
	
	TSS_data_min = get_TSS_data(data)
	peak_data_perc = classify(TSS_data_min)
	N = len(peak_data_perc)
	#print "N: %s"%N
	
	######################################
	####MAKE ROOM FOR RANDOM VALUES!!!####
	######################################
	#The random peaks are taken from the input file
	#At least 100x as many random peaks than there are actual peaks
	#But never more than 1M peaks
	
	base_name_ctrl = get_base_name(bam_ctrl_file)
	make_random_peak_list(peakargs, random_folder, aln_folder, base_name_ctrl, N)
	random_path = "%s%s"%(random_folder, "random_data.cod")
	
	#annotate regions
	#The random peaks need to be annotated in exactly the same way
	#as the actual peaks
	rand_fileargs = {
			'bam_test_file' : random_path,
			'annotation_file' : annotation_file,
			'seqpeak_folder' : random_folder,
			'neighbor_folder' : random_folder,
			'cisgenome_path' : fileargs['cisgenome_path']
			}
	neighbours(rand_fileargs, annoargs)
	
	#read random data
	random_peak_out = "%s%s"%(random_folder, "random_data.cod")
	header, random_data, base_name = open_master_file(
							random_folder, 
							random_peak_out
							)
	
	random_TSS_data_min = get_TSS_data(random_data)
	random_peak_data_perc = classify(random_TSS_data_min)
	
	
	########################
	#matplotlib stuff
	########################
	
	ind = ar(N)  # the x locations for the groups
	width = 0.35       # the width of the bars

	fig, ax = plt.subplots()
	#the graph will have double bars: 
	#red for the actual peaks
	#and yellow for the random peaks
	
	#print "######################"
	#print len(peak_data_perc)
	#print len(ind)
	#print peak_data_perc
	#print len(random_peak_data_perc)
	#print random_peak_data_perc
	
	bars_actual = ax.bar(ind, peak_data_perc, width, color='r')
	bars_random = ax.bar(ind+width, random_peak_data_perc, width, color='y')
	#print "################################"
	#print len(bars_actual)
	#print bars_actual
	#print len(bars_random)
	#print bars_random
	
	# add some text for labels, title and axes ticks
	plt.xticks(rotation=90)
	plt.tight_layout()
	
	fig = plt.gcf()
	fig.subplots_adjust(bottom = 0.25, left = 0.1, top = 0.8)
	
	ax.yaxis.grid(True)
	ax.set_ylabel('% of peaks')
	ax.set_xlabel('Distance to closest TSS (kb)')
	ax.set_title('Peak Distribution')
	ax.set_xticks(ind+width)
	ax.set_xticklabels( (
				'-200 to -100', '-100 to -50', '-50 to -20', 
				'-20 to -10', '-10 to -5', '-5 to -2', '-2 to -1', 
				'-1 to -0', '0 to 1', '1 to 2', '2 to 5', '5 to 10', 
				'10 to 20', '20 to 50','50 to 100', '100 to 200'
						) )

	ax.legend((bars_actual[0], bars_random[0]), 
				('Peaks', 'Random'),
				loc = 'upper center',
				bbox_to_anchor = (0.5,1.35))

	plot_name = "%s%s%s"%(
						neighbor_folder, base_name,
						'peak_distribution.png'
						)
						
	print "saving %s"%plot_name
	plt.savefig(plot_name)
	plt.close()


def get_TSS_data(data):
	"""
	Receive data of peaks and their distances to TSS
	And return the nearest TSS for each peak
	"""
	TSS_data_max = []
	for line in data:
		try:
			id_tss = [line[0],int(line[8])]
			TSS_data_max.append(id_tss)
		except ValueError:
			#sometime line[8] is a string
			#we dont need those...
			pass
			
	TSS_data_min = []
	temp_tss_list = []
	TSS_list = []
	
	for seq_id, TSS in TSS_data_max:
	#there are blocks for each peak with similar id numbers
	#from these blocks select the smallest absolute TSS
		temp_tss = TSS
		if not temp_tss_list:
			identifier = seq_id
			TSS_list.append(TSS)
			temp_tss_list.append(abs(temp_tss))
		elif identifier == seq_id:
			TSS_list.append(TSS)
			temp_tss_list.append(abs(temp_tss))
		elif identifier != seq_id:
			#print len(TSS_list)
			#print len(temp_tss_list)
			min_tss = min(temp_tss_list)
			tss_index = temp_tss_list.index(min_tss)
			TSS_data_min.append(TSS_list[tss_index])
			#reset the lists
			temp_tss_list = []
			TSS_list = []
			#
			identifier = seq_id
			TSS_list.append(TSS)
			temp_tss_list.append(abs(temp_tss))
	
	return TSS_data_min

		
		
def classify(data):
	"""
	Arrange the peaks based on their distances from the TSS sites
	"""
	
	minus200 = 0
	minus100 = 0
	minus50 = 0
	minus20 = 0
	minus10 = 0
	minus5 = 0
	minus2 = 0
	minus1 = 0
	plus200 = 0
	plus100 = 0
	plus50 = 0
	plus20 = 0
	plus10 = 0
	plus5 = 0
	plus2 = 0
	plus1 = 0
		
	for TSS in data:
		
		try:
			TSS = int(TSS)
					
			if TSS >= -200000 and TSS < -100000:
				minus200 += 1
			if TSS >= -100000 and TSS < -50000:
				minus100 += 1
			if TSS >= -50000 and TSS < -20000:
				minus50 += 1
			if TSS >= -20000 and TSS < -10000:
				minus20 += 1
			if TSS >= -10000 and TSS < -5000:
				minus10 += 1
			if TSS >= -5000 and TSS < -2000:
				minus5 += 1
			if TSS >= -2000 and TSS < -1000:
				minus2 += 1
			if TSS >= -1000 and TSS < 0:
				minus1 += 1
				
			if TSS >= 0 and TSS < 1000:
				plus1 += 1
			if TSS >= 1000 and TSS < 2000:
				plus2 += 1
			if TSS >= 2000 and TSS < 5000:
				plus5 += 1
			if TSS >= 5000 and TSS < 10000:
				plus10 += 1
			if TSS >= 10000 and TSS < 20000:
				plus20 += 1
			if TSS >= 20000 and TSS < 50000:
				plus50 += 1
			if TSS >= 50000 and TSS < 100000:
				plus100 += 1
			if TSS >= 100000 and TSS < 200000:
				plus200 += 1
			
			
		except ValueError:
			pass	
			
	peak_data = [
				minus200, minus100, minus50, minus20, minus10, 
				minus5, minus2, minus1, plus1, plus2, plus5,
				plus10, plus20, plus50, plus100, plus200
				]
	
	total_peaks = sum(peak_data)
	
	#the results will be presented in percentages
	peak_data_perc = [(float(peak)/float(total_peaks))*100 for peak in peak_data]
	
				
	return peak_data_perc


def peak_localizer(fileargs):
	"""
	Divide the peaks into 2 groups: 1;intergenic, 2;intragenic
	
	Based on each line in the data from the masterfile...
	Ascertain whether the physical locations in line[2] and line[3]...
	Are in- or outside the physical locations in line[13] and line[14]
	
	But make sure to only measure eahc peak once
	Each peak has only 1 seq_id (column 0)
	
	The results are saved in a pie chart using matplotlib
	"""
	#folder parameters
	bam_test_file = fileargs['bam_test_file']
	neighbor_folder = fileargs['neighbor_folder']
	
	header, data, base_name = open_master_file(
							neighbor_folder, 
							bam_test_file
							)
	
	intragenic = False
	identifier = None
	
	inter_count = 0
	intra_count = 0
	
	for line in data:
		
		try:
			seq_id = line[0]	#identifier
			ps = int(line[2])	#peak start
			pe = int(line[3])	#peak_end
			gs = int(line[13])	#gene transcription start site
			ge = int(line[14])	#gene transcription end site
			
			if identifier == seq_id:
				if ps > gs and pe < ge:
					intragenic = True
				else:
					pass
					
			elif identifier != seq_id:
				if ps > gs and pe < ge:
					intragenic = True
				if intragenic:
					intra_count += 1
				else:
					inter_count += 1
				identifier = seq_id
				intragenic = False
				
			elif identifier == None:
				identifier = seq_id
			
		except ValueError:
			pass
		
		#if ps > gs and pe < ge:
		#	intragenic += 1
			
		#if ps < gs and pe < gs:
		#	intergenic += 1
			
		#if ps > ge and pe > ge:
		#	intergenic +=1
	
	total = inter_count + intra_count
	inter_perc = float(inter_count) / total * 100
	intra_perc = float(intra_count) / total * 100
	
	
	#Prepare to display the calculated values with a pie chart
	labels = 'intergenic', 'intragenic'
	sizes = [inter_perc, intra_perc]
	colors = ['red', 'blue']
	
	plt.pie(sizes, labels = labels, colors = colors, autopct='%1.1f%%',
			startangle = 90)
	plt.axis('equal')
	
	pie_name = "%s%s%s"%(neighbor_folder, base_name,'peak_pie.png')
	print "saving %s"%pie_name
	plt.savefig(pie_name)
	plt.close()	
	
		

	
def peaks_per_chromosome(fileargs):
	"""
	As the function titel kind of gives away...
	this function counts the peaks per chromosome
	And saves the results in a text file and a graph
	
	The peaks need to be counted/arranged per 1Mb window
	Creating a so called density plot
	"""
	#folder parameters
	bam_test_file = fileargs['bam_test_file']
	neighbor_folder = fileargs['neighbor_folder']
	seqpeak_folder = fileargs['seqpeak_folder']
	
	#the window size is used to calculate the peak density
	#the number of peaks is counted for each window
	#default = 1Mb
	window = 1000000
	
	#read data
	base_name_test = get_base_name(bam_test_file)
	seqpeak_file = "%s%s_peak.cod"%(seqpeak_folder, base_name_test)
	
	#seqpeak_file = "%s%s"%(seqpeak_folder, "LPSIL102H4_dedupunique_peak.cod")
	
	#print seqpeak_file

	with open(seqpeak_file, 'r') as fo:
		
		fo.next()
		#remove all empty lines
		data = [line.rstrip() for line in fo]
		data = [line.split() for line in data if line]
	
	#place all peaks in a dictionary
	#with key = chromosome number
	#and value = peak start site	
	chrompeak_dict = {}
	for line in data:
		
		try:
			#remove the 'chr' and keep the integer
			#so that sorting can be done
			chrom = int(line[1][3:])
			peak_start = int(line[2])
			
			if not chrompeak_dict.has_key(chrom):
				chrompeak_dict[chrom] = [int(peak_start)]
			else:
				chrompeak_dict[chrom].append(int(peak_start))
			
		except ValueError:
			
			#in case there is no integer, but a X or Y
			chrom = line[1][3:]
			peak_start = int(line[2])
			
			if not chrompeak_dict.has_key(chrom):
				chrompeak_dict[chrom] = [int(peak_start)]
			else:
				chrompeak_dict[chrom].append(int(peak_start))
	
	#make a list of the keys (chromosomes) so that they can be 
	#called in order from the dictionary
	#
	#I know there is a sortedDictionary method but im to lazy to import
	key_chroms = sorted(chrompeak_dict.keys())

	#go through all peaks per chromosome
	counts_dict = {}
	for chrom in key_chroms:
		#chrom = 'chr' + str(chrom)
		#for each chromosome there will be a list of peak start sites
		peak_starts = sorted(chrompeak_dict[chrom])
		
		#A little list comprehension to lighten the mood
		windows = [peak_start/window for peak_start in peak_starts]
		
		counts = Counter(windows)
		counts_dict[chrom] = counts
		#print counts_dict

	
	#write the data to a text file, because of reasons
	fn = "%s/%s%s"%(neighbor_folder, base_name_test, '_peakdensity.txt')
	with open(fn, 'w') as fo:

		for chrom, counts in counts_dict.iteritems():
			chrom = 'chr' + str(chrom)
			fo.write("%s\n"%chrom)
				
			for window, tally in counts.iteritems():
				fo.write("%s\t%s\n"%(window, tally))
			fo.write("\n")
	
	#################################
	#######Matplotlib stuff##########
	#################################
	#write the data into a graph plot, because of more reasons
	plot_name = "%s/%s%s"%(neighbor_folder, base_name_test, '_peakdensity.png')
	
	x = [] # x coordinates
	y = [] # y coordinates
	v = [] # x coordinates for vertical lines
	xlabels = []
	
	i = 0
	for chrom, counts in counts_dict.iteritems():
		
		v.append(i)
		#create the locations for the labels
		for window, tally in counts.iteritems():
			#create the x and y coordinates
			x.append(i)
			y.append(tally)
			i += 1
		
		#create a label
		xlabels.append('chr' + str(chrom))

	v.append(max(x))
	
	#make coordinates for the xlabels
	#which should be in between 2 consecutive X values
	xlabel_locs = []
	for j in range(0, len(v)-1):
		#print j
		loc = (float(v[j]) + float(v[j+1]))/2.0
		xlabel_locs.append(loc)	

	#print v
	#print len(v)
	#print xlabel_locs
	#print len(xlabel_locs)
	#print xlabels
	#print len(xlabels)
	#print len(x)
	#print len(y)
	#print xlabels
	
	plt.figure(figsize=(20,5))
	plt.xlim(0, len(x))
	plt.ylim(0, max(y)+5)
	
	plt.title('peak density')
	plt.xlabel('chromosomes')
	plt.ylabel('peak counts per Mb')
	
	plt.xticks(xlabel_locs, xlabels, rotation = 'vertical')
	plt.subplots_adjust(bottom = 0.25)
	plt.vlines(v, 0, max(y)+5, colors = 'y')
	plt.scatter(x,y, s=5)
	plt.savefig(plot_name)
	plt.close()
	
	
	
