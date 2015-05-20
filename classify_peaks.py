#!/usr/bin/env

import os
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter


def open_master_file(neighbor_folder, bam_test_file):
	"""
	Use this function to read and parse the Master file
	"""
	
	base_name_test = os.path.splitext(bam_test_file)[0]
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


def peak_classifier(neighbor_folder, bam_test_file):
	"""
	classify all peaks based on distance from TSS
	but only for specific distances
	"""
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


def peak_distributor(neighbor_folder, bam_test_file):
	"""
	classify all peaks based on distance from TSS
	do the same with random regions
	
	display data in a bar graph
	
	"""
	
	#read data
	header, data, base_name = open_master_file(
											neighbor_folder, 
											bam_test_file
												)
	
	TSS_data = [line[8] for line in data]
	peak_data_perc = classify(TSS_data)

	#Start matplotlib stuff
	
	N = len(peak_data_perc)

	ind = np.arange(N)  # the x locations for the groups
	width = 0.35       # the width of the bars

	fig, ax = plt.subplots()
	rects1 = ax.bar(ind, peak_data_perc, width, color='r')
	
	#RESERVED FOR RANDOM VALUES
	random_regions = np.random.random_integers(-200000, 200000, len(data))
	random_data = classify(random_regions)
	rects2 = ax.bar(ind+width, random_data, width, color='y')

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

	ax.legend((rects1[0], rects2[0]), 
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


def peak_localizer(neighbor_folder, bam_test_file):
	"""
	Divide the peaks into 2 groups: 1;intergenic, 2;intragenic
	
	Based on each line in the data from the masterfile...
	Ascertain whether the physical locations in line[2] and line[3]...
	Are in- or outside the physical locations in line[13] and line[14]
	
	The results are saved in a pie chart using matplotlib
	"""
	
	header, data, base_name = open_master_file(
											neighbor_folder, 
											bam_test_file
												)
	
	intergenic = 0
	intragenic = 0
	
	for line in data:
		
		try:
			ps = int(line[2])	#peak start
			pe = int(line[3])	#peak_end
			gs = int(line[13])	#gene transcription start site
			ge = int(line[14])	#gene transcription end site
			
		except ValueError:
			pass
		
		
		if ps > gs and pe < ge:
			intragenic += 1
			
		if ps < gs and pe < gs:
			intergenic += 1
			
		if ps > ge and pe > ge:
			intergenic +=1
	
	total = intergenic + intragenic
	inter_perc = float(intergenic) / total * 100
	intra_perc = float(intragenic) / total * 100
	
	
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
	
		
def peaks_per_chromosome(seqpeak_folder, neighbor_folder, bam_test_file):
	"""
	As the function titel kind of gives away...
	this function counts the peaks per chromosome
	And saves the results in a text file
	"""
	#read data
	base_name_test = os.path.splitext(bam_test_file)[0]
	seqpeak_file = "%s%s_peak.cod"%(seqpeak_folder, base_name_test)

	with open(seqpeak_file, 'r') as fo:
		
		fo.next()
		#remove all empty lines
		data = [line.rstrip() for line in fo]
		data = [line.split() for line in data if line]
			
	chrom_list = []
	for line in data:
		
		try:
			chrom = line[1]
			chrom_list.append(chrom)
			
		except ValueError:
			pass
	
	#count all occurences of a chromosome number in a list
	chrom_counts = Counter(chrom_list)
	
	#extract the produced data from the dictionary and store in a list
	#so that the output can but saved in a sorted order
	cc = []
	for k,v in chrom_counts.iteritems():
		cc.append([k,v])
	sorted_cc = sorted(cc)
	
	#write the data to a text file
	fn = "%s/%s%s"%(neighbor_folder, base_name_test, '_peaksperchrom.txt')
	with open(fn, 'w') as fo:

		for chrom, count in sorted_cc:
			fo.write("%s\t%s\n"%(chrom,count))	
	
	
	
	
	


	
