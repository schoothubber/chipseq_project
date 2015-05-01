import os, fnmatch
import matplotlib.pyplot as plt
import numpy as np



def peak_classifier(neighbor_folder):
	"""
	"""
	pattern = "*.cod"
	neighbour_files = [file for file in os.listdir(neighbor_folder) if fnmatch.fnmatch(file, pattern)]
	
	#classify genes based on TSS vicinity :2kb, 5kb, 20kb, 100kb
	vics = [2000, 5000, 20000, 100000]
	#filter each line based on distance of TSS
	#distance to TSS is column 8 (counting starts from 0)		
	for neighbour_file in neighbour_files:
		
		neighbour_file = "%s%s"%(neighbor_folder, neighbour_file)
		#open 4 separate files for each vic
		base_name = os.path.splitext(neighbour_file)[0]
		vic_files = []
		for vic in vics:
			vic_file = "%s_%skb.txt"%(base_name, str(vic/1000))
			vic_files.append(vic_file)
			open (vic_file, 'w').close()
				#do nothing yet
			
		
		#open a master file
		with open(neighbour_file, 'r') as fo:
			
			#skip the first line which containes the headers
			next(fo)
			#remove all empty lines
			data = [line.rstrip() for line in fo]
			data = [line for line in data if line]
			
			#start iteration through all lines
			for line in data:
				info = line.split()
				
				#compare each vic for each line
				#and write it to the corresponding file
				for vic in vics:
					try:
						if abs(int(info[8])) <= vic:
							vic_file = "%s_%skb.txt"%(
													base_name, 
													str(vic/1000)
														)
							with open(vic_file, 'a') as fo:
								fo.write("%s"%line)
								fo.write("\n")
							
					except ValueError:
						pass



def peak_distributor(neighbor_folder):
	"""
	classify all peaks based on distance from TSS
	do the same with random regions
	
	display data in a bar graph
	
	"""
	pattern = "*.cod"
	neighbour_files = [file for file in os.listdir(neighbor_folder) if fnmatch.fnmatch(file, pattern)]
	

	for neighbour_file in neighbour_files:
		
		neighbour_file = "%s%s"%(neighbor_folder, neighbour_file)
		
		with open(neighbour_file, 'r') as fo:
			
			#skip the first line which containes the headers
			next(fo)
			#remove all empty lines
			data = [line.rstrip() for line in fo]
			data = [line.split() for line in data if line]
			
			data = [line[8] for line in data]

		peak_data_perc = classify(data)

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
		fig.subplots_adjust(bottom = 0.25, left = 0.1)
		
		ax.yaxis.grid(True)
		ax.set_ylabel('% of peaks')
		ax.set_xlabel('Distance to closest TSS (kb)')
		ax.set_title('Peak Distribution')
		ax.set_xticks(ind+width)
		ax.set_xticklabels( (
					'-200 to -100', '-100 to -50', '-50 to -20', '-20 to -10',
					 '-10 to -5', '-5 to -2', '-2 to -1', '-1 to -0', '0 to 1',
					 '1 to 2', '2 to 5', '5 to 10', '10 to 20', '20 to 50',
					 '50 to 100', '100 to 200'
							) )

		ax.legend((rects1[0], rects2[0]), ('Peaks', 'Random'))

		plot_name = "%s%s"%(neighbor_folder,'peak_distribution.png')
		plt.savefig(plot_name)	
		
		
def classify(data):
	"""
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
	peak_data_perc = [(float(peak)/float(total_peaks))*100 for peak in peak_data]
	
				
	return peak_data_perc
	
