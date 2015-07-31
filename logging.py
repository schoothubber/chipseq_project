


def parameter_logger(fileargs, peakargs, annoargs, motifargs):
	"""
	The parameters need to be stored in the log_folder
	
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
		'sambamba_path' : args.sam,
		'cisgenome_path' : args.cis,
		'homer_path' : args.hom,
		'genome_path' : args.g
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
	"""
	
	fn = '%s/paramaters.log'%fileargs['log_folder']
	
	
	with open(fn, 'w') as fo:
	
		#write used files and folders
		fo.write("#parameters used for:\n")
		fo.write("bam test file: %s\n"%fileargs["bam_test_file"])
		fo.write("bam control file: %s\n"%fileargs["bam_ctrl_file"])
		fo.write("\n")
		fo.write("main output folder: %s\n"%fileargs['main_output'])
		fo.write("ALN folder: %s\n"%fileargs['aln_folder'])
		fo.write("seqpeak folder: %s\n"%fileargs['seqpeak_folder'])
		fo.write("Peak annotation folder: %s\n"%fileargs['neighbor_folder'])
		fo.write("Motif analysis folder: %s\n"%fileargs['motif_folder'])
		fo.write("Random data folder: %s\n"%fileargs['random_folder'])
		fo.write("\n")
		
		#write used tools
		fo.write("#tools used:\n")
		fo.write("sambamba path: %s\n"%fileargs['sambamba_path'])
		fo.write("cisgenome path: %s\n"%fileargs['cisgenome_path'])
		fo.write("homer path: %s\n"%fileargs['homer_path'])
		fo.write("\n")
		
		#Write the seqpeak parameters
		fo.write("#peak calling parameters used:\n")
		fo.write("seqpeak_readextension:\t%s\n"%peakargs["seqpeak_readextension"])
		fo.write("seqpeak_binsize:\t%s\n"%peakargs["seqpeak_binsize"])
		fo.write("seqpeak_halfwinsize:\t%s\n"%peakargs["seqpeak_halfwinsize"])
		fo.write("seqpeak_standardize:\t%s\n"%peakargs["seqpeak_standardize"])
		fo.write("seqpeak_maxpeakgap:\t%s\n"%peakargs["seqpeak_maxpeakgap"])
		fo.write("seqpeak_minreglen:\t%s\n"%peakargs["seqpeak_minreglen"])
		fo.write("seqpeak_boundary:\t%s\n"%peakargs["seqpeak_boundary"])
		fo.write("\n")
		
		#write the neighbor parameters
		fo.write("#gene annotation parameters used:\n")
		fo.write("neighbor_species:\t%s\n"%annoargs["neighbor_species"])
		fo.write("neighbor_gdistance:\t%s\n"%annoargs["neighbor_gdistance"])
		fo.write("neighbor_upstreamgenes:\t%s\n"%annoargs["neighbor_upstreamgenes"])
		fo.write("neighbor_downstreamgenes:\t%s\n"%annoargs["neighbor_downstreamgenes"])
		fo.write("\n")
		
		#write motif parameters
		fo.write("#de novo motif analysis parameters used:\n")
		fo.write("motif_genome:\t%s\n"%motifargs["motif_genome"])
		fo.write("motif_size:\t%s\n"%motifargs["motif_size"])
		fo.write("motif_length:\t%s\n"%motifargs["motif_length"])
		fo.write("motif_number:\t%s\n"%motifargs["motif_number"])
		fo.write("motif_mismatch:\t%s\n"%motifargs["motif_mismatch"])
		fo.write("\n")















		
