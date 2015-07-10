
from prepare_data import get_base_name




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
		'log_folder' : log_folder
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
	
	bam_test = get_base_name(fileargs["bam_test_file"])
	bam_ctrl = get_base_name(fileargs["bam_ctrl_file"])
	
	with open(fn, 'w') as fo:
		
		fo.write("#parameters used for:\n")
		fo.write("bam test file: %s\n"%bam_test)
		fo.write("bam control file: %s\n"%bam_ctrl)
		fo.write("\n")
		
		#Write the seqpeak parameters
		fo.write("#Peak calling parameters used:\n")
		fo.write("seqpeak_readextension:\t%s\n "%peakargs["seqpeak_readextension"])
		fo.write("seqpeak_binsize:\t%s\n "%peakargs["seqpeak_binsize"])
		fo.write("seqpeak_halfwinsize:\t%s\n "%peakargs["seqpeak_halfwinsize"])
		fo.write("seqpeak_standardize:\t%s\n "%peakargs["seqpeak_standardize"])
		fo.write("seqpeak_maxpeakgap:\t%s\n "%peakargs["seqpeak_maxpeakgap"])
		fo.write("seqpeak_minreglen:\t%s\n "%peakargs["seqpeak_minreglen"])
		fo.write("seqpeak_boundary:\t%s\n "%peakargs["seqpeak_boundary"])
		fo.write("\n")
		
		#write the neighbor parameters
		fo.write("#Gene annotation parameters used:\n")
		fo.write("neighbor_species:\t%s\n "%annoargs["neighbor_species"])
		fo.write("neighbor_gdistance:\t%s\n "%annoargs["neighbor_gdistance"])
		fo.write("neighbor_upstreamgenes:\t%s\n "%annoargs["neighbor_upstreamgenes"])
		fo.write("neighbor_downstreamgenes:\t%s\n "%annoargs["neighbor_downstreamgenes"])
		fo.write("\n")
		
		#write motif parameters
		fo.write("#De novo motif analysis parameters used:\n")
		fo.write("motif_genome:\t%s\n "%motifargs["motif_genome"])
		fo.write("motif_size:\t%s\n "%motifargs["motif_size"])
		fo.write("motif_length:\t%s\n "%motifargs["motif_length"])
		fo.write("motif_number:\t%s\n "%motifargs["motif_number"])
		fo.write("motif_mismatch:\t%s\n "%motifargs["motif_mismatch"])
		fo.write("\n")















		
