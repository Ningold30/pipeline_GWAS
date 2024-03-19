#!/usr/bin/env python

##########################
#import
##########################
import pandas as pd
import os
import argparse
import numpy as np
import subprocess

##########################
#define arguments
##########################
#munging summstats
def munge_results(input_file,input_dir):
	print("--Munge set to True")
	munge_script_path= "/mnt/backedup/home/nathanI/Scripts/LDSC/munge.sh"
	subprocess.run(['qsub', '-v', f'file1={input_file},dirR={input_dir}', munge_script_path])
	subprocess.run(['qstat','-u nathanI'])

#plotting sum stats
def plotting_man_qqplot(input_file,input_dir,title):
	print("--Plot set to True")
	plot_script_path= "/mnt/backedup/home/nathanI/Scripts/pipeline_GWAS/qsub_ManQQplot.sh"
	subprocess.run(['qsub', '-v', f'file1={input_file},dirR={input_dir},title={title}', plot_script_path])
	print("plotting")
	subprocess.run(['qstat','-u nathanI'])
	

#combining, filtering  and processing sumstats
def combine_and_filter(input_folder
			, output_file
			, file_chr_22
			, min_af=0.01
			, min_info=0.3
			, num_chrom=22
			, munge=True
			, ManQQplot=False
			, ManQQplot_title=None
			):
	# Create empty DataFrame
	combined_df = pd.DataFrame()
	# deconsruct the file name
	prefix,suffix = file_chr_22.split("22")
	# Loop through the specified number of chromosomes
	for chrom in range(21, num_chrom + 1):
		print("for chromosome",chrom)
		# Construct the filename using pre and suffix
		filename = os.path.join(input_folder, f"{prefix}{chrom}{suffix}")
		print(f"Reading file: {filename}")
		# Read the GWAS summary statistics file for each chromosome
		df = pd.read_csv(filename, sep = "\s+")
		# Filter rows based on AF and INFO
		print("Line count before filtering",df.shape[0])
		df = df[(df['A1FREQ'] >= min_af) & (df['A1FREQ'] <= 1-min_af) & (df['INFO'] >= min_info)]
		print(f"Line count after filtering:, AF {min_af}, INFO {min_info}", df.shape[0])
		# convert log10 pval to pvalue
		df['p_value'] = 10 ** (-df['LOG10P'])
		#ammend current df to combined_df
		combined_df = pd.concat([combined_df, df], ignore_index=True)
		print("Line count of combined:", combined_df.shape[0])
		print(combined_df.head)
	#save
	output_name= os.path.join(input_folder, output_file)
	combined_df.to_csv(output_name, index=False, sep = " ", na_rep ="NA")
	print(f"File saved to: {output_name}")
	
	
	if munge:
		#format for LDSC
		LDSC = combined_df.loc[:, ["ID","ALLELE1","ALLELE0","BETA","p_value","N"]]
		LDSC.columns = ["SNP","A1","A2","BETA","P","N"]
		LDSC_output = f"{output_name}.LDSC"		
		LDSC.to_csv(f"{LDSC_output}",index=False, sep = " ", na_rep = "NA")
		print(f"LDSC files saved here: {LDSC_output}")
		LDSC_path, LDSC_file = os.path.split(LDSC_output)
		munge_results(LDSC_file, LDSC_path )
		
	if ManQQplot:
		plotting_path,plotting_file = os.path.split(output_name)
		plotting_man_qqplot(plotting_file,plotting_path,ManQQplot_title)
		

###########################
#Parse arguments	
###########################
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="This Script combines and filters GWAS summary statistics across chromsomes, allows you to munge the summary statistics, and allows you to plot a Mannhattan plot and qqplot using the --ManQQplot and --ManQQplot_title")
	parser.add_argument("--input_folder", required=True, help="folder containing GWAS summary statistics")
	parser.add_argument("--output_file", required=True, help="name for output file")
	parser.add_argument("--file_chr_22", required=True, help="the full file name of chromsome 22 - used to work out naming convention of input and output files")
	parser.add_argument("--min_af", required=False, type=float, default=0.01, help="minimum allele frequency, default = 0.01")
	parser.add_argument("--min_info", required=False, type=float, default=0.3, help="minimum info score, default = 0.3")
	parser.add_argument("--num_chrom", required=False, type=int, default=22, help="number of chromosomes, default = 22" )
	parser.add_argument("--munge", required=False, action="store_true", default=True, help="do you want to munge the sumstats for LDSC, default = True")
	parser.add_argument("--ManQQplot", required=False, action="store_true", default=False, help="do you want to plot a Manhattan and QQ plt, default = False, must be accompanied with the --ManQQplot_title flag")
	parser.add_argument("--ManQQplot_title", required=False, help="provide a title for the plots")

	# Parse the command-line arguments
	args = parser.parse_args()
  
	# Call the function to combine and filter the data
	combine_and_filter(args.input_folder, args.output_file, args.file_chr_22, args.min_af,  args.min_info, args.num_chrom, args.munge, args.ManQQplot)

