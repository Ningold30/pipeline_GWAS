#!/usr/bin/env python
#import
import pandas as pd
import os
import argparse
import numpy as np

#define arguments
def combine_and_filter(input_folder
                       , output_file
                       , file_chr_22
                       , min_af=0.01
                       , min_info=0.3
                       , num_chrom=22
                      ):
	# Create an empty DataFrame to store the combined data
	combined_df = pd.DataFrame()
	# deconsruct the file name
	prefix,suffix = file_chr_22.split("22")
	# Loop through the specified number of chromosomes
	for chrom in range(21, num_chrom + 1):
		print("for chromosome",chrom)
		# Construct the filename using the provided prefix and suffix
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
	#save data frame to output_file
	output_name= os.path.join(input_folder, output_file)
	combined_df.to_csv(output_name, index=False, sep = " ", na_rep ="NA")
	print(f"File saved to: {output_name}")
	#format for LDSC
	LDSC = combined_df.loc[:, ["ID","ALLELE1","ALLELE0","BETA","p_value","N"]]
	LDSC.columns = ["SNP","A1","A2","BETA","P","N"]
	LDSC.to_csv(f"{output_name}.LDSC",index=False, sep = " ", na_rep = "NA")
	print(f"LDSC File saved to: {output_name}.LDSC")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Combine and filter GWAS summary statistics.")
	parser.add_argument("--input_folder", required=True, help="folder containing GWAS summary statistics")
	parser.add_argument("--output_file", required=True, help="name for output file")
	parser.add_argument("--file_chr_22", required=True, help="the full file name of chromsome 22 - used to work out naming convention of input and output files")
	parser.add_argument("--min_af", required=False, type=float, default=0.01, help="minimum allele frequency, default = 0.01")
	parser.add_argument("--min_info", required=False, type=float, default=0.3, help="minimum info score, default = 0.3")
	parser.add_argument("--num_chrom", required=False, type=int, default=22, help="number of chromosomes, default = 22" )

	# Parse the command-line arguments
	args = parser.parse_args()
  
	# Call the function to combine and filter the data
	combine_and_filter(args.input_folder, args.output_file, args.file_chr_22, args.min_af,  args.min_info, args.num_chrom)

