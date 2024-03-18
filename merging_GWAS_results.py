#!/usr/bin/env python
#import
import pandas as pd
import os
import argparse
import numpy as np

#define arguments
def combine_and_filter(input_folder
                       , output_file
                       , file_prefix
                       , file_suffix
                       , min_af=0.01
                       , min_info=0.3
                       , num_chrom=22
                      ):
	# Create an empty DataFrame to store the combined data
	combined_df = pd.DataFrame()
	# Loop through the specified number of chromosomes
	for chrom in range(21, num_chrom + 1):
		print("for chromosome",chrom)
		# Construct the filename using the provided prefix and suffix
		filename = os.path.join(input_folder, f"{file_prefix}{chrom}{file_suffix}")
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

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Combine and filter GWAS summary statistics.")

	parser.add_argument("--input_folder", required=True, help="folder containing GWAS summary statistics")
	parser.add_argument("--output_file", required=True, help="name for output file")
	parser.add_argument("--file_prefix", required=True, help="prefix before the chromosome number in the GWAS file names")
	parser.add_argument("--file_suffix", required=True, help="suffix after the chromosome number in the GWAS file names")
	parser.add_argument("--min_af", required=False, type=float, default=0.01, help="minimum allele frequency")
	parser.add_argument("--min_info", required=False, type=float, default=0.3, help="minimum info score")
	parser.add_argument("--num_chrom", required=False, type=int, default=22, help="number of chromosomes" )

	# Parse the command-line arguments
	args = parser.parse_args()
  
	# Call the function to combine and filter the data
	combine_and_filter(args.input_folder, args.output_file, args.file_prefix, args.file_suffix, args.min_af,  args.min_info, args.num_chrom)

#these are my hanges to the script.
