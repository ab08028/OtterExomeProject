import os
import sys
import argparse
import csv
import gzip

def parse_args():
	"""
	Parse command-line arguments
	"""
	parser = argparse.ArgumentParser(description="This script takes in a VCF files before any filtering on variants is done. The input VCF should contain all the sites (both variants and nonvariants). This script outputs the value of QD for each site. ")

	parser.add_argument(
            "--VCF", required=True,
            help="REQUIRED. Path to the VCF file. Should be gzipped.")

	parser.add_argument(
            "--outfile", required=True,
            help="REQUIRED. Path to the output file.")

	args = parser.parse_args()
	return args

def main():
	args = parse_args()

	outfile = open(args.outfile, "w")
	scaff=str(args.scaffold)

	with gzip.open(args.VCF, "r") as VCF:
		for line in VCF:
			if not line.startswith("#"):
				line = line.rstrip("\n")
				line = line.split("\t")
				scaffold= line[0]
				info_col = line[7]
				if info_col != ".":
					info_col_sep = info_col.split(";")
					for i in info_col_sep:
						i_sep = i.split("=")
						if i_sep[0] == "QD":
							print >>outfile, i_sep[1]

sys.exit(main())

#sys.exit()
