import os
import sys
import argparse
import csv
import gzip
# dummy file for testing:
#filepath="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering/sandbox/dummyVCF.forSandbox.allSites_5_passingFilters.vcf.gz"
def parse_args():
    """
    Parse command-line arguments
    """
    parser = argparse.ArgumentParser(description="This script takes in a VCF files before any filtering on variants is done. The input VCF should contain all the sites (both variants and nonvariants). This script outputs the value of QD for each site. ")

    parser.add_argument(
            "--VCF", required=True,
            help="REQUIRED. Path to the VCF file. Should be gzipped.")

    parser.add_argument(
            "--scaffold", required=True,
            help="REQUIRED. Full name of scaffold")

    parser.add_argument(
            "--outfile", required=True,
            help="REQUIRED. Path to the output file.")

    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    outfile = open(args.outfile, "w")
    scaff=str(args.scaffold)
    outfile.write("QUAL\tDP\tQD\n") 
    with gzip.open(args.VCF, "r") as VCF:
    #with gzip.open(filepath,"r") as VCF:
        for line in VCF:
            if not line.startswith("#"):
                line = line.rstrip("\n")
                line = line.split("\t")
                scaffold= line[0]
                info_col = line[7]
                # if the scaffold isn't the scaff you've chosen, skip rest of loop
                if scaffold!=scaff:
                   continue
                # only work with lines that contain QD (snps)         
                #elif info_col != "." # don't need this because if it contains QD it won't be "." -- then you split that line
                if "QD" in info_col: 
                    myqual=line[5]
                    myinfo=line[7]
                    # split info fields:
                    # instead of iterating through each one, make a dict.
                    infoFields=dict(s.split('=') for s in myinfo.split(";"))
                    myQD=infoFields["QD"]
                    myDP=infoFields["DP"]
                    outfile.write("\t".join(str(x) for x in (myqual, myDP, myQD)))
                    outfile.write("\n")
    VCF.close()
    outfile.close()
sys.exit(main())

