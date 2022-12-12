#!/usr/bin/env python

from argparse import (ArgumentParser, FileType)

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Convert Tandem Repeat Finder (TRF) dat file to bed format with repeat units for microsatellite genotyping')
    parser.add_argument(
        '-i', '--trf', type=str, required=True,
        help='Input file produced by Tandem Repeat Finder (TRF)')
    parser.add_argument(
        '-o', '--out', type=str, required=True,
        help='Output file containing genomic locations and repeat units of microsatellites.')

    return parser.parse_args() 

def main():
    # Parse command line arguments
    args = parse_args()
    trffile = args.trf
    outfile = args.out

    with open(outfile, 'w') as out:
        chrom = ""
        with open(trffile, 'r') as trf:
            for line in trf:
                splitline = line.split()
                if line.startswith("@"):
                    chrom = line.split()[0]
                else:
                    start = splitline[0]
                    end = splitline[1]
                    period = splitline[2]
                    count = splitline[3]
                    ssr = splitline[13]
                    out.write('\t'.join([chrom,start,end,period,count,ssr]) + '\n')

if __name__ == '__main__':
    main()
                    