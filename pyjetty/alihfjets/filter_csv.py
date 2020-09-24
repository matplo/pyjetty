#!/usr/bin/env python

import csv
import os
import argparse


def main():
	parser = argparse.ArgumentParser(description='filter csv file with D0 counts; file format is in rows: fname,count', prog=os.path.basename(__file__))
	parser.add_argument('fname', help='single root file or a file with a list of files to process', type=str, default='')
	parser.add_argument('--nmin', help='minimum number of counts in a file', type=int, default=-1, required=False)
	parser.add_argument('--show-counts', help='show the counts next to the file name', action='store_true', default=False, required=False)
	parser.add_argument('--stop-sum', help='stop when the sum reaches a value', type=int, default=0, required=False)
	parser.add_argument('--show-sum', help='show the total n in the selected files', action='store_true', default=False, required=False)
	args = parser.parse_args()

	total_so_far = 0
	with open(args.fname, newline='') as csvfile:
		reader = csv.DictReader(csvfile)
		for row in reader:
			nD0s = int(row['ND0_cand_events'])
			if  nD0s >= args.nmin:
				total_so_far = total_so_far + nD0s
				if args.show_counts:
					print(row['fname'], nD0s)
				else:
					print(row['fname'])
				if args.stop_sum > 0:
					if total_so_far > args.stop_sum:
						break
	if args.show_sum:
		print('- total', total_so_far)


if __name__ == '__main__':
	main()