#!/usr/bin/env python

import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import os
import argparse
import pyjetty.rootutils as rootutils


def main():
	parser = argparse.ArgumentParser(description='list branches in a tree', prog=os.path.basename(__file__))
	parser.add_argument('-t', '--tree-path', help='path to tree fname/subfolder/tree', default=None, type=str, required=True)
	parser.add_argument('-c', '--column', help='print a column not a list', default=False, action='store_true')
	parser.add_argument('-s', '--select', help='select branches comma separated on if <selection> in <branchname>', default=None, type=str)
	args = parser.parse_args()	
	with rootutils.QuietWarning():
		tu = rootutils.TreeUtils(args.tree_path)
		bnames = tu.bnames
	if args.select:
		bnames_to_print = []
		for b in bnames:
			_sels = args.select.split(',')
			for _s in _sels:
				if _s in b:
					bnames_to_print.append(b)
	else:
		bnames_to_print = bnames

	if args.column:
		for b in bnames_to_print:
			print(b)
	else:
		print(bnames_to_print)


if __name__ == '__main__':
	main()