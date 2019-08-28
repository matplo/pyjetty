#!/usr/bin/env python

import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import os
import argparse
import rootutils


def main():
	parser = argparse.ArgumentParser(description='list branches in a tree', prog=os.path.basename(__file__))
	parser.add_argument('-t', '--tree-path', help='path to tree fname/subfolder/tree', default=None, type=str, required=True)
	args = parser.parse_args()	
	with rootutils.QuietWarning():
		tu = rootutils.TreeUtils(args.tree_path)
		bnames = tu.bnames
	print (bnames)


if __name__ == '__main__':
	main()