#!/usr/bin/env python3

from pyjetty.mputils import MPBase, perror, pinfo, pwarning, pindent, pdebug, treewriter
import pyjetty.rootutils as rootutils
import argparse
import os
import uproot
import pandas as pd
import tqdm
import ROOT
ROOT.gROOT.SetBatch(True)

class SplitDataFileIO(MPBase):
	def __init__(self, **kwargs):
		super(SplitDataFileIO, self).__init__(**kwargs)
		self.fout = None

	def get_list_of_trees(self, fname):
		with rootutils.QuietWarning():
			tu = rootutils.list_trees_dict(fname)
		return [t.decode("utf-8") for t in tu]

	def get_event_id(self, df):
		self.ev_id = df['ev_id'].values[0]
		self.run_number = df['run_number'].values[0]
		#pindent('processing', 'run_number=', self.run_number, 'ev_id=', self.ev_id, 'internal iev=', self.iev, 'cyclical iev=', self.iev_count)
		if self.iev_count >= self.max_events or self.iev == 0:
			#open a new file
			if self.fout:
				self.fout.Write()
				self.fout.Close()
				self.td = None
				self.fout = None
			outfname = self.file_output.replace('.root', '_{}.root'.format(self.nfile))
			pinfo('opening new file at', self.iev, 'file number', self.nfile, 'file name', outfname)
			self.fout = ROOT.TFile(outfname, 'recreate')
			self.nfile = self.nfile + 1
			self.iev_count = 0
			self.trees = {}

			# these are branches
			# for sdf in df:

			for sdf in self.df:
				foldername = os.path.dirname(sdf.split(';')[0])
				tname = os.path.basename(sdf.split(';')[0])
				self.fout.cd()
				if self.td is None:
					self.td = ROOT.TDirectoryFile(foldername, foldername)
				self.td.cd()
				# tw = ROOT.TNtuple(tname, tname, ':'.join(self.df[sdf].columns))
				self.trees[tname] = treewriter.RTreeWriter(tree_name=tname, fout=self.td)

		# fill for loc of iev and run number for each tree

		# for sdf in self.df:
		# 	pinfo(self.df[sdf])

		for sdf in self.df:
			tname = os.path.basename(sdf.split(';')[0])
			_df_sel = self.df[sdf].loc[(self.df[sdf]['run_number']==self.run_number) & (self.df[sdf]['ev_id']==self.ev_id)]
			#if len(_df_sel) < 1:
			#	pwarning('no entries for', sdf, len(_df_sel), self.run_number, self.ev_id)
			#else:
			#	pinfo('rows for', sdf, len(_df_sel))
			for index, row in _df_sel.iterrows():
				#print(row)
				for c in _df_sel.columns:
					# pindent(self.trees[tname])
					# pindent(tname, index, c, '=', row[c])
					val = row[c]
					self.trees[tname].fill_branch(c, row[c])
				# pindent(self.trees[tname].tree.GetName(), 'fill')
				self.trees[tname].fill_tree()

		self.iev = self.iev + 1
		self.iev_count = self.iev_count + 1

	def fill_evid_rno(self, df):
		#pindent(df)
		pass

	def split_file(self, file_input, file_output, nevents, tolerance=0.1):
		pinfo('spliting file', file_input, 'max nevents:', nevents, 'tolerance:', tolerance)
		tlist = self.get_list_of_trees(file_input)
		pinfo('list of trees:', tlist)
		self.df_grouped = {}
		self.df = {}
		for tname in tlist:
			t = uproot.open(file_input)[tname]
			self.df[tname] = t.pandas.df()
			self.df_grouped[tname] = self.df[tname].groupby(['run_number','ev_id'])
			pindent('tree', tname, 'Nrows in pandas:', len(self.df[tname]), 'Nevents:', len(self.df_grouped[tname]))
		# get max events order
		df_names_sorted = sorted(self.df_grouped, key=lambda s: len(self.df_grouped[s]), reverse=True)
		for dfname in df_names_sorted:
			pindent(dfname, len(self.df_grouped[dfname]))
		self.iev = 0
		self.iev_count = 0
		self.max_events = nevents
		self.file_output = file_output
		self.nfile = 0
		self.fout = None
		self.td = None
		self.df_grouped[df_names_sorted[0]].apply(self.get_event_id)
		# for sname in df_names_sorted:
		# 	pinfo('processing for', sname)
		# 	self.df_grouped[sname].apply(self.get_event_id)

	def __del__(self):
		if self.fout:
			self.fout.Write()
			self.fout.Close()

def main():
	parser = argparse.ArgumentParser(description='split a root file', prog=os.path.basename(__file__))
	parser.add_argument('-f', '--fname', help='path to a root file', default=None, type=str, required=True)
	parser.add_argument('-o', '--fout', help='output a root file - basename', default=None, type=str, required=True)
	parser.add_argument('-n', '--n-max-events', help='maximum number of events in a file', default=500, type=int, required=True)
	args = parser.parse_args()	
	sp = SplitDataFileIO()
	sp.split_file(args.fname, args.fout, args.n_max_events)

if __name__ == '__main__':
	main()
