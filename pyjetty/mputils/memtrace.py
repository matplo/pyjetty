from pyjetty.mputils import MPBase, pinfo, perror
import os
import psutil
import sys
import array

import ROOT
ROOT.gROOT.SetBatch(True)

class MemTrace(MPBase):
	_instance = None
	_reset_once_done = False

	def __init__(self):
		raise RuntimeError('Call instance() instead')

	@classmethod
	def instance(cls):
		if cls._instance is None:
			pinfo('Creating new MemTrace instance')
			cls._instance = cls.__new__(cls)
		return cls._instance

	def _reset_once(self, **kwargs):
		if self._reset_once_done == False:
			self.configure_from_args(output_name='memtrace.root', fout=None)
			self.event_number = 0
			self.event_tree_name = 'mt'
			super(MemTrace, self).__init__(**kwargs)
			self._reset_once_done = True
			self.process = psutil.Process(os.getpid())
			self.trees = {}

	def reset(self, **kwargs):
		if self._reset_once_done == False:
			self._reset_once(**kwargs)
		self.configure_from_args(**kwargs)
		self.reset_output()

	def reset_output(self):
		if self.fout is None:
			self.fout = ROOT.TFile(self.output_name, 'recreate')
		else:
			self.fout.Close()
			self.fout = ROOT.TFile(self.output_name, 'recreate')

	def snapshot(self, label='mem'):
		self._reset_once()
		try:
			self.fout.cd()
		except:
			self.reset_output()
		try: 
			self.trees[label]
		except KeyError:
			self.fout.cd()
			self.trees[label] = ROOT.TNtuple(label, label, 'rss:vms')
			pinfo('MemTrace - new tuple {}'.format(label), file=sys.stderr)
		rss = self.process.memory_info().rss
		vms = self.process.memory_info().vms
		# n = self.trees[label].GetEntries()
		self.trees[label].Fill(rss, vms)
		self._write_()
		
	def write_as_graph(self, tree):
		data_x = []
		data_yrss = []
		data_yvms = []
		to_Gb = 1. / (1024.*1024*1024.)
		last_x   = -1
		last_rss = 0
		last_vms = 0
		for i,e in enumerate(tree):
			if e.rss != last_rss or e.vms != last_vms:
				data_x.append(last_x)
				data_yrss.append(last_rss)
				data_yvms.append(last_vms)
				data_x.append(i)
				data_yrss.append(e.rss * to_Gb)
				data_yvms.append(e.vms * to_Gb)
			last_x   = i
			last_rss = e.rss * to_Gb
			last_vms = e.vms * to_Gb
		data_x.append(last_x+1)
		data_yrss.append(last_rss)
		data_yvms.append(last_vms)
		data_x.append(last_x+2)
		data_yrss.append(0.0)
		data_yvms.append(0.0)

		f_x = array.array('f', data_x)
		f_rss = array.array('f', data_yrss)
		f_vms = array.array('f', data_yvms)

		gr_rss = ROOT.TGraph(len(data_x), f_x, f_rss)
		gr_rss.SetName('gr_rss')
		gr_rss.SetTitle('rss')
		gr_rss.Write()

		gr_vms = ROOT.TGraph(len(data_x), f_x, f_vms)
		gr_vms.SetName('gr_vms')
		gr_vms.SetTitle('vms')
		gr_vms.Write()

	def _write_(self):
		self.fout.cd()
		for t in self.trees:
			self.trees[t].Write()
		self.fout.Write()
		self.fout.Purge()

	def write(self, write_graphs = True):
		self.fout.cd()
		for t in self.trees:
			self.trees[t].Write()
			if write_graphs:
				self.write_as_graph(self.trees[t])
		self.fout.Write()			
		self.fout.Close()
		pinfo('MemTrace.write - write & close {}'.format(self.output_name))
		# except:
		# 	perror('MemTrace.write - unable to write & close {}'.format(self.output_name))

	def __del__(self):
		try:
			self.fout.Write()			
			self.fout.Close()
		except:
			None

		self._instance = None
