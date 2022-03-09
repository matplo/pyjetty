#!/usr/bin/env python

from __future__ import print_function

import logging
from tempfile import mkstemp
import os
import sys
import configparser
import argparse


# https://stackoverflow.com/questions/23662280/how-to-log-the-contents-of-a-configparser
class ConfigLogger(object):
    def __init__(self, log):
        self.__log = log

    def __call__(self, config):
        # self.__log.info("Config:")
        config.write(self)

    def write(self, data):
        line = data.strip()
        self.__log.info(line)


class RivetPlotDescription(object):
	def __init__(self, list_of_strings=[]):
		self.config = configparser.ConfigParser()
		self.tmp_fd, self.tmp_path = mkstemp()
		with open(self.tmp_path, 'w') as f:
			f.write('[RivetPlot]\n')
			for l in list_of_strings:
				if '=' not in l:
					continue
				if '\n' in l:
					f.write(l)
				else:
					l = l + '\n'
		self.config.read(self.tmp_path)
		os.close(self.tmp_fd)
		# print('[w] using', self.tmp_path)


class RivetHistogramDescription(object):
	def __init__(self, list_of_strings=[], number=None):
		self.config = configparser.ConfigParser()
		self.tmp_fd, self.tmp_path = mkstemp()
		self.section_name = 'Histogram'
		if number is not None:
			self.section_name += '{}'.format(number)
		with open(self.tmp_path, 'w') as f:
			f.write('[{}]\n'.format(self.section_name))
			for l in list_of_strings:
				if '=' not in l:
					continue
				if '\n' in l:
					f.write(l)
				else:
					l = l + '\n'
		self.cols = []
		read_data = False
		data = []
		for l in list_of_strings:
			if '=' in l:
				read_data = False
				continue
			if l[0:3] == '# x':
				_l = l.replace('# ', '')
				self.cols = [s.strip() for s in _l.split('\t')]
				read_data = True
				continue
			if '#' in l:
				read_data = False
				continue
			if read_data:
				data_line = [s.strip() for s in l.split('\t')]
				data.append(data_line)
		self.table = {}
		for ic, c in enumerate(self.cols):
			self.table[c] = [float(d[ic]) for d in data]
			#print(self.cols)
			#print(self.table)
		self.config.read(self.tmp_path)
		# self.config[self.section_name].add_section('data')
		# self.config[self.section_name]['data'] = self.table
		# self.config['data'] = self.table
		# self.config.add_section(self.table)
		for ic, c in enumerate(self.cols):
			#self.config[self.section_name][c] = str([float(d[ic]) for d in data])
			self.config[self.section_name][c] = ', '.join([d[ic] for d in data])
		os.close(self.tmp_fd)
		# print('[w] using', self.tmp_path)


class RivetPlot(object):
	def __init__(self, fname):
		self.fname = fname
		self.plot_section = self.read_plots()
		self.plot_config = RivetPlotDescription(self.plot_section[0])
		# print(self.plot_config.config)
		self.histo_sections = self.read_histograms()
		self.histograms = []
		for ih, hs in enumerate(self.histo_sections):
			hscfg = RivetHistogramDescription(hs, ih)
			self.histograms.append(hscfg)
		self.config_logger = None
		logging.basicConfig(stream=sys.stdout, level=logging.INFO)
		self.config_logger = ConfigLogger(logging)
		self.tmp_fd, self.tmp_path = mkstemp()
		self.foutname = self.tmp_path
		self.write_config_file(self.foutname, silent=True)
		self.config = configparser.ConfigParser()
		self.config.read(self.foutname)

	def dump(self):
		self.config_logger(self.plot_config.config)
		for h in self.histograms:
			self.config_logger(h.config)
	
	def write_config_file(self, foutname=None, silent=False):
		self.foutname = foutname
		if self.foutname is None:
			self.foutname = self.fname + '.cfg'
		with open(self.foutname, 'w') as config_file:
			self.plot_config.config.write(config_file)
			for h in self.histograms:
				h.config.write(config_file)
		if silent is False:
			logging.info('written config file ' + self.foutname)
    
	def read_sections(self, begin_string='# BEGIN PLOT', end_string='# END PLOT'):
		with open(self.fname) as f:
			lines = f.readlines()
		_read = False
		splots = []
		_spl = None
		for l in lines:
			iplot = len(splots)
			if begin_string in l:
				_descr = True
				_spl = []
				continue
			if end_string in l:
				if _descr:
					_read = True
					splots.append(_spl)
				_descr = False
				_spl = None
				continue 
			if _spl is not None:
				_spl.append(l)
		return splots

	def read_plots(self):
		return self.read_sections(begin_string='# BEGIN PLOT', end_string='# END PLOT')

	def read_histograms(self):
		return self.read_sections(begin_string='# BEGIN HISTO', end_string='# END HISTO')


def test():
	dir = '/Users/ploskon/tmp/angularities/fromJames/rivet-plots-10/ALICE_2021_I1891385'
	fname = 'd01-x01-y01.dat'
	dfname = os.path.join(dir, fname)

	rp = RivetPlot(dfname)
	# print(rp.plot_section)
	# print(rp.histo_sections)
	# rp.dump()
	rp.write_config_file()

	read_back_config = configparser.ConfigParser()
	read_back_config.read(rp.foutname)
	# read_back_config.write(sys.stdout)
	# read_back_config.write(ConfigLogger(logging))

# configparser sucks compared to ConfigObj..
def get_line_as_list_of_floats(sline):
	return [float(s) for s in sline.split(',')]

def write_graphs(fname):
	import ROOT
	from array import array 
	rp = RivetPlot(fname)
	rp.write_config_file()
	fout = ROOT.TFile(fname + '.root', 'recreate')
	for s in rp.config.sections():
		if 'Histogram' in s:
			y = array('f', get_line_as_list_of_floats(rp.config[s]['val']))
			n = len(y)
			xl 	= array('f', get_line_as_list_of_floats(rp.config[s]['xlow']))
			xh 	= array('f', get_line_as_list_of_floats(rp.config[s]['xhigh']))
			eyl = array('f', get_line_as_list_of_floats(rp.config[s]['errminus']))
			eyh = array('f', get_line_as_list_of_floats(rp.config[s]['errplus']))
			x 	= array('f', [xl[i] + (xh[i] - xl[i])/2. for i in range(n)])
			exl = array('f', [(xh[i] - xl[i])/2. for i in range(n)])
			exh = array('f', [(xh[i] - xl[i])/2. for i in range(n)])
			gr 	= ROOT.TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
			gr.SetName(rp.config[s]['title'].replace(' ', '_'))
			gr.SetTitle(rp.config[s]['title'])
			fout.cd()
			gr.Write()
	fout.Close()
	logging.info('written file {}'.format(fout.GetName()))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='project histograms', prog=os.path.basename(__file__))
	parser.add_argument('-i', '--input', help='input file',
                     default='d01-x01-y01.dat', type=str, required=True)
	args = parser.parse_args()
	if os.path.exists(args.input):
		write_graphs(args.input)

