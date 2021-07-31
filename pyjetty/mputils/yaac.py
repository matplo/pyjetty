#!/usr/bin/env python

import argparse
import os
import subprocess
import string
import tqdm
import yaml
import os

## this is YAAC - yet another alien copy

## # this is example configuration file to download some MC
## run yaac.py twice:
## 1) with --list to build the list of files to download
## 2) without --list to actually copy - note existing files on local disk are NOT overwritten; warning will be printed

## # Configuration to download LHC18a4a2 sim
## # (to be used by download_data.py)
## 
## period: 'LHC18a4a2_fast' # pass3
## parent_dir: 'sim'
## year: '2018'
## train_name: 'HF_TreeCreator'
## train_PWG: 'AOD235/PWGHF'
## train_number: '635_20210528-0103_child_1'
## runlist: [282016, 282021, 282025, 282031, 282050, 282051, 282078, 282098, 282099, 282118, 282120, 282122, 282123, 282125, 282126, 282127, 282146, 282147, 282189, 282206, 282224, 282227, 282229, 282230, 282247, 282302, 282303, 282304, 282305, 282306, 282307, 282309, 282312, 282313, 282314, 282340, 282341, 282342, 282343]
## 
## output_dir: '/mnt/rstorage/alice/sim/LHC18/LHC18a4a2/635/child_1'
## 
## # these are taken as some default values by yaac.py
## nthreads: 20
## file_pattern: 'AnalysisResults.root'
## pt_hat_bins: ['/']
## or could be
## pt_hat_bins: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]


def exec_cmnd(cmnd_with_args):
	vret = []
	try:
		cproc = subprocess.run(cmnd_with_args.split(' '), shell=False, check=False, stdout=subprocess.PIPE)
		if cproc.stdout:
			sds = cproc.stdout.decode("utf-8").replace('\n', 'EOL').split('EOL')
			vret = [s for s in sds if len(s) > 0]
	except:
		pass
	return vret


import threading

class AlienCopyFile(threading.Thread):
	def __init__(self, name, alien_fname, local_fname):
		threading.Thread.__init__(self)
		self.name = name
		self.alien_fname = alien_fname
		self.local_fname = local_fname
		self.stdout = None
		self.stderr = None

	def run(self):
		cmnd = 'alien_cp alien://{} file://{}'.format(self.alien_fname, self.local_fname)
		cproc = subprocess.run(cmnd.split(' '), shell=False, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if cproc.stdout:
			self.stdout = cproc.stdout.decode('utf-8')
		if cproc.stderr:
			self.stderr = cproc.stderr.decode('utf-8')
			os.remove(self.local_fname)


def str_run_number_from_file(alien_f, config):
	# guess the run number in the file as stand alone directory name
	parts = alien_f.split('/')
	sruns = [str(r) for r in config['runlist']]
	for _p in parts:
		if _p in sruns:
			return _p
	return None


def do_copy(config):
	copy_jobs = []
	copy_jobs_finished = []
	try:
		with open(config['file_list'], 'r') as flist:
			file_list = [_f.strip('\n') for _f in flist.readlines()]
	except:
		print('[e] no file list - expected in {}'.format(config['file_list']))
		return
	if len(file_list) < 1:
		print('[e] no files to copy - read from {}'.format(config['file_list']))
		return		
	print('[i] will use {} threads'.format(config['nthreads']))
	warnings = []
	skipped = []
	for alien_f in tqdm.tqdm(file_list):
		run_number = str_run_number_from_file(alien_f, config)
		if run_number is None:
			warnings.append('[w] skipped {} - missing run number?'.format(alien_f))
			continue
		local_fname = '/'.join([config['output_dir'], str(run_number), alien_f.strip('\n').split(config['train_number'])[1]])
		# print('cp {} to {}'.format(alien_f, local_fname))
		if os.path.isfile(local_fname):
			warnings.append(local_fname)
			continue
		os.makedirs(os.path.dirname(local_fname), exist_ok=True)
		cj = AlienCopyFile(name=alien_f.replace('/', '_'), alien_fname=alien_f, local_fname=local_fname)
		while len(copy_jobs) > config['nthreads']:
			for j in copy_jobs:
				j.join(1)
				if j.is_alive():
					pass
				else:
					if j.stderr:
						print('[e] copying', j.alien_fname, '\n', j.stderr)
					copy_jobs.remove(j)
		copy_jobs.append(cj)
		cj.start()
	if len(skipped) > 0:
		print('[w] files NOT copied - already existing:')
		for _f in skipped:
			print('    ', _f)
	for _w in warnings:
		print(_w)
	print('[i] copy done.')

def compile_basedir_list(config):
	d = ['/alice']
	for feature in ['parent_dir', 'year', 'period', 'pt_hat_bins[]', 'runlist[]', 'train_PWG', 'train_name', 'train_number']:
		if '[]'	in feature:
			d.append(feature)
		else:
			d.append(config[feature])
	bdir = '/'.join(d)
	print('[i] basedirs for:', bdir)

	dlist = []
	dlist.append(bdir)
	for feature in ['parent_dir', 'year', 'period', 'pt_hat_bins[]', 'runlist[]', 'train_PWG', 'train_name', 'train_number']:
		if '[]'	in feature:
			_tmp_dlist = []
			for _ifeat in config[feature.replace('[]','')]:
				for _dl in dlist:
					_d = _dl.replace(feature, '{}'.format(_ifeat))
					_tmp_dlist.append(_d)
			_ = [dlist.append(_d) for _d in _tmp_dlist]
	# cleanup
	dlist = [_d for _d in dlist if '[]' not in _d]
	return dlist
	#return '/alice/{pdir}/{year}/{period}/{pt_hat}/{run_number}/{train_PWG}/{train_}/{train_number}'.format(pdir=config['parent_dir'])

def is_in_subdir(dir, should_be_a_subdir):
	_d = os.path.dirname(dir)
	# print(should_be_a_subdir, _d, dir, (should_be_a_subdir in _d))
	return(should_be_a_subdir in _d)

def find_files(config):
	files = []
	bdirs = compile_basedir_list(config)
	for bdir in bdirs:
		_cmnd = 'alien_find {_dir} {_what}'.format(_dir=bdir, _what=config['file_pattern'])
		print('[i] executing:', _cmnd)
		_files = exec_cmnd(_cmnd)
		_ = [files.append(_f) for _f in _files if is_in_subdir(os.path.dirname(_f), config['train_number']) ]
	print('    found {} files.'.format(len(files)))
	return files

def main():
	parser = argparse.ArgumentParser(description='copy file list from alien', prog=os.path.basename(__file__))
	#parser.add_argument('input', default=None, help="list of files in a directory to be copied; alternatively could be /alice/data/2018/LHC18 with --list", type=str)
	parser.add_argument('--list', default=False, action="store_true")
	parser.add_argument('--nthreads', default=0, type=int, help='number of threads')
	parser.add_argument('-c', '--configFile', default='', 
									help='yaml config file',
									type=str, required=True)
	args = parser.parse_args()

	with open(args.configFile, 'r') as stream:
		config = yaml.safe_load(stream)

	# print('#', args)
	try:
		_fp = config['file_pattern']
	except:
		config['file_pattern'] = 'AnalysisResults.root'

	try:
		_fp = config['pt_hat_bins']
	except:
		config['pt_hat_bins'] = '/'

	try:
		_fp = config['file_list']
	except:
		config['file_list'] = args.configFile + '.FileList'

	try:
		_fp = config['nthreads']
	except:
		config['nthreads'] = 10

	if args.nthreads > 0:
		config['nthreads'] = args.nthreads

	print('config', config)

	if args.list:
		files_to_copy = find_files(config)
		with open(config['file_list'], 'w') as f:
			f.writelines([_f+'\n' for _f in files_to_copy])
		print('[i] file list written to', config['file_list'])
	else:
		print('[i] using file list', config['file_list'])
		do_copy(config)

if __name__ == '__main__':
	main()

