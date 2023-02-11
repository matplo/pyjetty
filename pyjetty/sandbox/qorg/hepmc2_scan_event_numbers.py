from __future__ import print_function
import tqdm
import argparse
import os
import sys
import pyhepmc_ng
import threading
import multiprocessing
import subprocess
import shlex
import re


def count_events(fname):
	count = None
	try:
		with open(fname+'.count', 'r') as f:
			_s = f.readlines()[-1]
		count = int(_s)
	except:
		pass
	if count:
		return count
	input_hepmc = pyhepmc_ng.ReaderAsciiHepMC2(fname)
	if input_hepmc.failed():
		print ("[error] unable to read from {}".format(fname))
		sys.exit(1)
	event_hepmc = pyhepmc_ng.GenEvent()
	pbar = tqdm.tqdm(desc='counting events... {}'.format(fname), leave=False)
	while not input_hepmc.failed():
		_ = input_hepmc.read_event(event_hepmc)
		if input_hepmc.failed():
			pbar.update(0)
			break
		pbar.update()
	n = pbar.n
	pbar.close()
	# print('\n[i] number of events', n, 'in', fname)
	try:
	    with open(fname+'.count', 'w') as f:
        	f.writelines(['{}'.format(n)])
	except:
		pass
	return n

def is_subscriptable(o):
	try:
		_ = o[0]
	except TypeError as te:
		return False
	return True


def exec_cmnd(cmnd, shell=False):
	_args = shlex.split(cmnd)
	try:
		p = subprocess.Popen(_args, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = p.communicate()
		rc = p.returncode
	except OSError as e:
		out = f'[e] failed to execute: f{_args}'
		if is_subscriptable(e):
			err = '- Error #{0} : {1}'.format(e[0], e[1])
		else:
			err = f'- Error {e}'
		rc = 255
	return out, err, rc


def count_events_with_grep(fname):
	cmnd = f"grep -e 'E\ [1-9]+' -f {fname}"
	print('[i] trying', cmnd)
	out, err, rc = exec_cmnd(cmnd, shell=True)
	print (fname, 'out', out)
	print (fname, 'err', err)
	print (fname, 'rc', rc)


def count_events_with_re(fname):
	if os.path.isfile(fname+'.count'):
		return
	regex = r"E\ [1-9]+"
	count = 0
	pbar = tqdm.tqdm(desc=fname)
	with open(fname, "r") as f:
		while True:
			l = f.readline()
			if not l:
				break
			m = re.search(regex, l)
			if m is None:
				continue
			else:
				count = count + 1
				pbar.update(1)
	pbar.close()
	print(fname, count, file=sys.stderr)
	try:
	    with open(fname+'.count', 'w') as f:
        	f.writelines(['{}'.format(count)])
	except:
		pass


def count_threads_alive(threads):
	_count = len([thr for thr in threads if thr.is_alive()])
	return _count


def scan_with_threading(fnames, args):
	threads = list()
	pbar = tqdm.tqdm(fnames, desc='threads')
	threads_alive = {}
	for fname in fnames:
		# x = threading.Thread(target=count_events, args=(fname,))
		# x = threading.Thread(target=count_events_with_grep, args=(fname,))
		x = threading.Thread(target=count_events_with_re, args=(fname,))
		threads.append(x)
		x.start()
		pbar.update(1)
		while count_threads_alive(threads) >= multiprocessing.cpu_count() / 2: # disk IO is inefficient
			_ = [thr.join(0.01) for thr in threads if thr.is_alive()]
	pbar.close()


def main():
	parser = argparse.ArgumentParser(description='pythia8 in python', prog=os.path.basename(__file__))
	parser.add_argument('-i', '--input', help='hepmc file input to analyze', default=None, type=str)
	parser.add_argument('-l', '--list', help='treat input file as a file containing list of files to process', action='store_true', default=False)
	parser.add_argument('-g', '--debug', help='print some extras', default=False, action='store_true')

	args = parser.parse_args()

	if args.list:
		print('[i] opening', args.input, 'as list of files...')
		with open(args.input, 'r') as f:
			file_list = [fn.strip('\n') for fn in f.readlines()]
		scan_with_threading(file_list, args)

	else:
		fname = args.input
		analyze_file(fname, args)


if __name__ == "__main__":
	main()
