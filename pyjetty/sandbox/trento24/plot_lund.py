#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ast  # For converting string representations of lists/dictionaries to actual lists/dictionaries if necessary

import seaborn as sns
import argparse
import os
import yaml


import re
def split_but_keep_quotes(s):
    # Regular expression pattern to match substrings in quotes or non-whitespace sequences
    pattern = r'("[^"]*"|\'[^\']*\'|\S+)'
    matches = re.findall(pattern, s)
    return matches


def flavor_select(args, row):
	if args.flavor is not None:
		if args.flavor == 'quark':
			return row['quark'] is True
		elif args.flavor == 'gluon':
			return row['gluon'] is True
	return True


class LundSelector(object):
	def __init__(self, args, lunds):
		self.args = args
		self.ktmax = None
		if args.max_kt is True:
			# Extract the kt, Delta, and pt values
			_kts = [e['kt'] for e in lunds]
			if len(_kts) > 0:
				self.ktmax = max(_kts)
			else:
				self.ktmax = None
		self.z_cut_sd_selected = -1
		if args.z_cut_sd is not None:
			for e in lunds:
				if e['z'] > args.z_cut_sd:
					self.z_cut_sd_selected = e['z']
					break
		self.lunds = []
		for lund in lunds:
			if self.accept(lund):
				self.lunds.append(lund)

	def accept(self, lund):
		if self.ktmax is not None:
			if lund['kt'] != self.ktmax:
				return False
		if self.args.nsplit is not None:
			if lund['nsplit'] != self.args.nsplit:
				return False
		if self.args.nsplit_min is not None:
			if lund['nsplit'] < self.args.nsplit_min:
				return False
		if self.args.nsplit_max is not None:
			if lund['nsplit'] > self.args.nsplit_max:
				return False
		if self.args.kt_min is not None:
			if lund['kt'] < self.args.kt_min:
				return False
		if self.args.kt_max is not None:
			if lund['kt'] > self.args.kt_max:
				return False
		if self.args.log_kt_min is not None:
			if np.log(lund['kt']) < self.args.log_kt_min:
				return False
		if self.args.log_kt_max is not None:
			if np.log(lund['kt']) > self.args.log_kt_max:
				return False
		if self.args.z_cut is not None:
			if lund['z'] < self.args.z_cut:
				return False
		if self.args.z_cut_sd is not None:
			if lund['z'] != self.z_cut_sd_selected:
				return False
		return True

class EECplot(object):
	def __init__(self, args):
		self.args = args
		self.make_plots()
  
	def make_plots(self):
		# read yaml file defined by args.input
		with open(self.args.input, 'r') as stream:
			try:
				yaml_data = yaml.safe_load(stream)
			except yaml.YAMLError as exc:
				print(exc)
		# loop over the yaml file and extract the data
		for plot, specs in yaml_data.items():
			print('PLOT:', plot, 'SPECS:', specs)
			fig_title = self.args.title
			plot_data = {}
			plot_titles = {}
			for hist_name, hist_specs in specs.items():
				parser = get_parser()
				print('- sub item', hist_name, ':', hist_specs)
				if hist_name == 'title':
					fig_title = hist_specs
					continue
				args = parser.parse_args(split_but_keep_quotes('{}'.format(hist_specs)))
				plot_data[hist_name] = self.extract_data(args)
				plot_titles[hist_name] = args.title

			# make the plot
			# Plot the histogram using Seaborn
			plt.figure(figsize=(10, 7))
			for k, d in plot_data.items():
				if args.kde_only:
					sns.kdeplot(data=d, x='log_rl_values', weights='weights_for_log_rl', fill=args.fill, label=plot_titles[k] + ' KDE')
				else:
					sns.histplot(data=d, x='log_rl_values', weights='weights_for_log_rl', bins=75, 
                 			stat='density', element='step', fill=args.fill,
                  		kde=((not d.empty) and args.kde), label=plot_titles[k])
			plt.xlabel(r'$\ln(R_{L})$')
			if args.kde_only:
				plt.ylabel('KDE EECs')
			else:
				plt.ylabel('Normalized EECs')
			plt.title(fig_title)

			plt.legend(title=fig_title, loc='upper left', ncol=1, frameon=False)

			if args.logRL:
				plt.xlim( (np.log(0.005), np.log(1)) )
				# Create a secondary x-axis
				ax1 = plt.gca()  # Get the current axes
				# Create the secondary x-axis
				ax2 = ax1.twiny()
				ax2.set_xlabel('RL (linear scale)')
				# Set the limits and ticks on the secondary axis to correspond to the primary log scale
				log_rl_min, log_rl_max = ax1.get_xlim()
				rl_ticks = np.linspace(np.exp(log_rl_min), np.exp(log_rl_max), num=4)  # Calculate linearly spaced RL values
				rl_ticks = np.array([0.01, 0.1, 0.2, 0.3, 0.4])
				log_rl_ticks = np.log(rl_ticks)  # Convert RL values back to their log scale for correct positioning
				# Update secondary axis ticks and labels
				ax2.set_xlim(log_rl_min, log_rl_max)
				ax2.set_xticks(log_rl_ticks)
				ax2.set_xticklabels([f'{tick:.2g}' for tick in rl_ticks])  # Format RL tick labels
			else:
				plt.xlim( 0.001, 1)

			plt.show()

	def extract_data(self, args):
		rl_values = []
		# Add your indented block of code here
		weights_values = []
		n_jets = 0
		eec_types = sorted(list(dict.fromkeys([int(x) for x in args.eec_types if f'{x}'.isdigit() and int(x) >= 0 and int(x) <= 2])))
		df = pd.read_parquet(args.input)
		# Extract kt, Delta (assuming it's there), and pt for each row
		for index, row in df.iterrows():
			if args.jet_ptmin < row['pt'] < args.jet_ptmax:  # Assuming 'pt' is directly accessible in the DataFrame
				if flavor_select(args, row) is False:
					continue
				if len(row['lunds']) == 0:
					print('[w] no lund splits?', len(row['lunds']))
					continue
				n_jets += 1
				if not args.lund_eec:
					rl_values.extend(row['eecs']['RL'])
					weights_values.extend(row['eecs']['weights'])
					continue
				lselect = LundSelector(args, row['lunds'])
				# Extract the kt, Delta, and pt values
				for lund in lselect.lunds:  # Assuming 'lunds' now contains the actual list of dictionaries
					eecs = lund['eecs']  # Accessing 'eecs' dictionary
					rl_list = eecs['RL']  # Assuming 'RL' is directly accessible
					weights_list = eecs['weights']  # Assuming 'weights' is directly accessible        
					# Extend the lists with the values extracted
					#rl_values.extend(rl_list)
					#weights_values.extend(weights_list)
					_ = [rl_values.append(rl) for i,rl in enumerate(rl_list) if eecs['type'][i] in eec_types]
					#weights_values.extend(weights_list)
					_ = [weights_values.append(rl) for i,rl in enumerate(weights_list) if eecs['type'][i] in eec_types]
		# Convert the lists into a DataFrame for Seaborn
		if args.logRL:
			data_for_plot = pd.DataFrame({
					'log_rl_values': np.log([rl for rl in rl_values if rl > 0]),
					'weights_for_log_rl': [weight/n_jets for rl, weight in zip(rl_values, weights_values) if rl > 0]
			})
		else:
			data_for_plot = pd.DataFrame({
					'log_rl_values': [rl for rl in rl_values if rl > 0],
					'weights_for_log_rl': [weight/n_jets for rl, weight in zip(rl_values, weights_values) if rl > 0]
			})
		return data_for_plot


def eec_plot(args):
	df = pd.read_parquet(args.input)
	# Initialize lists to store the extracted RL and corresponding weights
	rl_values = []
	weights_values = []

	n_jets = 0
	# Extract kt, Delta (assuming it's there), and pt for each row
	for index, row in df.iterrows():
		if args.jet_ptmin < row['pt'] < args.jet_ptmax:  # Assuming 'pt' is directly accessible in the DataFrame
			if flavor_select(args, row) is False:
				continue
			if len(row['lunds']) == 0:
				print('[w] no lund splits?', len(row['lunds']))
				continue
			n_jets += 1
			lselect = LundSelector(args, row['lunds'])
			# Extract the kt, Delta, and pt values
			for lund in lselect.lunds:  # Assuming 'lunds' now contains the actual list of dictionaries
				eecs = lund['eecs']  # Accessing 'eecs' dictionary
				rl_list = eecs['RL']  # Assuming 'RL' is directly accessible
				weights_list = eecs['weights']  # Assuming 'weights' is directly accessible        
        # Extend the lists with the values extracted
        #rl_values.extend(rl_list)
        #weights_values.extend(weights_list)
				_ = [rl_values.append(rl) for i,rl in enumerate(rl_list) if eecs['type'][i] in args.eec_types]
        #weights_values.extend(weights_list)
				_ = [weights_values.append(rl) for i,rl in enumerate(weights_list) if eecs['type'][i] in args.eec_types]

	## # Ensure RL values are greater than 0 before taking the log
	## log_rl_values = np.log([rl for rl in rl_values if rl > 0])
	## weights_for_log_rl = [weight/n_jets for rl, weight in zip(rl_values, weights_values) if rl > 0]

	# Convert the lists into a DataFrame for Seaborn
	data_for_plot = pd.DataFrame({
			'log_rl_values': np.log([rl for rl in rl_values if rl > 0]),
			'weights_for_log_rl': [weight/n_jets for rl, weight in zip(rl_values, weights_values) if rl > 0]
	})

	# Plot the histogram using Seaborn
	plt.figure(figsize=(10, 7))
	sns.histplot(data=data_for_plot, x='log_rl_values', weights='weights_for_log_rl', bins=50, kde=True)
	plt.xlabel('ln(RL)')
	plt.ylabel('Normalized Weights')
	plt.title(args.title)

	# Create a secondary x-axis
	ax1 = plt.gca()  # Get the current axes
	# Create the secondary x-axis
	ax2 = ax1.twiny()
	ax2.set_xlabel('RL (linear scale)')

	# Set the limits and ticks on the secondary axis to correspond to the primary log scale
	log_rl_min, log_rl_max = ax1.get_xlim()
	rl_ticks = np.linspace(np.exp(log_rl_min), np.exp(log_rl_max), num=4)  # Calculate linearly spaced RL values
	rl_ticks = np.array([0.01, 0.1, 0.2, 0.3, 0.4])
	log_rl_ticks = np.log(rl_ticks)  # Convert RL values back to their log scale for correct positioning

	# Update secondary axis ticks and labels
	ax2.set_xlim(log_rl_min, log_rl_max)
	ax2.set_xticks(log_rl_ticks)
	ax2.set_xticklabels([f'{tick:.2g}' for tick in rl_ticks])  # Format RL tick labels

	plt.show()
  
def lund_plot(args):  
	df = pd.read_parquet(args.input)
	# Initialize lists to store extracted kt, Delta, and pt values
	kts = []
	deltas = []
	pts = []
 
	skt = 'kt'
	if args.use_kappa:
		skt = 'kappa'
  
	n_jets = 0
	# Extract kt, Delta (assuming it's there), and pt for each row
	for index, row in df.iterrows():
		if args.jet_ptmin < row['pt'] < args.jet_ptmax:  # Assuming 'pt' is directly accessible in the DataFrame
			if flavor_select(args, row) is False:
				continue
			if len(row['lunds']) == 0:
				print('[w] no lund splits?', len(row['lunds']))
				continue
			n_jets += 1
			lselect = LundSelector(args, row['lunds'])
			# Extract the kt, Delta, and pt values
			for lund in lselect.lunds:  # Assuming 'lunds' now contains the actual list of dictionaries
				kts.append(lund[skt])
				deltas.append(lund['delta'])  # Assuming each lund has a 'Delta'
				pts.append(row['pt'])  # This will repeat pt values, one for each kt, Delta pair

	print(f'found {n_jets} jets with {len(kts)} kts and {len(deltas)} deltas and {len(pts)} pts')
	# Convert lists to a DataFrame for easier manipulation
	extracted_df = pd.DataFrame({'kt': kts, 'delta': deltas, 'pt': pts})

	# Calculate the logarithmic values
	extracted_df['log_kt'] = np.log(extracted_df['kt'])
	extracted_df['log_inv_delta'] = np.log(1 / extracted_df['delta'])
	weights = np.ones_like(extracted_df['log_kt']) / n_jets

	# Set the aesthetic style of the plots
	sns.set_style("whitegrid")

	# Plotting the 2D histogram (heat map) using seaborn
	plt.figure(figsize=(10, 7))
	sns.histplot(data=extracted_df, x='log_inv_delta', y='log_kt', weights=weights, bins=75, cbar=True, cbar_kws={'label': r'$n_{\text{splits}}$ per jet'})
	plt.xlabel(r'$\ln(1/\Delta)$')
	plt.ylabel(r'$\ln(k_{t})$')
	if args.use_kappa:
		plt.ylabel(r'$\ln(\kappa)$')
	plt.xlim(np.log(1/0.4), 6)  # Limiting x-axis
	plt.ylim(args.plot_log_kt_min, args.plot_log_kt_max)  # Limiting y-axis
	if args.use_kappa:
		plt.ylim(args.plot_log_kappa_min, args.plot_log_kappa_max)  # Limiting y-axis
	# plt.title('Gluons for 100 < $p_{T}$ < 120')
	plt.title(args.title)

	# Create a secondary x-axis
	ax1 = plt.gca()  # Get the current axes
	# Create the secondary x-axis
	ax2 = ax1.twiny()
	# ax2.set_xlabel('RL (linear scale)')
	# Set the limits and ticks on the secondary axis to correspond to the primary log scale
	log_rl_min, log_rl_max = ax1.get_xlim()
	print(log_rl_min, log_rl_max)
	print(1./np.exp(log_rl_min), 1./np.exp(log_rl_max))
	rl_ticks = np.linspace(1./np.exp(log_rl_min), 1./np.exp(log_rl_max), num=5)  # Calculate linearly spaced RL values
	rl_ticks = np.array([0.4, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005])
	# print(rl_ticks)
	log_rl_ticks = np.log(1/rl_ticks)  # Convert RL values back to their log scale for correct positioning

	# Update secondary axis ticks and labels
	ax2.set_xlim(log_rl_min, log_rl_max)
	ax2.set_xticks(log_rl_ticks)
	ax2.set_xticklabels([f'{tick:.3g}' for tick in rl_ticks])  # Format RL tick labels

	return plt

def get_parser():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	parser.add_argument('-i','--input', help='input filename parquet or yaml', default='', type=str, required=True)
	parser.add_argument('-v', '--verbose', help="be verbose", default=False, action='store_true')
	parser.add_argument('-o','--output', help='output filename', default=None, type=str)
	parser.add_argument('--flavor', help='jet flavor', default=None, type=str)
	parser.add_argument('--title', help='title for figure', default=None, type=str)
	parser.add_argument('--nsplit', help='which splitting', default=None, type=int)
	parser.add_argument('--nsplit-min', help='which splitting', default=None, type=int)
	parser.add_argument('--nsplit-max', help='which splitting', default=None, type=int)
	parser.add_argument('--kt-min', help='which kt splitting', default=None, type=float)
	parser.add_argument('--kt-max', help='which kt splitting', default=None, type=float)
	parser.add_argument('--max-kt', help='use only max-kt split', default=False, action='store_true')
	parser.add_argument('--z-cut', help='cut on lund z', default=None, type=float)
	parser.add_argument('--z-cut-sd', help='cut on lund z - take the first one - soft drop', default=None, type=float)
	parser.add_argument('--log-kt-min', help='which log(kt) splitting', default=None, type=float)
	parser.add_argument('--log-kt-max', help='which log(kt) splitting', default=None, type=float)
	parser.add_argument('--plot-log-kt-min', help='plotting limit log(kt) splitting', default=-4, type=float)
	parser.add_argument('--plot-log-kt-max', help='plotting limit log(kt) splitting', default=None, type=float)
	parser.add_argument('--plot-log-kappa-min', help='plotting limit log(kt) splitting', default=-8, type=float)
	parser.add_argument('--plot-log-kappa-max', help='plotting limit log(kt) splitting', default=-1, type=float)
	parser.add_argument('--jet-ptmin', help='minimum jet pt', default=20, type=float)
	parser.add_argument('--jet-ptmax', help='maximum jet pt', default=1e5, type=float)
	parser.add_argument('--jet-etamax', help='maximum jet eta', default=2.0, type=float)
	parser.add_argument('-q', '--quit', help="just draw", default=False, action='store_true')
	parser.add_argument('-K', '--use-kappa', help="use kappa for drawing not kt", default=False, action='store_true')
	parser.add_argument('--eec-types', help='what eec types 0,1,2', default=[0, 1, 2], type=list)
	parser.add_argument('--kde-only', help="just draw kde", default=False, action='store_true')
	parser.add_argument('--kde', help="just add kde", default=False, action='store_true')
	parser.add_argument('--fill', help="use a fill", default=False, action='store_true')
	parser.add_argument('--lund-eec', help="eecs in lund splittings", default=False, action='store_true')
	parser.add_argument('--logRL', help="use log axis for EECs", default=True, action='store_true')
	return parser

def main():
	parser = get_parser()
	args = parser.parse_args()
 
	_eec_types = sorted(list(dict.fromkeys([int(x) for x in args.eec_types if f'{x}'.isdigit() and int(x) >= 0 and int(x) <= 2])))
	args.eec_types = _eec_types

	args_default = parser.parse_args('-i {args.input}'.split())
	if args.title is None:
		stitle = []
	for k, val in vars(args).items():
		if k == 'quit' or k == 'output':
			continue
		if getattr(args, k) != getattr(args_default, k):
			print(f'[i] {k}: {getattr(args, k)}')
			stitle.append(f'{k}: {getattr(args, k)}')
		else:
			print(f'[i] {k}: (default) {getattr(args, k)}')

	if args.title is None:
		args.title = 'lund plane ' + ' '.join(stitle)

	# plt = lund_plot(args)
	# if args.quit is True:
	# 	if args.output is None:
	# 		args.output = args.title.replace(' ', '_').replace(':', '_').replace('(', '_').replace(')', '_') + '.png'
	# 	plt.savefig(args.output, bbox_inches='tight', dpi=300)
	# 	plt.close()
	# else:
	# 	plt.show()

	# eec_plot(args)
	eec_plot = EECplot(args)
 
if __name__ == '__main__':
	main()