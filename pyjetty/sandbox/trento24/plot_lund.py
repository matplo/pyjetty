#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ast  # For converting string representations of lists/dictionaries to actual lists/dictionaries if necessary

import seaborn as sns
import argparse
import os


def flavor_select(args, row):
	if args.flavor is not None:
		if args.flavor == 'quark':
			return row['quark'] is True
		elif args.flavor == 'gluon':
			return row['gluon'] is True
	return True


def lund_plot(args):
  
	df = pd.read_parquet(args.input)

	# Example DataFrame loading step
	# df = pd.read_csv('path_to_your_file.csv')  # You've already loaded your data

	# Assuming 'lunds' is a string representation of a list of dictionaries, convert it
	# df['lunds'] = df['lunds'].apply(ast.literal_eval)  # Only if 'lunds' is a string that needs conversion

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
			ktmax = None
			if args.max_kt is True:
				# Extract the kt, Delta, and pt values
				_kts = [e['kt'] for e in row['lunds']]
				if len(_kts) > 0:
					ktmax = max(_kts)
				else:
					ktmax = None
			z_cut_sd_selected = -1
			if args.z_cut_sd is not None:
				for e in row['lunds']:
					if e['z'] > args.z_cut_sd:
						z_cut_sd_selected = e['z']
			# Extract the kt, Delta, and pt values
			for entry in row['lunds']:  # Assuming 'lunds' now contains the actual list of dictionaries
				if args.max_kt:
					if entry['kt'] != ktmax:
						continue
				if args.nsplit is not None:
					if entry['nsplit'] != args.nsplit:
						continue
				if args.nsplit_min is not None:
					if entry['nsplit'] < args.nsplit_min:
						continue
				if args.nsplit_max is not None:
					if entry['nsplit'] > args.nsplit_max:
						continue
				if args.kt_min is not None:
					if entry['kt'] < args.kt_min:
						continue
				if args.kt_max is not None:
					if entry['kt'] > args.kt_max:
						continue
				if args.log_kt_min is not None:	
					if np.log(entry['kt']) < args.log_kt_min:
						continue
				if args.log_kt_max is not None:
					if np.log(entry['kt']) > args.log_kt_max:
						continue
				if args.z_cut is not None:
					if entry['z'] < args.z_cut:
						continue
				if args.z_cut_sd is not None:
					if entry['z'] != z_cut_sd_selected:
						continue
				# print(entry['kt'], ktmax, len(kts), len(deltas), len(pts))
				kts.append(entry[skt])
				deltas.append(entry['delta'])  # Assuming each entry has a 'Delta'
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
	plt.xlabel(r'$\log(1/\Delta)$')
	plt.ylabel(r'$\log(k_{t})$')
	if args.use_kappa:
		plt.ylabel(r'$\log(\kappa)$')
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


def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	parser.add_argument('-i','--input', help='input filename', default='', type=str, required=True)
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
 
	args = parser.parse_args()

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

	plt = lund_plot(args)
	if args.quit is True:
		if args.output is None:
			args.output = args.title.replace(' ', '_').replace(':', '_').replace('(', '_').replace(')', '_') + '.png'
		plt.savefig(args.output, bbox_inches='tight', dpi=300)
		plt.close()
	else:
	  plt.show()

if __name__ == '__main__':
	main()