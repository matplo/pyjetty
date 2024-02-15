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

def get_parser():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=None)
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
	parser.add_argument('--plot-log-kt-min', help='plotting limit log(kt) splitting', default=np.log(0.05), type=float)
	parser.add_argument('--plot-log-kt-max', help='plotting limit log(kt) splitting', default=np.log(30), type=float)
	parser.add_argument('--plot-log-kappa-min', help='plotting limit log(kt) splitting', default=np.log(0.001), type=float) #-8
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
	parser.add_argument('--logRL', help="use log axis for EECs", default=False, action='store_true')
	parser.add_argument('--logx', help="use log axis for EEC x-axis", default=False, action='store_true')
	parser.add_argument('--lund-plot', help="create a lund plot along eec", default=False, action='store_true')
	parser.add_argument('--abspid', help='what eec types 0,1,2', default=None, type=str)
	parser.add_argument('--n-lund-norm', help='normalize to number of lund splits', default=False, action='store_true')
	parser.add_argument('--q-adjust', help='enhance quarks by a factor - f >= 0 - default is 1.', default=1., type=float)
	parser.add_argument('--g-adjust', help='enhance gluons by a factor - f >= 0 - default is 1.', default=1., type=float)
	parser.add_argument('--palette', help='what eec types 0,1,2', default=None, type=str)
	parser.add_argument('--kde-multiple', help='for kde plot {“layer”, “stack”, “fill”}', default='layer', type=str)
	parser.add_argument('--line-styles', help='what line styles - sns :,--,-.,-', default='-', type=str)
	parser.add_argument('--line-colors', help='what line styles - b : blue.; g : green.; r : red.; c : cyan.; m : magenta.; y : yellow.; k : black.; w : white.', default='k', type=str)
	parser.add_argument('--rlxpt', help='scale RL with pT of the object', default=False, action='store_true')
	return parser

def args_from_sconfig(sconfig, dummy_input=False):
	parser = get_parser()
	if dummy_input:
		if '-i ' not in sconfig and '--input' not in sconfig:
			sconfig = '-i dummy.parquet ' + sconfig
	args = parser.parse_args(split_but_keep_quotes(sconfig))
 
	_eec_types = sorted(list(dict.fromkeys([int(x) for x in args.eec_types if f'{x}'.isdigit() and int(x) >= 0 and int(x) <= 2])))
	args.eec_types = _eec_types

	if args.abspid is not None:
		_abspid = sorted(list(dict.fromkeys([abs(int(x)) for x in args.abspid.split(',') if f'{x}'.isdigit()])))
		args.abspid = _abspid

	args_default = parser.parse_args('-i {args.input}'.split())

	stitle = []
	for k, val in vars(args).items():
		if k == 'quit' or k == 'output':
			continue
		if getattr(args, k) != getattr(args_default, k):
			stitle.append(f'{k}: {getattr(args, k)}')

	if args.title is None:
		args.title = 'lund plane ' + ' '.join(stitle)

	schanged = []
	sdefaults = []
	print('[i] settings summary:')
	for k, val in vars(args).items():
		if k == 'quit' or k == 'output':
			continue
		if getattr(args, k) != getattr(args_default, k):
			schanged.append(f'    {k}: {getattr(args, k)}')
		else:
			sdefaults.append(f'    {k}: (default) {getattr(args, k)}')
	if args.verbose:
		for s in schanged:
			print(s)
		for s in sdefaults:
			print(s)

	return args, args_default

def accept_pt(args, row):
	return (args.jet_ptmin < row['pt'] < args.jet_ptmax)

def accept_abspid(args, row):
	if args.abspid is None:
		return True
	return (row['abspid'] in args.abspid)

import itertools
def get_cycle_n(l, n):
	cc = itertools.cycle(l)
	rv = []
	for i in range(n):
		rv.append(next(cc))
	return rv

def get_line_styles(args, n):
	lstyles = ['-', '--', '-.', ':']
	# rv = ['-'] * n
	rv = get_cycle_n(lstyles, n)
	i = 0
	for s in args.line_styles.replace('"','').replace("'",'').split(','):
		if i < len(rv) and len(s)	> 0:
			rv[i] = s
			i += 1
	rv = list(reversed(rv))
	print(rv)
	return rv

def get_line_colors(args, n):
	lcolors = ['k', 'b', 'g', 'r', 'c', 'm']
	rv = get_cycle_n(lcolors, n)
	print(rv)
	i = 0
	for s in args.line_colors.replace('"','').replace("'",'').split(','):
		if i < len(rv) and len(s)	== 1:
			rv[i] = s
			i += 1
	rv = list(reversed(rv))
	print(rv)
	return rv

######################### LUND_PLOT #########################

def lund_plot(args):  
	df = pd.read_parquet(args.input)
	# Initialize lists to store extracted kt, Delta, and pt values
	kts = []
	deltas = []
	pts = []
 
	skt = 'kt'
	if args.use_kappa:
		skt = 'kappa'
  
	n_jets = len([index for index, row in df.iterrows() if accept_pt(args, row)])
	# Extract kt, Delta (assuming it's there), and pt for each row
	for index, row in df.iterrows():
		if accept_pt(args, row) and accept_abspid(args, row):  # Assuming 'pt' is directly accessible in the DataFrame
			if flavor_select(args, row) is False:
				continue
			if len(row['lunds']) == 0:
				# print('[w] no lund splits?', len(row['lunds']))
				continue
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
	plt.figure(figsize=(14, 7))
	#if args.kde:
	#	# sns.kdeplot(data=extracted_df, x='log_inv_delta', y='log_kt', weights=weights, cmap="Greys", fill=True, cbar=True) # "Reds" "Blues"
	#	sns.kdeplot(data=extracted_df, x='log_inv_delta', y='log_kt', weights=weights, cmap="Blues", fill=True, cbar=True) # "Reds" "Blues"
	#else:
	#	sns.histplot(data=extracted_df, x='log_inv_delta', y='log_kt', weights=weights, bins=75, cbar=True, cbar_kws={'label': r'$n_{\text{splits}}$ per jet'})
	sns.jointplot(
        data=extracted_df, x='log_inv_delta', y='log_kt', weights=weights,
        cmap="Blues", fill=True, cbar=True,
        kind="kde", height=6
		)
	plt.xlabel(r'$\ln(1/\Delta)$')
	plt.ylabel(r'$\ln(k_{t})$')
	if args.use_kappa:
		plt.ylabel(r'$\ln(\kappa)$')
	plt.xlim(np.log(1/0.4), 5.5)  # Limiting x-axis
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
	print(1./np.exp(log_rl_min), 1./np.exp(log_rl_max))
	rl_ticks = np.linspace(1./np.exp(log_rl_min), 1./np.exp(log_rl_max), num=5)  # Calculate linearly spaced RL values
	rl_ticks = np.array([0.4, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005])
	# print(rl_ticks)
	log_rl_ticks = np.log(1/rl_ticks)  # Convert RL values back to their log scale for correct positioning
	# Update secondary axis ticks and labels
	ax2.set_xlim(log_rl_min, log_rl_max)
	ax2.set_xticks(log_rl_ticks)
	ax2.set_xticklabels([f'{tick:.3g}' for tick in rl_ticks])  # Format RL tick labels

	ax3 = ax1.twinx()
	log_rl_min, log_rl_max = ax1.get_ylim()
	rl_ticks = np.linspace(np.exp(log_rl_min), np.exp(log_rl_max), num=5)  # Calculate linearly spaced RL values
	rl_ticks = []
	if args.use_kappa:
		yticks = [0.001, 0.01, 0.05, 0.1]
	else:
		yticks = [0.1, 0.5, 1, 2, 5, 10, 20]
	for y in yticks:
		if np.log(y) > log_rl_min and np.log(y) < log_rl_max:
			rl_ticks.append(y)
	# print(rl_ticks)
	log_rl_ticks = np.log(rl_ticks)  # Convert RL values back to their log scale for correct positioning
	# Update secondary axis ticks and labels
	ax3.set_ylim(log_rl_min, log_rl_max)
	ax3.set_yticks(log_rl_ticks)
	ax3.set_yticklabels([f'{tick:.3g}' for tick in rl_ticks])  # Format RL tick labels
	if args.output:
		plt.savefig(args.output.replace('.png', '_lund.png'), dpi=300, transparent=True, bbox_inches='tight')
	plt.show()
	return plt

def flavor_adjust(args, row):
		_f = 1.
		if row['quark'] is True:
			_f = args.q_adjust
		if row['gluon'] is True:
			_f = args.g_adjust
		if _f < 0:
			_f = 0.
		return _f

def fix_quote_title(title):
	if title[0] == '"' and title[-1] == '"':
		return title[1:-1]
	return title

######################### EECPLOT #########################

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
			plot_args = {}
			common_args = self.args
			for hist_name, hist_specs in specs.items():
				print('- sub item', hist_name, ':', hist_specs)
				if hist_name == 'title':
					fig_title = hist_specs
					continue
				if hist_name == 'args':
					common_args, common_args_default = args_from_sconfig(hist_specs, dummy_input=True)
					continue
				# args = parser.parse_args(split_but_keep_quotes('{}'.format(hist_specs)))
                # sconfig = split_but_keep_quotes('{}'.format(hist_specs))
				args, args_default = args_from_sconfig(hist_specs)
				plot_data[hist_name] = self.extract_data(args, common_args)
				plot_args[hist_name] = args
			print('-', fig_title)
			# make the plot
			# Plot the histogram using Seaborn
			plt.figure(figsize=(10, 7))
			result_df = pd.concat(
                [df.assign(EECs=plot_args[k].title) for k, df in plot_data.items()],
                ignore_index=True
                )

			xmin = 0.005
			xmax = 1.0
			if args.rlxpt or common_args.rlxpt:
				xmin = 0.35
				xmax = 1e2

			# https://seaborn.pydata.org/generated/seaborn.color_palette.html#seaborn.color_palette
			pal = None
			if common_args.palette:
				pal = sns.color_palette(common_args.palette)
			if common_args.kde_only:
				self.plot = sns.kdeplot(data=result_df, x='rl', weights='ws', hue='EECs', palette=pal, fill=common_args.fill,
	                	multiple=common_args.kde_multiple)
			else:
				self.plot = sns.histplot(data=result_df, x='rl', weights='ws', hue='EECs', binrange=(np.log(xmin), np.log(xmax)),
                 		bins=50,common_bins=True,kde=((not result_df.empty) and common_args.kde),fill=common_args.fill, element='step', palette=pal,
                    log_scale=(common_args.logRL or common_args.logx))

			if args.logRL or common_args.logRL:
				plt.xlabel(r'$\ln(R_{L})$')
				if args.rlxpt or common_args.rlxpt:
					plt.xlabel(r'$\ln(R_{L} x p_{T})$')
			else:
				plt.xlabel(r'$R_{L}$')
				if args.rlxpt or common_args.rlxpt:
					plt.xlabel(r'$R_{L} ^{x} p_{T})$')

			if common_args.kde_only:
				plt.ylabel('KDE EECs')
			else:
				plt.ylabel('Normalized EECs')
			plt.title(fig_title)
   
			# https://stackoverflow.com/questions/70089199/how-to-set-a-different-linestyle-for-each-hue-group-in-a-kdeplot-displot
			if common_args.fill:
				handles = self.plot.legend_.legendHandles[::-1]
				nplots = len(self.plot.collections)
				for line, ls, lc, handle in zip(self.plot.collections, 
						get_line_styles(common_args, nplots), 
						get_line_colors(common_args, nplots),
						handles):
						line.set_linestyle(ls)
						line.set_color(lc)
						handle.set_ls(ls)
						handle.set_color(lc)
			else:
				handles = self.plot.legend_.legendHandles[::-1]
				nplots = len(self.plot.lines)
				for line, ls, lc, handle in zip(self.plot.lines, 
						get_line_styles(common_args, nplots), 
						get_line_colors(common_args, nplots),
						handles):
						line.set_linestyle(ls)
						line.set_color(lc)
						handle.set_ls(ls)
						handle.set_color(lc)
	    # plt.legend(title=fig_title, loc='upper left', ncol=1, frameon=False)
			if args.logRL:
				plt.xlim( (np.log(xmin), np.log(xmax)) )
				if args.rlxpt or common_args.rlxpt:
					pass
				else:
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
				plt.xlim( xmin, xmax)

			if args.logx or common_args.logx:
				plt.xscale('log')
			plot_output_fname = None
			if self.args.output:
				plot_output_fname = self.args.output.split('.')[0] + '_' + fig_title + '.png'
				plt.savefig(plot_output_fname, dpi=300, transparent=True, bbox_inches='tight')

			plt.show()
			ilund = 0
			for k, pargs in plot_args.items():
				ilund += 1
				if pargs.lund_plot:
					pargs.output = plot_output_fname.split('.png')[0] + '{}.png'.format(ilund)
					plt_lund = lund_plot(pargs)

	########################################### EXTRACT EEC DATA
	def extract_data(self, args, common_args):
		rl_values = []
		# Add your indented block of code here
		weights_values = []
		eec_types = sorted(list(dict.fromkeys([int(x) for x in args.eec_types if f'{x}'.isdigit() and int(x) >= 0 and int(x) <= 2])))
		df = pd.read_parquet(args.input)
		# Extract kt, Delta (assuming it's there), and pt for each row
		n_jets = len([index for index, row in df.iterrows() if accept_pt(args, row)])
		n_lunds = 0
		# print('- extract data with', args)
		print('- number of jets', n_jets)
		for index, row in df.iterrows():
			if accept_pt(args, row) and accept_abspid(args, row):
				if flavor_select(args, row) is False:
					continue
				quark_glue_adjust = flavor_adjust(args, row)
				# print(row['quark'], row['gluon'], quark_glue_adjust)
				if len(row['lunds']) == 0:
					# print('[w] no lund splits?', len(row['lunds']))
					continue
				rlpTscale = 1.
				if args.rlxpt or common_args.rlxpt:
					rlpTscale = row['pt']
				if not args.lund_eec:
					# rl_values.extend(row['eecs']['RL'])
					_ = [rl_values.append(rl * rlpTscale) for rl in row['eecs']['RL']] # if eecs['type'][i] in eec_types] type -1 here
					# weights_values.extend(row['eecs']['weights'])
					_ = [weights_values.append(quark_glue_adjust * w) for w in row['eecs']['weights']] # if eecs['type'][i] in eec_types]
					continue
				lselect = LundSelector(args, row['lunds'])
				n_lunds += len(lselect.lunds)
				# Extract the kt, Delta, and pt values
				for lund in lselect.lunds:  # Assuming 'lunds' now contains the actual list of dictionaries
					eecs = lund['eecs']  # Accessing 'eecs' dictionary
					rl_list = eecs['RL']  # Assuming 'RL' is directly accessible
					weights_list = eecs['weights']  # Assuming 'weights' is directly accessible
					rlpTscale = 1.
					if args.rlxpt or common_args.rlxpt:
						rlpTscale = lund['pt']
					# Extend the lists with the values extracted
					#rl_values.extend(rl_list)
					#weights_values.extend(weights_list)
					_ = [rl_values.append(rl * rlpTscale) for i,rl in enumerate(rl_list) if eecs['type'][i] in eec_types]
					#weights_values.extend(weights_list)
					_ = [weights_values.append(quark_glue_adjust * w) for i,w in enumerate(weights_list) if eecs['type'][i] in eec_types]
		norm = n_jets
		if args.n_lund_norm:
			norm = n_lunds
    # Convert the lists into a DataFrame for Seaborn
		if args.logRL:
			data_for_plot = pd.DataFrame({
					'rl': np.log([rl for rl in rl_values if rl > 0]),
					'ws': [weight/norm for rl, weight in zip(rl_values, weights_values) if rl > 0],
                    'EECs' : fix_quote_title(args.title)
      })
		else:
			data_for_plot = pd.DataFrame({
					'rl': [rl for rl in rl_values if rl > 0],
					'ws': [weight/norm for rl, weight in zip(rl_values, weights_values) if rl > 0],
                    'EECs' : fix_quote_title(args.title)
			})
		return data_for_plot

