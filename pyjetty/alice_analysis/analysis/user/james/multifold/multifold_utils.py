#!/usr/bin/env python

"""
  Analysis utilities for multifold analysis.
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import ctypes

# Data analysis and plotting
import numpy as np
import pandas as pd

try:
    import ROOT
except ImportError:
    pass

# Base class
from pyjetty.alice_analysis.analysis.user.substructure import analysis_utils_obs

################################################################
class MultiFoldUtils(analysis_utils_obs.AnalysisUtils_Obs):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, **kwargs):
        super(MultiFoldUtils, self).__init__(**kwargs)

    #---------------------------------------------------------------
    # Return formatted observable key
    #---------------------------------------------------------------
    def obs_key_formatted(self, observable_info, observable, i_subobservable, tlatex=False):

        s = observable_info[observable]['xtitle']

        obs_setting = observable_info[observable]['obs_settings'][i_subobservable]
        subobs_label_list = self.formatted_subobs_label(observable)
        if subobs_label_list:
            subobs_label = subobs_label_list
            s += f' {subobs_label} = {obs_setting}'

        grooming_setting = observable_info[observable]['obs_grooming_settings'][i_subobservable]
        if grooming_setting:
            s += f' {self.formatted_grooming_label(grooming_setting)}'

        # Convert from tlatex to latex
        if not tlatex:
            s = f'${s}$'
            s = s.replace('#it','')
            s = s.replace('} {','},\;{')
            s = s.replace('#','\\')
            s = s.replace('SD',',\;SD')
            s = s.replace(', {\\beta} = 0', '')
            s = s.replace('{\Delta R}','')
            s = s.replace('Standard_WTA','\mathrm{Standard-WTA}')
            s = s.replace('{\\lambda}_{{\\alpha}},\;{\\alpha} = ','\lambda_')

        return s

    #---------------------------------------------------------------
    # Mask all results according to specified conditions
    #   variables = [[variable1], [variable2, variable3], ...]
    #   cuts = [ [[variable1_min,variable1_max]], [[variable2_min,variable2_max], [variable3_min,variable3_max]], ...]
    #
    # Each entry in the variables/cuts list corresponds to a set of cuts that will be simultanously applied
    #
    # The optional argument `mask_data_type_dict` is a dict that specifies which data_type to use to determine the mask
    #   e.g. mask_data_type_dict={'sim_truth': 'sim_truth', 'sim_det': 'sim_truth'}  will determine the mask for 'sim_truth' 
    #        for both 'sim_det' and 'sim_truth'.
    #   Any keys not specified will use their own data type to determine the mask.
    #
    # The input dict 'results' is expected to be in the format 
    #   - results[data_type][obs_key] = np.array([...])
    #       where data_type is e.g.: 'data', 'mc_det_matched', 'mc_truth_matched'
    #       and obs_key is self.observable_info[observable]['obs_key'][i] (plus 'pt_hat_scale_factors' for MC)
    #
    # A dictionary will be returned, containing the different cut combinations specified, e.g.:
    #   result_dict['data'][obs_key][f'{variable1}{variable1_min}-{variable1_max}] = result1
    #   result_dict['data'][obs_key][f'{variable2}{variable2_min}-{variable2_max}_{variable3}{variable3_min}-{variable3_max}] = result2
    # The dictionary will include both data and MC and det/truth levels.
    #
    # The reason for this is that we want to support both:
    #   - Returning multiple different cuts on the same variable (e.g. pt=[20,40,60,80])
    #   - Cutting on multiple variables simultaneously
    #---------------------------------------------------------------
    def apply_cut(self, results, n_jets, variables_list, cuts_list, mask_data_type_dict=None, n_jets_max=100000000):

        if len(variables_list) != len(cuts_list):
            raise ValueError(f'variables_list has different length than cuts_list! {variables_list} vs. {cuts_list}')

        results_dict = self.recursive_defaultdict()
        
        # Loop through both data and MC
        for data_type in results.keys():

            if mask_data_type_dict and data_type in mask_data_type_dict.keys():
                if n_jets[mask_data_type] != n_jets[data_type]:
                    raise ValueError(f'Data type mask does not correspond to the appropriate number of jets: {n_jets[mask_data_type]} vs. {n_jets[data_type]}')
                mask_data_type = mask_data_type_dict[data_type]
            else:
                mask_data_type = data_type

            # Loop over all cut combinations
            for variables, cuts in zip(variables_list, cuts_list):

                # Loop through all cuts in a given combination and construct mask
                total_mask = np.zeros(n_jets[data_type], dtype=bool) # For simplicity: cut_mask[i]=True means we will discard index i
                cut_key = ''
                for i in range(len(variables)):
                    variable = variables[i]
                    cut_min, cut_max = cuts[i]
                    cut_mask = (results[mask_data_type][variable] < cut_min) | (results[mask_data_type][variable] > cut_max)
                    total_mask = (total_mask | cut_mask)

                    if i>0:
                        cut_key += '_'
                    cut_key += f'{variable}{cut_min}-{cut_max}'

                # Now that we have the complete mask for this combination, apply it to all arrays
                for obs_key,result in results[data_type].copy().items():
                    results_dict[data_type][obs_key][cut_key] = result[~total_mask][:n_jets_max]
                if 'mc' in data_type:
                    results_dict[data_type]['pt_hat_scale_factors'][cut_key] = results[data_type]['pt_hat_scale_factors'].copy()[~total_mask][:n_jets_max]
                    
        return results_dict

    #---------------------------------------------------------------
    # Construct dataframe and ndarray of observables for each data_type
    #---------------------------------------------------------------
    def convert_results_to_ndarray(self, results, observables, observable_info):
        print('  Converting data to ndarray...')

        df_dict = {}
        ndarray_dict = {}
        n_jets_dict = {}
        cols_dict = {}
        for data_type in results.keys():
            df_dict[data_type] = pd.DataFrame()
            for observable in observables:
                for i in range(observable_info[observable]['n_subobservables']):
                    obs_key = observable_info[observable]['obs_keys'][i]
                    df_dict[data_type][obs_key] = results[data_type][obs_key]

            # Convert to ndarray
            ndarray_dict[data_type] = df_dict[data_type].to_numpy()
            n_jets_dict[data_type] = ndarray_dict[data_type].shape[0]

            # Store columns
            cols_dict[data_type] = list(df_dict[data_type].columns)

        # Check that the columns are ordered the same for all data_types
        for i,data_type in enumerate(results.keys()):
            if i == 0:
                columns_reference = cols_dict[data_type]
            else:
                if cols_dict[data_type] != columns_reference:
                    raise ValueError(f'Columns are not the same! {cols_dict}')

        print('  Done.')
        print()

        return ndarray_dict, cols_dict, n_jets_dict

    #---------------------------------------------------------------
    # Convert unfolding_results to format of results dict.
    #
    # The purpose of this is so that we can use the same machinery to apply 
    # cuts and make plots before and after unfolding.
    #
    # unfolding_results stores ndarrays of shape (n_jets, n_observables):
    #   - unfolding_results['nature_det'] = ndarray 
    #   - unfolding_results['sim_det'] = ndarray
    #   - unfolding_results['sim_truth'] = ndarray
    # as well as 1D arrays:
    #   - unfolding_results['weights_det_iteration{i}'] = 1darray
    #   - unfolding_results['weights_truth_iteration{i}'] = 1darray
    #   - unfolding_results['pt_hat_scale_factors'] = 1darray 
    # and the list of observables:
    #   - unfolding_results['columns'] = list of obs_keys 
    #
    # The results will be stored in a dict of 1D numpy arrays:
    #   - self.results[data_type][obs_key] = np.array([...])
    #       where data_type is: 'data', 'mc_det_matched', 'mc_truth_matched'
    #       and obs_key is self.observable_info[observable]['obs_key'][i] 
    #       (plus 'pt_hat_scale_factors' and 'weights_det/truth_iteration{i}' for MC)
    #
    # We use this dictionary of 1D numpy arrays as opposed to e.g. a dataframe
    # because we the different data types (data, sim_det, sim_truth) can have
    # different size.
    #---------------------------------------------------------------
    def convert_unfolding_results_to_dict(self, unfolding_results):

        results = self.recursive_defaultdict()
        obs_keys = unfolding_results['columns']

        for data_type in unfolding_results.keys():

            if data_type == 'columns':
                continue

            # If 1D array, store it in sim det/truth as appropriate
            # (store separately under sim det/truth so that cuts can easily be made independently)
            if len(unfolding_results[data_type].shape) == 1:
                if data_type == 'pt_hat_scale_factors':
                    results['sim_det'][data_type] = unfolding_results[data_type]
                    results['sim_truth'][data_type] = unfolding_results[data_type]
                elif 'weights_det' in data_type:
                        results['sim_det'][data_type] = unfolding_results[data_type]
                elif 'weights_truth' in data_type:
                    results['sim_truth'][data_type] = unfolding_results[data_type]
            
            # If ndarray, break into 1d arrays
            elif len(unfolding_results[data_type].shape) == 2:
                for i,obs_key in enumerate(obs_keys):
                    results[data_type][obs_key] = unfolding_results[data_type][:,i]

        n_jets = {}
        for data_type,val in results.items():
            if isinstance(val, np.ndarray):
                n_jets[data_type] = val.shape[0]
            else:
                n_jets[data_type] = val[list(val.keys())[0]].shape[0]

        return results, n_jets

    # ---------------------------------------------------------------
    # Get bin array (for each pt bin) from hepdata file
    # ---------------------------------------------------------------
    def bins_from_hepdata(self, block, jetR):

        if 'hepdata' in block:

            # Open the HEPData file
            hepdata_filename = block['hepdata']['file']        
            f = ROOT.TFile(hepdata_filename, 'READ')

            # The list of tables contains one entry per pt bin
            tables = block['hepdata']['tables'][jetR]
            h_name = block['hepdata']['hname']
            bins_list = []
            for table in tables:

                # Get the histogram, and return the bins
                if table:
                    dir = f.Get(table)
                    h = dir.Get(h_name)
                    bins = np.array(h.GetXaxis().GetXbins())

                    # For Soft Drop observables, we need to exclude the "untagged" bin so that it will become underflow
                    if 'SoftDrop' in block:
                        bins = bins[1:]

                else:
                    bins = np.array([])

                bins_list.append(bins)

            f.Close()

        else: 
            bins_list = None

        return bins_list

    # ---------------------------------------------------------------
    # Get tgraph (for each pt bin) from hepdata file
    # ---------------------------------------------------------------
    def tgraph_from_hepdata(self, block, jetR):

        if 'hepdata' in block:

            # Open the HEPData file
            hepdata_filename = block['hepdata']['file']        
            f = ROOT.TFile(hepdata_filename, 'READ')

            # The list of tables contains one entry per pt bin
            tables = block['hepdata']['tables'][jetR]
            g_name = block['hepdata']['gname']
            tgraph_list = []
            for table in tables:

                # Get the TGraph
                if table:
                    dir = f.Get(table)
                    g = dir.Get(g_name)
                else:
                    g = None

                tgraph_list.append(g)

            f.Close()

        else: 
            tgraph_list = None

        return tgraph_list

    #---------------------------------------------------------------
    # Divide a histogram by a tgraph, point-by-point
    #---------------------------------------------------------------
    def divide_histogram_by_tgraph(self, h, g, include_tgraph_uncertainties=True):

        # Truncate tgraph to range of histogram bins
        g_truncated = self.truncate_tgraph(g, h)
        if not g_truncated:
            return None

        # Clone tgraph, in order to return a new one
        g_new = g_truncated.Clone(f'{g_truncated.GetName()}_divided')

        nBins = h.GetNbinsX()
        for bin in range(1, nBins+1):

            # Get histogram (x,y)
            h_x = h.GetBinCenter(bin)
            h_y = h.GetBinContent(bin)
            h_error = h.GetBinError(bin)

            # Get TGraph (x,y) and errors
            gx, gy, yErrLow, yErrUp = self.get_gx_gy(g_truncated, bin-1)

            #print(f'h_x: {h_x}')
            #print(f'gx: {gx}')
            #print(f'h_y: {h_y}')
            #print(f'gy: {gy}')

            if not np.isclose(h_x, gx):
                print(f'WARNING: hist x: {h_x}, graph x: {gx} -- will not plot ratio')
                return None

            new_content = h_y / gy

            # Combine tgraph and histogram relative uncertainties in quadrature
            if gy > 0. and h_y > 0.:
                if include_tgraph_uncertainties:
                    new_error_low = np.sqrt( pow(yErrLow/gy,2) + pow(h_error/h_y,2) ) * new_content
                    new_error_up = np.sqrt( pow(yErrUp/gy,2) + pow(h_error/h_y,2) ) * new_content
                else:
                    new_error_low = h_error/h_y * new_content
                    new_error_up = h_error/h_y * new_content
            else:
                new_error_low = (yErrLow/gy) * new_content
                new_error_up = (yErrUp/gy) * new_content

            g_new.SetPoint(bin-1, h_x, new_content)
            g_new.SetPointError(bin-1, 0, 0, new_error_low, new_error_up)

        return g_new

    #---------------------------------------------------------------
    # Divide a tgraph by a tgraph, point-by-point: g1/g2
    # NOTE: Ignore uncertainties on denominator
    #---------------------------------------------------------------
    def divide_tgraph_by_tgraph(self, g1, g2):

        # Clone tgraph, in order to return a new one
        g_new = g1.Clone(f'{g1.GetName()}_divided')

        if g1.GetN() != g2.GetN():
            sys.exit(f'ERROR: TGraph {g1.GetName()} has {g1.GetN()} points, but {g2.GetName()} has {g2.GetN()} points')

        for i in range(0, g1.GetN()):

            # Get TGraph (x,y) and errors
            g1_x = ctypes.c_double(0)
            g1_y = ctypes.c_double(0)
            g1.GetPoint(i, g1_x, g1_y)
            y1ErrLow = g1.GetErrorYlow(i)
            y1ErrUp  = g1.GetErrorYhigh(i)
            g1x = g1_x.value
            g1y = g1_y.value

            g2_x = ctypes.c_double(0)
            g2_y = ctypes.c_double(0)
            g2.GetPoint(i, g2_x, g2_y)
            g2x = g2_x.value
            g2y = g2_y.value

            if not np.isclose(g1x, g2x):
                sys.exit(f'ERROR: TGraph {g1.GetName()} point {i} at {g1x}, but {g2.GetName()} at {g2x}')

            new_content = g1y / g2y
            new_error_low = y1ErrLow/g1y * new_content
            new_error_up = y1ErrUp/g1y * new_content

            g_new.SetPoint(i, g1x, new_content)
            g_new.SetPointError(i, 0, 0, new_error_low, new_error_up)
        return g_new

    #---------------------------------------------------------------
    # Truncate data tgraph to histogram binning range
    #---------------------------------------------------------------
    def truncate_tgraph(self, g, h):

        #print('truncate_tgraph')
        #print(h.GetName())
        #print(np.array(h.GetXaxis().GetXbins()))

        # Create new TGraph with number of points equal to number of histogram bins
        nBins = h.GetNbinsX()
        g_new = ROOT.TGraphAsymmErrors(nBins)
        g_new.SetName(f'{g.GetName()}_truncated')

        h_offset = 0
        for bin in range(1, nBins+1):

            # Get histogram (x)
            h_x = h.GetBinCenter(bin)

            # Get TGraph (x,y) and errors
            gx, gy, yErrLow, yErrUp = self.get_gx_gy(g, bin-1)

            #print(f'h_x: {h_x}')
            #print(f'gx: {gx}')
            #print(f'gy: {gy}')

            # If traph is offset from center of the bin, center it
            xErrLow = g.GetErrorXlow(bin-1)
            xErrUp = g.GetErrorXhigh(bin-1)
            if xErrLow > 0 and xErrUp > 0:
                x_min = gx - xErrLow
                x_max = gx + xErrUp
                x_center = (x_min + x_max)/2.
                if h_x > x_min and h_x < x_max:
                    if not np.isclose(gx, x_center):
                        gx = x_center

            # If tgraph starts below hist (e.g. when hist has min cut), try to get next tgraph point
            g_offset = 0
            while gx+1e-8 < h_x and g_offset < g.GetN()+1:
                g_offset += 1
                gx, gy, yErrLow, yErrUp = self.get_gx_gy(g, bin-1+g_offset)
            #print(f'new gx: {gx}')

            # If tgraph started above hist (see below) and we exhausted the tgraph points, skip
            if h_offset > 0 and np.isclose(gx, 0):
                continue

            # If tgraph starts above hist, try to get next hist bin
            h_offset = 0
            while gx-1e-8 > h_x and h_offset < nBins+1:
                h_offset += 1
                h_x = h.GetBinCenter(bin+h_offset)
                #print(f'h_x: {h_x}')
                #print(f'gx: {gx}')

            if not np.isclose(h_x, gx):
                print(f'WARNING: hist x: {h_x}, graph x: {gx}')
                return None

            g_new.SetPoint(bin-1, gx, gy)
            g_new.SetPointError(bin-1, 0, 0, yErrLow, yErrUp)
            #print()

        return g_new

    #---------------------------------------------------------------
    # Get points from tgraph by index
    #---------------------------------------------------------------
    def get_gx_gy(self, g, index):

        g_x = ctypes.c_double(0)
        g_y = ctypes.c_double(0)
        g.GetPoint(index, g_x, g_y)
        yErrLow = g.GetErrorYlow(index)
        yErrUp  = g.GetErrorYhigh(index)

        gx = g_x.value
        gy = g_y.value

        return gx, gy, yErrLow, yErrUp

    #---------------------------------------------------------------
    # Create a single output subdirectory
    #---------------------------------------------------------------
    def create_output_subdir(self, output_dir, name):

        output_subdir = os.path.join(output_dir, name)
        if not os.path.isdir(output_subdir):
            os.makedirs(output_subdir)

        return output_subdir