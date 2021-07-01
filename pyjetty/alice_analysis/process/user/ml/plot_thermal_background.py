#!/usr/bin/env python3

"""
Plot thermal background, to determine parameter values
"""

import os
import sys
import argparse

# Data analysis and plotting
import numpy as np
from matplotlib import pyplot as plt

# Analysis utilities
from pyjetty.alice_analysis.process.base import thermal_generator
from pyjetty.alice_analysis.process.base import process_base
from pyjetty.alice_analysis.process.base import process_utils

# Base class
from pyjetty.alice_analysis.process.base import common_base

################################################################
class PlotThermalBackground(common_base.CommonBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, config_file='', input_file='', output_dir='', **kwargs):
        super(common_base.CommonBase, self).__init__(**kwargs)
        
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            
        # Initialize utils class
        self.utils = process_utils.ProcessUtils()
        
        print(self)
        print()

    #---------------------------------------------------------------
    # Initialize thermal model parameters
    # ALICE measures dNch/deta~1700 for 0-10% at 5.02 TeV
    # So for full jets, could estimate dNch/deta~2500
    #---------------------------------------------------------------
    def plot_thermal_background(self):
    
        self.n_events = 1000
        self.eta_max = 0.5
        N_avg_list = [2000., 2500., 3000.]
        beta_list = [0.35, 0.4, 0.5]
        
        for N_avg in N_avg_list:
            for beta in beta_list:
                self.plot_delta_pt(N_avg, beta)
                
    #---------------------------------------------------------------
    # Generate thermal background and compute delta-pt by random cone method
    #---------------------------------------------------------------
    def plot_delta_pt(self, N_avg, beta):
   
        self.thermal_generator = thermal_generator.ThermalGenerator(N_avg=N_avg,
                                                                    sigma_N=500,
                                                                    beta=beta,
                                                                    eta_max=self.eta_max)
        
        # Loop through events
        self.delta_pt_random_cone = []
        self.mean_pt = []
        for i in range(self.n_events):
            fj_particles_background = self.thermal_generator.load_event()

            # Compute delta-pt by random cone method
            self.delta_pt_RC(fj_particles_background, N_avg)
            
        # Plot and save
        mean = np.round(np.mean(self.delta_pt_random_cone),2)
        sigma = np.round(np.std(self.delta_pt_random_cone),2)
        plt.hist(self.delta_pt_random_cone,
                 np.linspace(-50, 50, 100),
                 histtype='stepfilled',
                 label = rf'$\mathrm{{mean}} = {mean},\;\sigma = {sigma}$',
                 linewidth=2,
                 linestyle='-',
                 alpha=0.5)
        plt.xlabel(r'$\delta p_{T}$', fontsize=14)
        plt.yscale('log')
        legend = plt.legend(loc='best', fontsize=14, frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, f'delta_pt_random_cone_N{N_avg}_beta{beta}.pdf'))
        plt.close()
        
        print(f'mean pt (N={N_avg},beta={beta}): {np.round(np.mean(self.mean_pt),2)}')
           
    #---------------------------------------------------------------
    # Compute delta-pt by random cone method
    #---------------------------------------------------------------
    def delta_pt_RC(self, fj_particles_background, N_avg):
    
        event_pt = 0.
        cone_pt = 0.
        R_cone = 0.4
        eta_random = np.random.uniform(-self.eta_max+R_cone, self.eta_max-R_cone)
        phi_random = np.random.uniform(0, 2*np.pi)
        for particle in fj_particles_background:
            event_pt += particle.pt()
            delta_R = self.utils.delta_R(particle, eta_random, phi_random)
            if delta_R < R_cone:
                cone_pt += particle.pt()
        rho = event_pt / (2*self.eta_max *2*np.pi)
        delta_pt = cone_pt - rho*np.pi*R_cone*R_cone
        
        self.delta_pt_random_cone.append(delta_pt)
        self.mean_pt.append(event_pt/N_avg)

##################################################################
if __name__ == '__main__':

    # Define arguments
    parser = argparse.ArgumentParser(description='Process pp AA')
    parser.add_argument('-o', '--outputDir', action='store',
                        type=str, metavar='outputDir',
                        default='./TestOutput',
                        help='Output directory for output to be written to')

    # Parse the arguments
    args = parser.parse_args()

    print('Configuring...')
    print('ouputDir: \'{0}\"'.format(args.outputDir))

    analysis = PlotThermalBackground(output_dir=args.outputDir)
    analysis.plot_thermal_background()
