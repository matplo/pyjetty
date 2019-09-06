#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import argparse
import os
import numpy as np
import pandas
import array 

import pyhepmc_ng
import ROOT
import pu14_reader

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

#--------------------------------------------------------------
def main():
  parser = argparse.ArgumentParser(description='jetscape in python', \
                                   prog=os.path.basename(__file__))
  parser.add_argument('-i', '--input', help='input file', \
                      default='low', type=str, required=True)
  parser.add_argument('--nev', help='number of events', \
                      default=1000, type=int)
  args = parser.parse_args()	

  plot_hadrons = False
  plot_final_state_partons = True
  plot_nth_partons = False
  include_thermal_background = False
  
  # Initialize histogram dictionary
  hDict = initializeHistograms()
  
  # Use pyhepmc_ng to parse the HepMC file
  input_hepmc = pyhepmc_ng.ReaderAscii(args.input)
  if input_hepmc.failed():
    print ("[error] unable to read from {}".format(args.input))
    sys.exit(1)

  # jet finder
  fj.ClusterSequence.print_banner()
  print()
  jetR = 0.4
  jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
  jet_selector = fj.SelectorPtMin(50.0) & fj.SelectorAbsEtaMax(2)

  # Jet re-clustering definition, for primary Lund plane
  jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, jetR)
  lund_gen = fjcontrib.LundGenerator(jet_def_lund)

  # Load thermal background: array of dataframes, one per event
  if include_thermal_background:
    thermal_background_array = get_thermal_background()
  
  # Loop through events
  all_jets_hadron = []
  all_jets_parton = []
  kt_shower_list = []
  event_hepmc = pyhepmc_ng.GenEvent()
  pbar = tqdm.tqdm(range(args.nev))
  while not input_hepmc.failed():
    ev = input_hepmc.read_event(event_hepmc)
    if input_hepmc.failed():
      nstop = pbar.n
      pbar.close()
      print('End of HepMC file at event {} '.format(nstop))
      break

    if plot_hadrons:
      hadrons = get_hadrons(event_hepmc)
      jets_hadron = find_jets(jet_def, jet_selector, hadrons)
      all_jets_hadron.extend(jets_hadron)
    if plot_final_state_partons:
      parton_list = get_final_partons(event_hepmc, hDict)
    elif plot_nth_partons:
      # Get the first parton in the shower (it is apparently the child of the initiating parton)
      first_parton = []
      first_parton.append(get_first_parton(event_hepmc, hDict).children[0])
      #print('-----------------------------')
      #print('First parton E: {}'.format(first_parton[0].momentum.e))

      n = 3
      parton_list = get_nth_partons(first_parton, n)
      #print('Total number of partons: {}'.format(len(parton_list)))      
      
    fj_particles = []
    for parton in parton_list:
      fj_particles.append(fj.PseudoJet(parton.momentum.px, parton.momentum.py, parton.momentum.pz, parton.momentum.e))

    # Append thermal background
    if include_thermal_background:
      thermal_background_df = thermal_background_array[0]
      thermal_background_fjparticles = get_fjparticles(thermal_background_df)
      fj_particles.extend(thermal_background_fjparticles)

    jets_parton = find_jets(jet_def, jet_selector, fj_particles)
    all_jets_parton.extend(jets_parton)
    
    # Analyze shower history
    max_kt = get_max_kt_shower(event_hepmc, hDict)
    kt_shower_list.append(max_kt)

    pbar.update()

    if pbar.n >= args.nev:
      pbar.close()
      print('{} event limit reached'.format(args.nev))
      break

  print('Constructing histograms...')

  if plot_hadrons:
    n_jets_hadron = len(all_jets_hadron)
    print('n_jets_hadron: {}'.format(n_jets_hadron))
  n_jets_parton = len(all_jets_parton)
  print('n_jets_parton: {}'.format(n_jets_parton))
        
  # Fill histogram
  if plot_hadrons:
    [fill_jet_histogram(hDict, jet) for jet in all_jets_hadron]

    # Fill Lund diagram
    # Note: l in lunds, l is the list of splittings in a given jet (following hardest splitting)
    lunds_hadron = [lund_gen.result(jet) for jet in all_jets_hadron]
    [fill_lund_histogram(hDict, "hLundHadron", splitting_list) for splitting_list in lunds_hadron]
    hDict['hLundHadron'].Scale(1./n_jets_hadron, "width")
  
  lunds_parton = [lund_gen.result(jet) for jet in all_jets_parton]
  [fill_lund_histogram(hDict, "hLundParton", splitting_list) for splitting_list in lunds_parton]
  hDict['hLundParton'].Scale(1./n_jets_parton, "width")

  
  for splitting_list in lunds_parton:
    max_kt_recluster = 0
    for split in splitting_list:
      kt = split.kt()
      if kt > max_kt_recluster:
        max_kt_recluster = kt
    hDict['hMaxKT_2D'].Fill(np.log(max_kt_recluster), np.log(kt))

  
  plot_histograms(hDict)

#---------------------------------------------------------------
def get_max_kt_shower(event_hepmc, hDict):

  # Loop over all vertices in the HepMC graph, and compute kt for those which are 1-->2 splittings
  max_kt = 0
  for vertex in event_hepmc.vertices:
    parent = vertex.particles_in
    children = vertex.particles_out
    
    # Check if we have a splitting 
    if len(parent) == 1 and len(children) == 2:

      soft_child = None
      if children[0].momentum.pt() < children[1].momentum.pt():
        soft_child = children[0]
      else:
        soft_child = children[1]

      pt_parent = parent[0].momentum.pt()
      pt_child = soft_child.momentum.pt()
      theta_parent = parent[0].momentum.theta()
      theta_child = soft_child.momentum.theta()
      phi_parent = parent[0].momentum.phi()
      phi_child = soft_child.momentum.phi()
        
      z = pt_child / pt_parent
      delta_R = np.sqrt( (theta_parent - theta_child)**2 + (phi_parent - phi_child)**2 )
      kt = delta_R * pt_child
      
      if kt > max_kt:
        max_kt = kt

  hDict['hMaxKT'].Fill(np.log(max_kt))
  return max_kt
      
#---------------------------------------------------------------
# Use swig'd function to create a vector of fastjet::PseudoJets from numpy arrays of px,py,pz,e
def get_thermal_background():

  # Get thermal particle four-vectors from PU14 file
  file = '/lustre/emmi/emmi05/thermal/Mult1700/ThermalEventsMult1700PtAv0.90_0.pu14'
  n_particles = 1700
  print('\n Loading PU14 background with multiplicity {} from: {}'.format(n_particles, file))
  df_array = pu14_reader.read_file(file, n_particles, n_events_max=10)
  print('Done.')
  return df_array

#---------------------------------------------------------------
# Get fastjet particles from dataframe of four-vectors
def get_fjparticles(df):

  px_array = df['px'].values
  py_array = df['py'].values
  pz_array = df['pz'].values
  m_array = df['m'].values
  fj_particles = fjext.vectorize_px_py_pz_m(px_array, py_array, pz_array, m_array)

  return fj_particles
  
#---------------------------------------------------------------
# Get a list of the children of a given parton.
# If the parton has no children, return the parton itself.
def get_children(parton):

  # First, check if the parton has no children -- if so, return the parton
  # Then, check if the parent has the same energy as the child -- if so skip down to the child
  children_list = []
  children = parton.children
  if len(children) is 0:
    children_list.append(parton)
  elif len(children) is 1 and abs(children[0].momentum.e - parton.momentum.e) < 1e-5:
    #print('parent E == child E ... move one level down')
    return get_children(children[0])
  else:
    children_list.extend(children)

  return children_list
    
#---------------------------------------------------------------
def get_nth_partons(parton_list, n):

  #print('get_nth_partons from {} partons at level {}'.format(len(parton_list), n))
  
  if n is 0:
    return parton_list
  else:
    children_list = []
    for parton in parton_list:
      children = get_children(parton)
      children_list.extend(children)
      
      #print('parent E: {}'.format(parton.momentum.e))
      #for child in children:
        #print('child E: {}'.format(child.momentum.e))

    return get_nth_partons(children_list, n-1)
  
#---------------------------------------------------------------
def initializeHistograms():

  hDict = {}

  lbins = logbins(1., 500, 50)
  hJetPt = ROOT.TH1F('hJetPt', 'hJetPt', 50, lbins)
  hJetPt.GetXaxis().SetTitle('p_{T,jet}')
  hJetPt.GetYaxis().SetTitle('dN/dp_{T}')
  hDict['hJetPt'] = hJetPt

  hStatus = ROOT.TH2D('hStatus', 'hStatus', 100, 0, 100, 100, 0, 100)
  hStatus.GetXaxis().SetTitle('HepMC particle status')
  hStatus.GetYaxis().SetTitle('n_children')
  hDict['hStatus'] = hStatus

  hStatusParton = ROOT.TH2D('hStatusParton', 'hStatusParton', 100, 0, 100, 2, 0, 2)
  hStatusParton.GetXaxis().SetTitle('HepMC particle status')
  hStatusParton.GetYaxis().SetTitle('is_parton')
  hDict['hStatusParton'] = hStatusParton
        
  hVertices = ROOT.TH2D('hVertices', 'hVertices', 1000, 0, 1000, 1000, 0, 1000)
  hVertices.GetXaxis().SetTitle('n_vertices')
  hVertices.GetXaxis().SetTitle('n_particles')
  hDict['hVertices'] = hVertices

  hLundHadron = ROOT.TH2D("hLundHadron", "hLundHadron", 80, 0, 8, 60, 0, 6)
  hLundHadron.GetXaxis().SetTitle('ln(1/#Delta R)')
  hLundHadron.GetYaxis().SetTitle('ln(k_{T})')
  hDict['hLundHadron'] = hLundHadron

  hLundParton = ROOT.TH2D("hLundParton", "hLundParton", 80, 0, 8, 60, 0, 6)
  hLundParton.GetXaxis().SetTitle('ln(1/#Delta R)')
  hLundParton.GetYaxis().SetTitle('ln(k_{T})')
  hDict['hLundParton'] = hLundParton

  #hMaxKT = ROOT.TH2D("hMaxKT", "hMaxKT", 60, 0, 6, 60, 0, 6)
  hMaxKT = ROOT.TH1D("hMaxKT", "hMaxKT", 60, 0, 6)
  hMaxKT.GetXaxis().SetTitle('ln[max(k_{T})] shower')
  #hMaxKT.GetYaxis().SetTitle('ln[max(k_{T})] re-cluster')
  hDict['hMaxKT'] = hMaxKT
  
  for key, val in hDict.items():
    val.Sumw2()
  
  return hDict
  
#--------------------------------------------------------------
def get_hadrons(hepmc_event):

  fjparts = []
  hadrons = []
  for vertex in hepmc_event.vertices:
    vertex_time = vertex.position.t
    if abs(vertex_time - 100) < 1e-3:
      hadrons = vertex.particles_out

  for hadron in hadrons:
    psj = fj.PseudoJet(hadron.momentum.px, hadron.momentum.py, hadron.momentum.pz, hadron.momentum.e)
    fjparts.append(psj)

  return fjparts

#--------------------------------------------------------------
def get_final_partons(hepmc_event, hDict):

  partons = []
  n_vertices = len(hepmc_event.vertices)
  n_particles = len(hepmc_event.particles)
  hDict['hVertices'].Fill(n_vertices, n_particles)
  for particle in hepmc_event.particles:
    status = particle.status
    n_children = len(particle.children)
    hDict['hStatus'].Fill(status, n_children)

    parent_vertex = particle.production_vertex
    end_vertex = particle.end_vertex
    parent_vertex_time = parent_vertex.position.t 
    is_parton = False
    is_final_parton = False
    if abs(parent_vertex_time - 100) > 1e-3:
      is_parton = True
    if is_parton and not end_vertex:
      is_final_parton = True
      
    hDict['hStatusParton'].Fill(status, is_parton)
    # It seems that status=0 means the particle is a parton
    # Not sure what the nonzero status codes mean (62, 83, 84) 
    
    if is_final_parton:
      partons.append(particle)

  return partons

#--------------------------------------------------------------
def get_first_parton(hepmc_event, hDict):

  for particle in hepmc_event.particles:
    parent_vertex = particle.production_vertex
    parent_vertex_time = parent_vertex.position.t
    n_parents = len(particle.parents)
    is_parton = False
    is_first_parton = False
    if abs(parent_vertex_time - 100) > 1e-3:
      is_parton = True
    if is_parton and n_parents is 0:    
      is_first_parton = True

    if is_first_parton:
      return particle

  return None

#--------------------------------------------------------------
def find_jets(jet_def, jet_selector, fjparts):

  jets = jet_selector(jet_def(fjparts))
  return jets

#--------------------------------------------------------------
def fill_jet_histogram(hDict, jet):

  hDict['hJetPt'].Fill(jet.perp())

#--------------------------------------------------------------
def fill_lund_histogram(hDict, hist_name, splitting_list):

  for split in splitting_list:
    z = split.z()
    kt = split.kt()
    delta_R = split.Delta() 
    if delta_R < 1e-3:
      print('delta_R = 0!')
      continue
    hDict[hist_name].Fill(np.log(1/delta_R), np.log(kt))

#--------------------------------------------------------------
def plot_histograms(hDict):

  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)
  
  output_filename = "hJetPt.pdf"
  plotHist(hDict['hJetPt'], output_filename, 'width', setLogy=True)

  output_filename = "hLundHadron.pdf"
  plotHist(hDict['hLundHadron'], output_filename, "colz")

  output_filename = "hLundParton.pdf"
  plotHist(hDict['hLundParton'], output_filename, "colz")
        
  fout = ROOT.TFile("AnalysisResults.root", "recreate")
  fout.cd()
  for key, val in hDict.items():
    val.Write()
  fout.Close()

#---------------------------------------------------------------
def plotHist(h, outputFilename, drawOptions = "", setLogy = False, setLogz = False):

  h.SetLineColor(1)
  h.SetLineWidth(1)
  h.SetLineStyle(1)

  c = ROOT.TCanvas("c","c: hist",600,450)
  c.cd()
  ROOT.gPad.SetLeftMargin(0.15)
  if setLogy:
    c.SetLogy()
  if setLogz:
    c.SetLogz()
  ROOT.gPad.SetLeftMargin(0.15)
        
  h.Draw(drawOptions)
  c.SaveAs(outputFilename)
  c.Close()
  
#--------------------------------------------------------------
def logbins(xmin, xmax, nbins):
  lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
  arr = array.array('f', lspace)
  return arr
                
#--------------------------------------------------------------
if __name__ == '__main__':
  main()
