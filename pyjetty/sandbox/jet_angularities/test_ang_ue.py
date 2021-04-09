#!/usr/bin/env python

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import fastjet as fj
import fjext
import fjcontrib
import fjtools

import pythia8
import pythiafjext
import pythiaext
from heppy.pythiautils import configuration as pyconf

# from tqdm.notebook import tqdm
from tqdm import tqdm
import argparse
import os
import sys


def get_args_from_settings(ssettings):
    sys.argv=[' '] + ssettings.split()
    parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly')
    pyconf.add_standard_pythia_args(parser)
    parser.add_argument('--output', default="test_ang_ue.root", type=str)
    parser.add_argument('--user-seed', help='pythia seed', default=1111, type=int)
    args = parser.parse_args()
    return args

def matched(j1, j2):
    # j2i = fjtools.matched_Ry(j1, j2)
    mpt = fjtools.matched_pt(j1, j2)
    if mpt > 0.5:
        return True, j1, j2, fjext.lambda_beta_kappa(j1, 1.0, 1.0, 1.0), fjext.lambda_beta_kappa(j2, 1.0, 1.0, 1.0)
    return False

def fill_matched(j1s, j2s, tj_no_pup, tj_pup, tj_delta, jet_R0):
    for j1 in j1s:
        tj_no_pup.Fill(j1.perp(), j1.eta(), j1.phi(), 
                       fjext.lambda_beta_kappa(j1, 1.0, 1.0, jet_R0),
                       fjext.lambda_beta_kappa(j1, 2.0, 1.0, jet_R0),
                       fjext.lambda_beta_kappa(j1, 3.0, 1.0, jet_R0))
        for j2 in j2s:
            mpt = fjtools.matched_pt(j1, j2)
            tj_delta.Fill(j1.perp(), j1.eta(), j1.phi(), 
                          fjext.lambda_beta_kappa(j1, 1.0, 1.0, jet_R0),
                          fjext.lambda_beta_kappa(j1, 2.0, 1.0, jet_R0),
                          fjext.lambda_beta_kappa(j1, 3.0, 1.0, jet_R0),
                          j2.perp(), j2.eta(), j2.phi(),
                          fjext.lambda_beta_kappa(j2, 1.0, 1.0, jet_R0),
                          fjext.lambda_beta_kappa(j2, 2.0, 1.0, jet_R0),
                          fjext.lambda_beta_kappa(j2, 3.0, 1.0, jet_R0), 
                          mpt)
    for j1 in j2s:
        tj_pup.Fill(j1.perp(), j1.eta(), j1.phi(),
                    fjext.lambda_beta_kappa(j1, 1.0, 1.0, jet_R0),
                    fjext.lambda_beta_kappa(j1, 2.0, 1.0, jet_R0),
                    fjext.lambda_beta_kappa(j1, 3.0, 1.0, jet_R0))

def main():
    mycfg = []
    ssettings = "--py-ecm 5000 --py-minbias --user-seed=100000"
    args = get_args_from_settings(ssettings)
    pythia_mb = pyconf.create_and_init_pythia_from_args(args, mycfg)

    mycfg = []
    ssettings = "--py-ecm 5000 --user-seed=100000 --nev 100000"
    args = get_args_from_settings(ssettings)
    pythia_hard = pyconf.create_and_init_pythia_from_args(args, mycfg)


    max_eta_hadron=1
    parts_selector_h = fj.SelectorAbsEtaMax(max_eta_hadron)
    jet_R0 = 0.4
    jet_selector = fj.SelectorPtMin(20.0) & fj.SelectorPtMax(200.0) & fj.SelectorAbsEtaMax(max_eta_hadron - 1.05 * jet_R0)
    n_pileup = 1 #5

    # print the banner first
    fj.ClusterSequence.print_banner()
    print()
    # set up our jet definition and a jet selector
    jet_R0 = 0.4
    jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
    print(jet_def)

    fout = ROOT.TFile(args.output, 'recreate')
    fout.cd()
    tj_delta = ROOT.TNtuple("tj_delta", "tj_delta", "pt:eta:phi:L11:L21:L31:ptm:etam:phim:L11m:L21m:L31m:mpt")
    tj_no_pup = ROOT.TNtuple("tj_no_pup", "tj_no_pup", "pt:eta:phi:L11:L21:L31")
    tj_pup = ROOT.TNtuple("tj_pup", "tj_pup", "pt:eta:phi:L11:L21:L31")
    hmult_hard = ROOT.TH1F("hmult_hard", "hmult_hard", 300, 0, 300)
    hmult_pup = ROOT.TH1F("hmult_pup", "hmult_pup", 300, 0, 300)
    hpt_acc_hard = ROOT.TProfile2D("hpt_acc_hard", "hpt_acc_hard;#eta;#varphi", 50, -1, 1, 50, 0, ROOT.TMath.Pi() * 2.)
    hpt_acc_pup = ROOT.TProfile2D("hpt_acc_pup", "hpt_acc_pup;#eta;#varphi", 50, -1, 1, 50, 0, ROOT.TMath.Pi() * 2.)

    for n in tqdm(range(args.nev)):
        if not pythia_hard.next():
            continue
        parts_pythia_h = pythiafjext.vectorize_select(pythia_hard, [pythiafjext.kFinal, pythiafjext.kCharged], 0, False)
        parts_pythia_h_selected = parts_selector_h(parts_pythia_h)

        parts_pileup = None
        for ipile in range(n_pileup):
            while not pythia_mb.next():
                continue
            parts_pythia_h_ue = pythiafjext.vectorize_select(pythia_mb, [pythiafjext.kFinal, pythiafjext.kCharged], 10000, False)
            parts_pythia_h_selected_ue = parts_selector_h(parts_pythia_h_ue)
            if parts_pileup is None:
                parts_pileup = parts_pythia_h_selected_ue
            else:
                parts_pileup += parts_pythia_h_selected_ue

        mult_hard = len(parts_pythia_h_selected)
        mult_ue = len(parts_pileup)
        
        jets_h = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_h_selected)))
        jets_h_w_ue = fj.sorted_by_pt(jet_selector(jet_def(parts_pileup + parts_pythia_h_selected)))

        if len(jets_h) < 1:
            continue
            
        fill_matched(jets_h, jets_h_w_ue, tj_no_pup, tj_pup, tj_delta, jet_R0)

        hmult_hard.Fill(mult_hard)
        hmult_pup.Fill(mult_ue)

        _tmp = [hpt_acc_hard.Fill(p.eta(), p.phi(), p.perp()) for p in parts_pythia_h_selected]
        _tmp = [hpt_acc_pup.Fill(p.eta(), p.phi(), p.perp()) for p in parts_pileup]

    pythia_hard.stat()
    pythia_mb.stat()

    fout.Write()
    fout.Close()
    print ('[i] written ', fout.GetName())


if __name__ == '__main__':
    main()
