'''
Print data points from ROOT file
'''

import ROOT

filename = '/home/james/fFinalResults.root'
f = ROOT.TFile(filename, 'read')

filename_systematics = '/home/james/fSystematics.root'
f_sys = ROOT.TFile(filename_systematics, 'read')

#print(f.ls())
print(f_sys.ls())

beta=2
hname = 'hMain_zg_R0.4_zcut01_B{}_60.0-80.0'.format(beta)
hname_sys = 'hSystematic_zg_RegParam1_R0.4_zcut01_B{}_60.0-80.0_avg_new'.format(beta)

h = f.Get(hname)
h_sys = f_sys.Get(hname_sys)

for i in range(1, h.GetNbinsX()+1):
    low = h.GetXaxis().GetBinLowEdge(i)
    high = h.GetXaxis().GetBinUpEdge(i)
    content = h.GetBinContent(i)
    error_stat = h.GetBinError(i)
    error_sys = h_sys.GetBinContent(i)*content/100.
    print('{}-{}     {} +/- {} (stat)  +/- {} (sys)'.format(low, high, content, error_stat, error_sys))
