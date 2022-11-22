# Quick script to print p_value of ratio to test whether it is consistent with a constant fit (i.e. no modification)
# Here we ignore uncertainty correlations

import ROOT

result_labels = ['zg_0-10', 'zg_10-30', 'theta_g_0-10', 'theta_g_10-30']
tables = [3,6,9,12]

filename_base = 'HEPData-1627454721-v1-Table_{}.root'
gname = 'Graph1D_y1'

for i,result in enumerate(result_labels):

    result_label = result_labels[i]
    table = tables[i]
    filename = filename_base.format(table)

    f = ROOT.TFile(filename, 'r')
    dir = f.Get(f'Table {table}')
    g = dir.Get(gname)

    g.Fit("pol0", "q", "", 0., 1.)
    fit = g.GetFunction("pol0")
    fit.FixParameter(0, 1)

    #print(f'{result}: p-val={fit.GetProb():.3f}, chi2/dof={fit.GetChisquare()/fit.GetNDF():.3f}')
    print(f'{result}: p-val={fit.GetProb():.3f}')
    f.Close()