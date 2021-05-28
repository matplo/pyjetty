I am attaching two ROOT files for the TG3 plots, one for the girth and one for the mass, that contain the data and theory curves to be added to the JETSCAPE figures. Following the sample figures I sent in the previous email, I list how everything is named in the files below. 

Let me know if any modifications or additional details are needed and if I can be any help making the figures. It might be that we need to reduce the amount of theory comparisons to make the figures readable but I included what I thought would be nice for the discussion in the text and can modify as needed.

Thanks in advance!
Laura 

Girth: 
Top panel:
Pb-Pb data points+stat errors: (TH1D*)h
Pb-Pb data points+sys errors: (TGraphAsymmErrors*)g
PYTHIA Perugia 2011: (TGraph*)g_pyt11
Bottom panel (ratio):
Pb-Pb/pythia ratio points+stat err: (TH1D*)h_ratio_pyt11
Pb-Pb/pythia ratio points+sys err: (TGraphAsymmErrors*)g_ratio_pyt11
Hybrid model, Lres=0, w/o wake:  (TGraphAsymmErrors*)g_ratio_hyb_wake0_lres0
Hybrid model, Lres=0, w/ wake:  (TGraphAsymmErrors*)g_ratio_hyb_wake1_lres0
Hybrid model, Lres=2, w/o wake:  (TGraphAsymmErrors*)g_ratio_hyb_wake0_lres2
Hybrid model, Lres=inf, w/o wake:  (TGraphAsymmErrors*)g_ratio_hyb_wake0_lresinf

Mass:
Pb-Pb data points+stat errors: (TH1D*)h
Pb-Pb data points+sys errors: (TGraphAsymmErrors*)g
PYTHIA Perugia 2011: (TGraph*)g_pyt11
Hybrid model, Lres=0, w/o wake:  (TGraphAsymmErrors*)g_hyb_wake0_lres0
Hybrid model, Lres=0, w/ wake:  (TGraphAsymmErrors*)g_hyb_wake1_lres0
Hybrid model, Lres=inf, w/o wake:  (TGraphAsymmErrors*)g_hyb_wake0_lresinf
Hybrid model, Lres=inf, w/ wake:  (TGraphAsymmErrors*)g_hyb_wake1_lresinf
JEWEL, with recoils:  (TGraphAsymmErrors*)g_jewel_rec
JEWEL, no recoils:  (TGraphAsymmErrors*)g_jewel_norec
