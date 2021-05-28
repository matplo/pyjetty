In terms of replotting, we're looking for three plots for our section: hadron Raa w/ models, jet Raa w/ models, and then a single plot with both hadron and jet Raa, all at 5 TeV. Most of the data is in an attached root file, except for the JEWEL points. The attached `info.txt` describes where all of the measured data are from (hepdata, dois, etc), most of which I've transcribed here. In any case, it's included for completeness. This is all thanks to Marta, who compiled and provided all of this. Everything is explained below.

The observables included in the root file are:
1. ALICE 0-5% charged hadron Raa: grALICEHad, grALICEHadSys: TGraphAsymmErrors. Source: doi:10.17182/hepdata.86210.v1/t8
2. CMS 0-5% charged hadron Raa: grCMSHad, grCMSHadSys: TGraphAsymmErrors. Source:  doi:10.17182/hepdata.77101.v1/t8
3. ALICE 0-10% R = 0.4 full jet Raa (which you of course already have, but are in the file in any case): grALICEJet, grALICEJetCorr, grALICEJetShape: TGraphErrors
4. ATLAS 0-10% R = 0.4 full jet Raa: grATLASJet, grATLASJetSys: TGraphAsymmErrors. Source: doi:10.17182/hepdata.84819.v1/t19

Then, the JEWEL hadron Raa points are included in RAA_jewel_lhc.txt . It just has bin centers and values. I think it's with recoils, but to be confirmed. However, in compiling today, we found out this may be at 2.76 TeV. So for now, let's use this, but it may be updated.

In terms of plots, we would like:

Hadron + Jet Raa (currently, Fig 7. Note that if the styles for data aren't consistent, then it may not be necessary to remake this figure):
#1, #2, #3, and #4

Hadron Raa:
#1, jetscape hadron Raa, JEWEL hadron Raa

Jet Raa:
ALICE R = 0.2 full jet Raa (ie. not #3), JEWEL w/ and w/out recoils, Hybrid model Lres = 0, 2/(pi T) (these aren't included in this message), jetscape jet Raa

This may be tweaked a little bit in the coming days, but this should be it, more or less.
