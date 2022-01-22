#!/usr/bin/env python

import ROOT

for r in [0.01, 0.1, 0.2, 0.3]:
  log1or = ROOT.TMath.Log(1./r)
  print('#line {x}, 0, {x}, 1, 1, 9, 1'.format(x=log1or))
