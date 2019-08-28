#!/usr/bin/env python
# coding: utf-8

# In[1]:


import fastjet as fj
import pythia8
from recursivetools import pyrecursivetools as rt
from pythiafjtools import pypythiafjtools as pyfj


# In[2]:


get_ipython().run_line_magic('matplotlib', 'widget')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
from tqdm import tnrange, tqdm_notebook


# In[3]:


def create_and_init_pythia(config_strings=[]):
    pythia = pythia8.Pythia()
    for s in config_strings:
        pythia.readString(s)
    for extra_s in ["Next:numberShowEvent = 0", "Next:numberShowInfo = 0", "Next:numberShowProcess = 0"]:
        pythia.readString(extra_s)
    if pythia.init():
        return pythia
    return None


# In[4]:


sconfig_pythia = [ "Beams:eCM = 8000.", "HardQCD:all = on", "PhaseSpace:pTHatMin = 20."]
pythia = create_and_init_pythia(sconfig_pythia)


# In[5]:


# set up our jet definition and a jet selector
jet_R0 = 0.4
jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
selector = fj.SelectorPtMin(20.0) & fj.SelectorPtMax(40.0) & fj.SelectorAbsEtaMax(1)
sd = rt.SoftDrop(0, 0.1, 1.0)


# In[6]:


# set up our jet definition and a jet selector
jet_R0 = 0.4
jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
jet_selector = fj.SelectorPtMin(20.0) & fj.SelectorPtMax(40.0) & fj.SelectorAbsEtaMax(1)
sd = rt.SoftDrop(0, 0.1, 1.0)


# In[7]:


all_jets = []
for iEvent in tqdm_notebook(range(100), 'event'):
    if not pythia.next(): continue
    parts = pyfj.vectorize(pythia, True, -1, 1, False)
    jets = jet_selector(jet_def(parts))
    all_jets.extend(jets)


# In[8]:


def deltas(jets, jets0):
    for i in range(len(jets)):
        yield jets0[i].perp() - jets[i].perp()


# In[9]:


get_ipython().run_cell_magic('time', '', 'all_sd_jets = [sd.result(j) for j in all_jets]\n\netas = [j.eta() for j in all_jets]\npts = [j.pt() for j in all_jets]\nsd_pts = [j.pt() for j in all_sd_jets]\nsd_delta_pt = [delta for delta in deltas(all_jets, all_sd_jets)]\n\nangs0 = [pyfj.angularity(j, 0.) for j in all_jets]\nsd_angs0 = [pyfj.angularity(j, 0.) for j in all_sd_jets]\nangs0_R0 = [pyfj.angularity(j, 0., jet_R0) for j in all_jets]\nsd_angs0_R0 = [pyfj.angularity(j, 0., jet_R0) for j in all_sd_jets]\n\nangs1 = [pyfj.angularity(j, 1.) for j in all_jets]\nsd_angs1 = [pyfj.angularity(j, 1.) for j in all_sd_jets]\nangs1_R0 = [pyfj.angularity(j, 1., jet_R0) for j in all_jets]\nsd_angs1_R0 = [pyfj.angularity(j, 1., jet_R0) for j in all_sd_jets]\n\nangs15 = [pyfj.angularity(j, 1.5) for j in all_jets]\nsd_angs15 = [pyfj.angularity(j, 1.5) for j in all_sd_jets]\nangs15_R0 = [pyfj.angularity(j, 1.5, jet_R0) for j in all_jets]\nsd_angs15_R0 = [pyfj.angularity(j, 1.5, jet_R0) for j in all_sd_jets]')


# In[10]:


fig, axes = plt.subplots(nrows=1, ncols=2, sharex=False, sharey=False)
ax0, ax1, = axes.flatten()
n, bins, patches = ax0.hist(pts, 25, density=1, facecolor='blue', alpha=0.3, label='anti-$k_{T}$ R=0.4')
n, bins, patches = ax0.hist(sd_pts, 25, density=1, facecolor='red', alpha=0.3, label='Soft Dropped (SD)')
# n, bins, patches = ax0.hist(sd_pts, 25, density=1, facecolor='red', alpha=0.3)
ax0.set_xlabel(r'$p_{T}$ (GeV)')
ax0.set_ylabel('Probability within $\hat{p_{T}} > 20$')
ax0.set_title(r'$\mathrm{PYTHIA\ jets}\ \sqrt{s}=8\ \mathrm{TeV}$ proton-proton')
ax0.grid(True)
ax0.legend(prop={'size': 10})
ax0.set_yscale('log')

n, bins, patches = ax1.hist(sd_delta_pt, 25, density=1, facecolor='green', alpha=0.3, label='$\Delta p_{T} = p_{T}^{SD} - p_{T}$')
ax1.legend(prop={'size': 10})
ax1.grid(True)
ax1.set_yscale('log')
fig.tight_layout()


# In[11]:


fig1, axes1 = plt.subplots(nrows=3, ncols=2, sharex=False, sharey=False)
ax10, ax11, ax12, ax13, ax14, ax15= axes1.flatten()

ax10.set_title(r'angularity $\alpha = 0$')
n, bins, patches = ax10.hist(angs0, 25, density=1, facecolor='blue', alpha=0.3)
n, bins, patches = ax10.hist(sd_angs0, 25, density=1, facecolor='red', alpha=0.3)

ax11.set_title(r'scaled by $R_{0}$')
n, bins, patches = ax11.hist(angs0_R0, 25, density=1, facecolor='blue', alpha=0.3)
n, bins, patches = ax11.hist(sd_angs0_R0, 25, density=1, facecolor='red', alpha=0.3)

ax12.set_title(r'angularity $\alpha = 1$')
n, bins, patches = ax12.hist(angs1, 25, density=1, facecolor='blue', alpha=0.3)
n, bins, patches = ax12.hist(sd_angs1, 25, density=1, facecolor='red', alpha=0.3)

ax13.set_title(r'scaled by $R_{0}$')
n, bins, patches = ax13.hist(angs1_R0, 25, density=1, facecolor='blue', alpha=0.3)
n, bins, patches = ax13.hist(sd_angs1_R0, 25, density=1, facecolor='red', alpha=0.3)

ax14.set_title(r'angularity $\alpha = 1.5$')
n, bins, patches = ax14.hist(angs15, 25, density=1, facecolor='blue', alpha=0.3)
n, bins, patches = ax14.hist(sd_angs15, 25, density=1, facecolor='red', alpha=0.3)

ax15.set_title(r'scaled by $R_{0}$')
n, bins, patches = ax15.hist(angs15_R0, 25, density=1, facecolor='blue', alpha=0.3)
n, bins, patches = ax15.hist(sd_angs15_R0, 25, density=1, facecolor='red', alpha=0.3)

fig1.tight_layout()


# In[ ]:





# In[ ]:




