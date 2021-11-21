import os
import pickle 
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd

# Load observable and labels
outputdir = '/rstorage/ml/egml/nsubjettiness/822239/hard_R0.4_pt[100.0, 125.0]_Rmax0'
filename = 'tau_10_11_14_14.pkl'
output_filename = os.path.join(outputdir, filename)
with open(output_filename, 'rb') as f:
    X = pickle.load(f)
    y_train = pickle.load(f)
print(f'observable: {X}')
print(f'label: {y_train}')
print(f'size: {X.shape}')

nonzero = np.count_nonzero(X)
print(f'nonzero elements: {nonzero}')

nonzero = (X > 1.e-30).sum()
print(f'nonzero elements: {nonzero}')

# Separate the PYTHIA/JEWEL samples and put them into a dataframe for plotting
xlabel = rf'$\mathcal{{O}}^{{\mathrm{{ML}}}}_{{N-\mathrm{{sub}}}}$ (4 terms)' 
ylabel = rf'$\frac{{1}}{{\sigma}} \frac{{d\sigma}}{{ d \mathcal{{O}}^{{\mathrm{{ML}}}}_{{N-\mathrm{{sub}}}} }}$'
jewel_indices = y_train
pythia_indices = 1 - y_train
observable_jewel = X[jewel_indices.astype(bool)]
observable_pythia = X[pythia_indices.astype(bool)]
df_jewel = pd.DataFrame(observable_jewel, columns=[xlabel])
df_pythia = pd.DataFrame(observable_pythia, columns=[xlabel])
df_jewel['generator'] = np.repeat(1, observable_jewel.shape[0])
df_pythia['generator'] = np.repeat(0, observable_pythia.shape[0])
df = df_jewel.append(df_pythia, ignore_index=True)
print(df.describe())

# Plot and save
bins = np.linspace(0., 1.e-11, 100)
#bins = np.linspace(np.amin(X), np.amax(X), 50)
h = sns.histplot(df, x=xlabel, hue='generator', stat='count', bins=bins, element='step', common_norm=False, log_scale=[False, True])
if h.legend_:
    h.legend_.set_title(None)
    plt.setp(h.get_legend().get_texts(), fontsize='14') # for legend text
plt.xlabel(xlabel, fontsize=12)
plt.ylabel(ylabel, fontsize=16)
plt.tight_layout()
plt.savefig(os.path.join(outputdir, 'tau_10_11_14_14.pdf'))
plt.close()