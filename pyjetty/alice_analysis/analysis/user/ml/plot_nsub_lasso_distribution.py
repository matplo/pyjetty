import os
import pickle 
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd

#---------------------------------------------------------------------
# Function to translate a given N,beta to its index in the 2D nsubjettiness array
def index(N,beta):
    index = 3*(N-1)
    if not np.isclose(beta,0.5):
        index += beta
    return index

#---------------------------------------------------------------------
# Function to form the product observable from the 2D nsubjettiness array
def product_observable(X_train, N_list, beta_list, coeff_list, eps=0.):

    index_list = [index(N_list[i], beta_list[i]) for i in range(len(N_list))]
    print(f'index_list: {index_list}')

    observable = np.power(X_train[:,index_list[0]] + eps, coeff_list[0])
    for i in range(1, len(index_list)):
        observable = np.multiply(observable, np.power(X_train[:,index_list[i]] + eps, coeff_list[i]))

    return observable

#---------------------------------------------------------------------
# Function to set the binning range for each product observable
def bins(X, N_list):

    if len(N_list) == 1:
        bins = np.linspace(0., 0.03, 100)
    elif len(N_list) == 4:
        bins = np.linspace(0., 1.e-11, 100)
    elif len(N_list) == 23:
        bins = np.linspace(0., 1.e-28, 100)
    else:
        bins = np.linspace(np.amin(X), np.amax(X), 100)

    return bins

#---------------------------------------------------------------------
# Function to plot a given product observable
def plot_product_observable(N_list, beta_list, coeff_list, eps=0.):
    print()
    print(f'Plotting product observable with {len(N_list)} terms...')

    # Load observable and labels
    outputdir = '/rstorage/ml/egml/nsubjettiness/822239/hard_R0.4_pt[100.0, 125.0]_Rmax0'
    filename = 'tau_N.pkl'
    output_filename = os.path.join(outputdir, filename)
    with open(output_filename, 'rb') as f:
        X_4terms = pickle.load(f)   #observable = r'$\tau_{10}^{1})^{0.152} (\tau_{11}^{1})^{0.335} (\tau_{14}^{1})^{1.382} (\tau_{14}^{2})^{2.13}$'
                                    #index_list = [28, 31, 40, 41]
        X_train = pickle.load(f)
        y_train = pickle.load(f)
    print(f'X_nsub size: {X_train.shape}')

    # Define the product observable that we want to compute
    X = product_observable(X_train, N_list, beta_list, coeff_list, eps=eps)
    print(X)

    nonzero = np.count_nonzero(X)
    print(f'nonzero elements: {nonzero} (method1)')
    nonzero = (X > 1.e-50).sum()
    print(f'nonzero elements: {nonzero} (method2)')

    # Separate the PYTHIA/JEWEL samples and put them into a dataframe for plotting
    xlabel = rf'$\mathcal{{O}}^{{\mathrm{{ML}}}}_{{N-\mathrm{{sub}}}}$ ({len(N_list)} terms)' 
    ylabel = rf'$\frac{{1}}{{\sigma}} \frac{{d\sigma}}{{ d \mathcal{{O}}^{{\mathrm{{ML}}}}_{{N-\mathrm{{sub}}}} }}$'
    jewel_indices = y_train
    pythia_indices = 1 - y_train
    observable_jewel = X[jewel_indices.astype(bool)]
    observable_pythia = X[pythia_indices.astype(bool)]
    df_jewel = pd.DataFrame(observable_jewel, columns=[xlabel])
    df_pythia = pd.DataFrame(observable_pythia, columns=[xlabel])
    df_jewel['generator'] = np.repeat('JEWEL', observable_jewel.shape[0])
    df_pythia['generator'] = np.repeat('PYTHIA8', observable_pythia.shape[0])
    df = df_jewel.append(df_pythia, ignore_index=True)
    print(df.describe())

    # Plot and save
    h = sns.histplot(df, x=xlabel, hue='generator', stat='count', bins=bins(X,N_list), element='step', common_norm=False, log_scale=[False, True])
    if h.legend_:
        h.legend_.set_title(None)
        plt.setp(h.get_legend().get_texts(), fontsize='14') # for legend text
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(outputdir, f'tau_{len(N_list)}_terms.pdf'))
    plt.close()

#---------------------------------------------------------------------
# Main
#---------------------------------------------------------------------

# 1 term
N_list = [14]
beta_list = [1]
coeff_list = [1]
plot_product_observable(N_list, beta_list, coeff_list)

# 4 terms
N_list = [10, 11, 14, 14]
beta_list = [1, 1, 1, 2]
coeff_list = [0.152, 0.335, 1.382, 2.13]
plot_product_observable(N_list, beta_list, coeff_list)

# 23 terms
N_list = [2, 4, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14  ]
beta_list = [0.5, 2, 0.5, 2, 0.5, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 0.5, 2]
coeff_list = [1.053, -0.323, -0.145, -0.15, -0.312, 0.403, -0.746, 0.738, -1.275, 0.85, -1.341, 0.873, -1.219, 0.906, -1.379, 0.92, -1.273, 1.09, -1.785, 0.489, -0.521, 3.476, 1.733]
plot_product_observable(N_list, beta_list, coeff_list, eps=1.e-10)  