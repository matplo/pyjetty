# This file takes the data generated from Yue Shi (many h5 files) and reformats it such that it can be used in the ML code
# The format at the end is similar to the quark/gluon test data set
# The data set at the end containes labeled pp and AA jets
# There will be an equal number of pp and AA jets and the data set is shuffled already
# The final format is (jets, particles, 4-vectors), where the 4-vector format is (px,py,pz,E)

import os
import h5py
import numpy as np

# ---------------------------------------------------------
# Yue Shi's files contain 5988 Pythia h5 files and 2176 Jewel h5 files
# Note: The number of jets in each file is different and there are generally more jets in Pythia files
# Choose how many h5 files to read in (here same for both Pythia and Jewel)

nh5 = 2176

# ---------------------------------------------------------
# First: read in Pythia data files and do some reformatting
# ---------------------------------------------------------

# Initialize data structure
data_pythia = np.zeros((1,800,4))
#print(data_pythia.shape)

print('Reading in Pythia h5 files')
for i,file in enumerate(os.listdir('/rstorage/ml/egml/pythia8/')):
    if file.endswith('.h5') and i<nh5:
        with h5py.File(os.path.join('/rstorage/ml/egml/pythia8/', file),'r') as hdf:
            # jetsubs_cons are the 4-vectors of particles inside the jets
            dataset_pythia0 = hdf.get('jetsubs_cons')
            
            # Change format: jet, particle, 4-vector
            dataset_pythia = np.transpose(np.array(dataset_pythia0),(0,2,1)) 
            
            # Reorder 4-vector from (px,py,pz,E) to (E,px,py,pz)
            dataset_pythia[:,:,[0,1,2,3]] = dataset_pythia[:,:,[3,0,1,2]]
            
            # Replace nan with zero
            where_are_nan = np.isnan(dataset_pythia)
            dataset_pythia[where_are_nan] = 0.
            
            # Zero pad such that they all have the same shape (.., 800, 4)
            dataset_pythia1 = np.pad(dataset_pythia,((0,0),(0,800-dataset_pythia.shape[1]),(0,0)),'constant',constant_values=(0))
            print(dataset_pythia1.shape,i,dataset_pythia.shape[1])

        # Concatenate the different arrays
        data_pythia = np.concatenate((data_pythia,dataset_pythia1),axis=0)

# Remove first entry again
# otherwise initialize data_jewel = None, and fill it with dataset_jewel1 if it is empty (or for i=0)
data_pythia_final = data_pythia[1:]
print('final data shape for Pythia: data and labels')
print(data_pythia_final.shape)

# Create labels: Pythia 0, Jewel 1
y_pythia = np.zeros(data_pythia_final.shape[0])
print(y_pythia.shape)

# ---------------------------------------------------------
# Second: read in Jewel data files and do some reformatting
# ---------------------------------------------------------

# Initialize data structure
data_jewel = np.zeros((1,800,4))
#print(data_jewel.shape)

print('Reading in Jewel h5 files')
for i,file in enumerate(os.listdir('/rstorage/ml/egml/jewel/')):
    if file.endswith('.h5') and i<nh5:
        with h5py.File(os.path.join('/rstorage/ml/egml/jewel/', file),'r') as hdf:
            # jetsubs_cons are the 4-vectors of particles inside the jets
            dataset_jewel0 = hdf.get('jetsubs_cons')
            
            # Change format: jet, particle, 4-vector
            dataset_jewel = np.transpose(np.array(dataset_jewel0),(0,2,1)) 
            
            # Reorder 4-vector from (px,py,pz,E) to (E,px,py,pz)
            dataset_jewel[:,:,[0,1,2,3]] = dataset_jewel[:,:,[3,0,1,2]]
            
            # Replace nan with zero
            where_are_nan = np.isnan(dataset_jewel)
            dataset_jewel[where_are_nan] = 0.
            
            # Zero pad such that they all have the same shape (.., 800, 4)
            dataset_jewel1 = np.pad(dataset_jewel,((0,0),(0,800-dataset_jewel.shape[1]),(0,0)),'constant',constant_values=(0))
            print(dataset_jewel1.shape,i,dataset_jewel.shape[1])

        # Concatenate the different arrays
        data_jewel = np.concatenate((data_jewel,dataset_jewel1),axis=0)

# Remove first entry again
data_jewel_final = data_jewel[1:]
print('final data shape for Jewel: data and labels')
print(data_jewel_final.shape)

# Create labels: Pythia 0, Jewel 1
y_jewel = np.ones(data_jewel_final.shape[0])
print(y_jewel.shape)

# -----------------------------------------------------------------------------
# Third: Write labeled, shuffled and equal sample size data to a single h5 file
# -----------------------------------------------------------------------------

# Write new data to file: Pythia & Jewel + labels
# 1. Get the number of jets, such that we get a data set with an equal number of pp and AA jets
number_of_jets = min(data_pythia_final.shape[0],data_jewel_final.shape[0])
print(number_of_jets)

# 2. Concatenate Pythia and Jewel data and labels
data_ppAA_jets = np.concatenate((data_pythia_final[:number_of_jets],data_jewel_final[:number_of_jets]),axis=0)
labels_ppAA_jets = np.concatenate((y_pythia[:number_of_jets],y_jewel[:number_of_jets]),axis=0)
print(data_ppAA_jets.shape,labels_ppAA_jets.shape)

# 3. Shuffle data set # check again!
idx = np.random.permutation(len(labels_ppAA_jets))
data_ppAA_jets1, labels_ppAA_jets1 = data_ppAA_jets[idx], labels_ppAA_jets[idx]

# 4. Write result to a new hdf5 file
print('writing to new file: data_ppAA_PythiaJewel_jets.h5')
with h5py.File('data_ppAA_PythiaJewel_jets.h5','w') as hdf:
    hdf.create_dataset('data_ppAA_PythiaJewel',data = data_ppAA_jets1)
    hdf.create_dataset('labels_ppAA_PythiaJewel',data = labels_ppAA_jets1)
    
print('done.')