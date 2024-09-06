import numpy as np
import h5py

f = h5py.File('EX1\Exercise1.mat', 'r')
data = f.get('Channel')
data = np.array(data)

print(data)