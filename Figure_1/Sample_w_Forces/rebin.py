# @author Maria Lesniewski
# @desc  This is a handy file to rebin trajectory data

import numpy as np

a = np.loadtxt("sample_trajectory.xy")
domain = [-3, 16]
bins = 5000

x = np.linspace(domain[0], domain[1], bins)

idx = np.digitize(a[:,1], x[:-1]) - 1

hist, _ = np.histogram(idx, bins=np.arange(0, len(x) + 1))

np.savetxt("rebinned_data_5000", np.column_stack((x, hist)))
