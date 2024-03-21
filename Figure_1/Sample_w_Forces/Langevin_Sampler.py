# @author Maria Lesniewski
# This is a script that takes an .xy file for force along with an initial condition
# We linearly interpolate between force points when not on a direct knotpoint inside the force.xy file
# We save a trajectory of the x values explored and bin them at the end to get a distribution
# Currently I think I am interested in 5000 points between [-3,16], so we'll follow suit with our grid 
# sampling method and attempt 100,000 per point I am interested in. 

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import sys
import math


# Function to linearly interpolate force information between knot points
def interpolate_force(Fx, f_domain, x0):
    print(x0)    
    f_idx = np.digitize(float(x0), f_domain) - 1 # this is the bin we'd read from the table
    print(f_idx)
    if f_idx == -1: 
        return -(x0-3)
    if f_idx >= len(f_domain) - 1:
        return -1/3*(x0-10)
    else:
        idx2 = f_idx + 1 # But the index we read floors the x value so we want to interp between

        slope = (Fx[idx2] - Fx[f_idx]) / (f_domain[idx2] - f_domain[f_idx])
        dx = x0 - float(f_domain[f_idx])

        return Fx[f_idx] + slope*dx


##########MAIN##################
# Script Usage
if len(sys.argv) < 3:
    print("Please provide the .xy force file and the initial condition for sampling")
    sys.exit(1)


force_file = sys.argv[1]
x0 = float(sys.argv[2])

# Sample Set Up
forces = np.loadtxt(force_file)
Fx = forces[:,1]
f_domain = forces[:,0]

target_domain = [-3,16]

bins_at_end = 500 
samples_per_bin = 100000 
num_frames_saved = bins_at_end * samples_per_bin
data_drawn = []
time = []
# Langevin sampling and trajectory saving
step_size = 0.005 
langevin_step_per_frame_saved = 10
xnew = 0

for i in range(num_frames_saved):
    for j in range(langevin_step_per_frame_saved):
        xnew = x0 + step_size * interpolate_force(Fx, f_domain, x0)  + math.sqrt(2*step_size) * np.random.normal(0,1)
        x0 = xnew
    data_drawn.append(xnew)
    time.append(i)

# Save Trajectory 
traj = np.column_stack((time, data_drawn))
np.savetxt("sample_trajectory.xy", traj)

# Bin instances to get distro
x_bins = np.linspace(target_domain[0], target_domain[1], bins_at_end)
sorted_indices = np.digitize(data_drawn, x_bins) - 1
hist, _ = np.histogram(sorted_indices, bins=np.arange(0, len(x_bins) +1))
distro = np.column_stack((x_bins, hist))

np.savetxt("trajectory_distro.xy", distro)
