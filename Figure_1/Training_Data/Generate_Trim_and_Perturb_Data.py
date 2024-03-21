# @author Maria Lesniewski
# @desc this script uses inverse CDF sampling to generate data from a bimodal distribution
#       it then drops the least sampled data points to simulate lacking data 
#       This partial data is saved, and then gaussian noise is thrown on each data point remaining
#       This new noisy distribution is then saved


import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
def f(x):
    # Edit this function according to your specific distribution
    return 5 * np.exp(-0.5 * (x - 3) ** 2) + np.exp(-0.5 * 1/3 * (x - 10) ** 2)

def dist(x):
    return f(x) / normalization_constant


# 1) Determin some function and the domain its not zero
domain = [-3, 16]  # Specify the domain for x

# 2) Obtain expression for normalization constant so that we have access to probabilities
normalization_constant, _ = integrate.quad(f, domain[0], domain[1])

# 3) Generate like 100k samples per gridpoint I am interested in. e.g. I am interested in 5k gridpoints from 0 to 15
end_length = 5000
fraction = int(end_length/3) # This is the fraction of the domain we'll throw away
samples_per_bin = 100000
num_objects = end_length * samples_per_bin

# 4) Getting gridpoint samples
uniform_probabilities = np.random.uniform(0,1, num_objects)
x_samples = np.linspace(domain[0], domain[1], num_objects)
x_bins = np.linspace(domain[0], domain[1], end_length)

accepted_samples = x_samples[uniform_probabilities <= dist(x_samples)]

# 5) Binning to find sampled distribution 
#hist, _ = np.histogram(accepted_samples, bins=np.arange(0, len(x_bins) + 1))
x_coarse_index = np.digitize(accepted_samples, x_bins[:-1]) - 1 # Obtaining the bin index that the samples were thrown in
hist, _ = np.histogram(x_coarse_index, bins=np.arange(0, len(x_bins) + 1))
histogram_data = np.column_stack((x_bins, hist))
np.savetxt("full_unperturbed_distribution.xy", histogram_data)

# 6) Now we need to drop the lowest sampled points to model a disparity in sampling
lowest_count_idx = np.argsort(hist)[:fraction]
x_w_lowest_counts = x_bins[lowest_count_idx]
filtered_bins = x_bins[~np.isin(x_bins, x_w_lowest_counts)]
filtered_hist = hist[~np.isin(np.arange(len(hist)), lowest_count_idx)]

filtered_data = np.column_stack((filtered_bins, filtered_hist))
np.savetxt("partial_unperturbed_distribution.xy", filtered_data)

# 7) Now the tricky part, we've got to drop the individual samples from our accepted data and sanity check
trimmed_accept = accepted_samples[~np.isin(x_coarse_index, lowest_count_idx)]
trimmed_coarse_indices = np.digitize(trimmed_accept, x_bins[:-1]) - 1
sanity_hist, _ = np.histogram(trimmed_coarse_indices, bins=np.arange(0, len(x_bins) +1))
sanity_data = np.column_stack((x_bins, sanity_hist))
np.savetxt("sanity_check_partial.xy", sanity_data)


# 8) Now that we dropped the objects in fine space we need to add noise to the accepted samples from our training data
perturbed_partial_data = trimmed_accept + np.random.normal(0,1, len(trimmed_accept)) # Add noise with mean zerio and dev 1
perturbed_coarse_idx = np.digitize(perturbed_partial_data, x_bins[:-1]) -1
noise_hist, _ = np.histogram(perturbed_coarse_idx, bins=np.arange(0, len(x_bins)+1))
final_data = np.column_stack((x_bins, noise_hist))
np.savetxt("perturbed_partial_distribution.xy", final_data)
