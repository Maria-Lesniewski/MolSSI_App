# @author Maria Lesniewski 
# @desc Usage : arg1 = a perturbed distribution sampled by gradients, formatted as an .xy file
#       This file deconvolutes a bimodal distribution by minimizing the KL divergence between the observed distribution
#       and the output file convolved with known noise
#       e.g. if h = f*g and I know h and g, minimize KL(h, f_output*g) 


import numpy as np
from scipy.optimize import minimize
from scipy import signal
import sys
from scipy.optimize import Bounds


# Noise Function
def gauss(x):
    return np.exp(-0.5*(x**2))

# Model 
def gauss_combo(domain, params):
    return params[0]*np.exp(-0.5/(params[1]**2)*(domain - np.ones_like(domain)*params[2])**2) + params[3]*np.exp(-0.5/(params[4]**2)*(domain - np.ones_like(domain)*params[5])**2)

# Kullbach Liebler without log0
def safe_KL(pu,qu):
    p = pu/(sum(pu))
    q = qu/(sum(qu))
    p_zero = (p <= 0)
    q_zero = (q <= 0)
    log_p_div_q = np.log(np.divide(p,q, out=np.ones_like(p), where=~(p_zero | q_zero )))
    i = p * log_p_div_q
    mask = log_p_div_q <= 0
    return sum(i)

# Construct model and take KL div
def min_objective(gauss_params, observed_domain, observed_prob):
    #a1, sig1, mu1, a2, sig2, mu2 = paramsorder
    model = gauss_combo(observed_domain, gauss_params)
    dx = observed_domain[1]-observed_domain[0]
    kernel_domain = np.linspace(-5,5, int(10/(dx)+1))
    kernel = gauss(kernel_domain)
     
  # Left KL Divergence [guess is consistent with prior]
    model_convolved = signal.convolve(model, kernel, 'same')
    b = safe_KL(observed_prob/(sum(observed_prob*dx)), model_convolved/(sum(model_convolved*dx)))
    print(b)
    return b 


#### MAIN ########
convolved_data = np.loadtxt(sys.argv[1]) # Pass the noisy distro.xy

x = convolved_data[:,0]
noisey_dist = convolved_data[:,1]

# Make a guess using properties of training data before perturbation [mode centers will be same after perturb]
parameter_init_guess = np.array([1,1,3,1,1,10])

# Minimize objective & save results
result = minimize(min_objective, parameter_init_guess, args=(x, noisey_dist), tol=0.001, bounds=Bounds(0.1,100))

print(result.x)
print(result.success)
print(result.message)
final = gauss_combo(x,result.x)
np.savetxt("Gaussians_Deconvolved.out", np.column_stack((x, final)))

#check convergence looks reasonable by replotting h = f*g
kernel_domain = np.linspace(-5,5, int(10/(x[1]-x[0])+1))
kernel = gauss(kernel_domain)
check = signal.convolve(final, kernel, 'same')
np.savetxt("check.out", np.column_stack((x,check)))
