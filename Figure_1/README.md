These scripts were used to create figure_1 in the proposal submitted

Training_Data: 
  Scripts for sampling from an analytic curve p(x) by grid to get a "fully sampled" training set
  dropping some portion of "hard to sample" points based on the analytic probability to create a "partially sampled" set 
  adding gaussian distributed noise to the partial data set to create a "perturbed" set
Learn_Forces:
  Script to take the derivative of the log of the distributions passed [using log 0 = 0 to get around award domain issues]
Sample_w_Forces:
  Script for uncorrected langevin monte carlo sampling using learned forces
Gaussian_Deconvolution:
  Script for deconvolution of sampled data using initial guess from perturbed distribution [= modes of true distro by convolution] & loglikelihood in 1D
