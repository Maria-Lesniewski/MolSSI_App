# @author Maria Lesniewski 
# @desc This script takes an input unnormalized distribution [stored as an xy file]
#       and calculates the gradients of the distribution. 
import sys
import numpy as np

def calculate_negative_log_y(y):
    # Check if y is non-zero to avoid division by zero
    nonzero_indices = np.nonzero(y)
    negative_log_y = np.zeros_like(y)
    negative_log_y[nonzero_indices] = -np.log(y[nonzero_indices])
    return negative_log_y

def calculate_negative_gradient(x, y):
    # Check if y is non-zero to avoid division by zero
    #nonzero_indices = np.nonzero(y)
    negative_gradient = np.zeros_like(y)
    #negative_gradient[nonzero_indices] = np.gradient(y[nonzero_indices], x[nonzero_indices])
    negative_gradient = -np.gradient(y, x)
    return negative_gradient

# Step 1: Read the input xy file from command-line argument
if len(sys.argv) < 2:
    print("Please provide the input filename as a command-line argument.")
    sys.exit(1)

input_file = sys.argv[1]
data = np.loadtxt(input_file)
x = data[:, 0]
y = data[:, 1]

# Step 2: Calculate the negative log of y and write to "potential.xy"
negative_log_y = calculate_negative_log_y(y)
potential_data = np.column_stack((x, negative_log_y))
np.savetxt("potential.xy", potential_data)

# Step 3: Calculate the negative gradient and write to "force.xy"
force_data = np.column_stack((x, calculate_negative_gradient(x, negative_log_y)))
np.savetxt("force.xy", force_data)

