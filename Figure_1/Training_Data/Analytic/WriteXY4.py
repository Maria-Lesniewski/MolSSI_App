# @Author Maria Lesniewski
# This is a handy dandy script for writing a function you want to plot
# to an xy file, simply alter my_func and domain as you want


import numpy as np

def write_fx(f, x, filename):
    y = f(x)
    out = np.column_stack((x,y))
    print(out)
    np.savetxt(filename, out, header="x y")

def my_fun(x):
    # Edit this as you want to print the function
    ones = np.ones_like(x)
    return np.log(5*np.exp(-0.5*(x-3*ones)**2) + np.exp(-0.5*1/3*(x-10*ones)**2))

#####Main######################################
filename = "unnorm_pdf.xy"
domain = [-3,16]
numpoints = 5000

x = np.linspace(domain[0], domain[1], numpoints)

write_fx(my_fun, x, filename)
