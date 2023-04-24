import numpy as np


#We define a function for computing the initial Density Matrix of any system
#For that we need a list A with all the states and AdjA its duals

def Compute_DOP(A, AdjA):
    DOP_row_size = len(A)
    DOP = np.zeros((DOP_row_size,DOP_row_size), dtype = np.float64)
    for i,j in zip(A, AdjA):
        DOP += np.outer(i,j)
    return DOP

A = [np.sqrt(1/2),np.sqrt(1/2)]


print(Compute_DOP(A,A))
