import numpy as np
import basis_set_exchange as bse

#We define a function for specifying the basis in which we are going to expand the space orbitals
# flag arguments indicates if we will use lowdin or other package:

def Define_Basis(Basis_Set, flag):
    basis_elements_array = []
    if (flag == "gauss"):
        bs_dict = bse.get_basis(Basis_Set, elements = ['C',8], header=False)
        return bs_dict
    if (flag == "lowdin"):
        return 0
    
#We give the Coefficients(C) of the QC method for the elements of the basis thus creating the
# spacial orbitals

def Construct_Orbitals(C):




#We define a function for computing the initial Density Matrix of any system
#For that we need a list A with all the states and AdjA its duals


def Compute_DOP(A, AdjA):
    DOP_row_size = len(A)
    DOP = np.zeros((DOP_row_size,DOP_row_size), dtype = np.float64)
    for i,j in zip(A, AdjA):
        DOP += np.outer(i,j)
    return DOP

#A = [np.sqrt(1/2),np.sqrt(1/2)]
#print(Compute_DOP(A,A))

print(Define_Basis("6-31G*", "test"))
