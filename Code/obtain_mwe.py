import numpy as np
import pandas as pd
import basis_set_exchange as bse


def GF(alpha, r, R_A):
    return 0

#Read NMO, NCGFC and Basis from ConfigFile
CFile_config_numbers= pd.read_csv("ConfigFile.txt", sep = " ", skiprows = lambda x: x not in [5,6,7])
NMO = CFile_config_numbers["NMO"].tolist()[0]
NCGF = CFile_config_numbers["NCGFC"].tolist()[0]
Basis = CFile_config_numbers["Basis"].tolist()[0]

#Read coordinates from ConfigFile
CFile_atom_coordinates= pd.read_csv("ConfigFile.txt", sep = " ", skiprows = lambda x: x not in [8,9,10,11])
Atom_1_xyz = CFile_atom_coordinates["C"].tolist()
Atom_2_xyz = CFile_atom_coordinates["O"].tolist()


#Read Orbitals information from ConfigFile
CFile_Orbitals = pd.read_csv("ConfigFile.txt", sep = " ", skiprows = 12)
MO1 = CFile_Orbitals["MO1"].tolist()
MO2 = CFile_Orbitals["MO2"].tolist()




initial_basis = bse.get_basis(Basis, elements = 1 ,fmt='gaussian94', header=False)
initial_basis_dict = bse.get_basis(Basis, elements = 1 , header=False)

coeff_CGF = initial_basis_dict['elements']['1']['electron_shells'][0]['coefficients']
expon_CGF = initial_basis_dict['elements']['1']['electron_shells'][0]['exponents']


#print(initial_basis_dict)

#for key in Coeff_d:
 # print(key)
