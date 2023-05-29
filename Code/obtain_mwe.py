import numpy as np
import pandas as pd
import basis_set_exchange as bse


initial_basis = bse.get_basis('STO-3G', elements = 1 ,fmt='gaussian94', header=False)
initial_basis_dict = bse.get_basis('STO-3G', elements = 1 , header=False)

coeff_CGF = initial_basis_dict['elements']['1']['electron_shells'][0]['coefficients']
expon_CGF = initial_basis_dict['elements']['1']['electron_shells'][0]['exponents']


#print(initial_basis_dict)
print(coeff_CGF)
print(expon_CGF)

#for key in Coeff_d:
 # print(key)
