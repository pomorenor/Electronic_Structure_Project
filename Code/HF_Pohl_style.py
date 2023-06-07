import numpy as np
from scipy import special
import math



######################################################
## First we declare the things associated with the ##
## gaussian functions                              ##
#####################################################


contraction_coefficients = [[1.0],[0.678914,0.430129],[0.444635,0.535328,0.154329]]
contraction_exponents_zeta_1 = [[0.270950],[0.152623,0.851819],[0.109818,0.405771,2.22766]]
contraction_length = [1,2,3]


def scale_exponents(STO_NG,contraction_exponents_zeta_1,new_zeta):

    def scale(alpha_zeta_uno, new_zeta):
        new_alpha = alpha_zeta_uno*new_zeta**2
        return new_alpha

    if (STO_NG == 1):
        return [scale(i,new_zeta) for i in contraction_exponents_zeta_1[0]]
    elif (STO_NG == 2):
        return [scale(i,new_zeta) for i in contraction_exponents_zeta_1[1]]
    elif (STO_NG == 3):
        return [scale(i,new_zeta) for i in contraction_exponents_zeta_1[2]]


##########################################################################

def CGF_of_each_orbital(orbitals_coefficients, contraction_coefficients):
    orbitals_coefficients[0] = contraction_coefficients
    orbitals_coefficients[1] = contraction_coefficients

    return orbitals_coefficients

def CGF_exponents_of_each_orbital(orbitals_exponents, contraction_scaled_exponents):
    orbitals_exponents[0] = contraction_scaled_exponents
    orbitals_exponents[1] = contraction_scaled_exponents

    return orbitals_exponents

def Centers(Center_array, R_A, R_B):
    Center_array[0] = R_A
    Center_array[1] = R_B

    return Center_array
######################################################################

def normalized_two_centers_gaussian_integral(alpha, R_A, beta, R_B):
    normalization = ((2*alpha)/np.pi)**0.75*((2*beta)/np.pi)**0.75
    K = (np.pi/(alpha + beta))**1.5
    integral_value = np.exp(((-alpha*beta)/(alpha + beta))*np.abs(R_A - R_B)**2)
    return normalization*K*integral_value

def kinetic_integral(alpha, R_A, beta, R_B):
    normalization = ((2.0*alpha)/np.pi)**0.75*((2.0*beta)/np.pi)**0.75
    integral_first_part = ((alpha*beta)/(alpha + beta))*(3.0-((2.0*alpha*beta*np.abs(R_A-R_B)**2)/(alpha+beta)))*(np.pi/(alpha+beta))**1.5
    integral_exp =  np.exp(((-alpha*beta)/(alpha+beta))*np.abs(R_A-R_B)**2)

    return normalization*integral_first_part*integral_exp

def nuclear_attraction_integral(alpha, R_A, beta, R_B, nuclei_Z, nuclei_coords):

    R_p = (alpha*R_A + beta*R_B)/(alpha + beta)
    normalization = ((2.0*alpha)/np.pi)**0.75*((2.0*beta)/np.pi)**0.75
    integral_exp =  np.exp(((-alpha*beta)/(alpha+beta))*np.abs(R_A-R_B)**2)
    erf_part = special.erf((alpha + beta)*np.abs(R_p - nuclei_coords)**2)

    return ((-2.0*np.pi*nuclei_Z)/(alpha + beta))*normalization*integral_exp*erf_part


def kinetic_integral_with_CGF(mu, nu, contraction_exponents_of_orbitals, centros, contraction_length,contraction_coefficients_of_orbitals):
    T = 0.0
    for p in range(0, contraction_length):
        for q in range(0,contraction_length):
            T += contraction_coefficients_of_orbitals[mu][p]*contraction_coefficients_of_orbitals[nu][q]*kinetic_integral(contraction_exponents_of_orbitals[mu][p], centros[mu], contraction_exponents_of_orbitals[nu][q], centros[nu])
    return T

def overlap_integral_with_CGF(mu, nu, contraction_exponents_of_orbitals, centros, contraction_length,contraction_coefficients_of_orbitals):
    S = 0.0
    for p in range(0,contraction_length):
        for q in range(0,contraction_length):
            S += contraction_coefficients_of_orbitals[mu][p]*contraction_coefficients_of_orbitals[nu][q]*normalized_two_centers_gaussian_integral(contraction_exponents_of_orbitals[mu][p], centros[mu], contraction_exponents_of_orbitals[nu][q], centros[nu])
    return S

def nuclear_attraction_integral_with_CGF(mu, nu, contraction_exponents_of_orbitals, centros, contraction_length,contraction_coefficients_of_orbitals, nucleus_z, nuclei_coords):
    V = 0.0
    for p in range(0, contraction_length):
        for q in range(0, contraction_length):
            V += contraction_coefficients_of_orbitals[mu][p]*contraction_coefficients_of_orbitals[nu][q]*nuclear_attraction_integral(contraction_exponents_of_orbitals[mu][p], centros[mu], contraction_exponents_of_orbitals[nu][q], centros[nu], nucleus_z, nuclei_coords)
    return V

orbitals_coefficients = [[],[]]
orbital_exponents = [[],[]]
CENTERS = [[],[]]

R_A = 0.0
R_B = 1.4

scaled_exponents = scale_exponents(3, contraction_exponents_zeta_1, 1.24)

construct_initial_orbitals = CGF_of_each_orbital(orbitals_coefficients, contraction_coefficients[2])
construct_initial_orbitals_exponents = CGF_exponents_of_each_orbital(orbital_exponents, scaled_exponents)
array_of_centers = Centers(CENTERS, R_A, R_B)

S_12 = overlap_integral_with_CGF(1,1,construct_initial_orbitals_exponents, array_of_centers, 3,construct_initial_orbitals)

T_11 = kinetic_integral_with_CGF(1,1,construct_initial_orbitals_exponents, array_of_centers, 3,construct_initial_orbitals)
V1_11 = nuclear_attraction_integral_with_CGF(1,1,construct_initial_orbitals_exponents, array_of_centers, 3,construct_initial_orbitals,2,1.4)

print(scaled_exponents)
print(construct_initial_orbitals_exponents)
print(construct_initial_orbitals)
print(array_of_centers)
print(S_12)
print(T_11)
print(V1_11)
######################################################
## Now we write the code to compute the integrals ###
######################################################



####################################################
## We write the code to ensamble de matrices #######
###################################################
