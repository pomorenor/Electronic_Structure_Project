import numpy as np
from scipy import special
import math
import matplotlib.pyplot as plt


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
##  Now we write the necessary functions for computing the integrals##
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


def F_0(t):
    if (t<1e-6):
        return 1.0-t/3.0
    else:
        return 0.5*(np.pi/t)**0.5*special.erf(t**(0.5))


def nuclear_attraction_integral(alpha, R_A, beta, R_B, nuclei_Z, nuclei_coords):

    R_C = nuclei_coords
    R_p =(alpha*R_A + beta*R_B)/(alpha + beta)
    normalization = (2.0*alpha/np.pi)**0.75*(2.0*beta/np.pi)**0.75
    exp_part = np.exp(((-alpha*beta)/(alpha + beta))*np.abs(R_A-R_B)**2)
    F0 = F_0((alpha+beta)*np.abs(R_p - R_C)**2)
    factor = -2.0*np.pi/(alpha + beta)
    return normalization*factor*nuclei_Z*exp_part*F0

def four_centers_integral(alpha, R_A, beta, R_B,  gamma, R_C, delta, R_D):

    R_P = (alpha*R_A + beta*R_B)/(alpha + beta)
    R_Q = (gamma*R_C + delta*R_D)/(gamma + delta)
    normalization = (2.0*alpha/np.pi)**0.75*(2.0*beta/np.pi)**0.75*(2.0*gamma/np.pi)**0.75*(2.0*delta/np.pi)**0.75
    factor = 2*np.pi**2.5/((alpha + beta)*(gamma + delta)*(alpha + beta + gamma + delta)**0.5)
    exp_part = np.exp(((-alpha*beta)/(alpha + beta))*np.abs(R_A-R_B)**2 -  ((gamma*delta)/(gamma + delta))*np.abs(R_C-R_D)**2)
    F0 = F_0(((alpha + beta)*(gamma + delta)/(alpha + beta + gamma + delta))*np.abs(R_P -R_Q)**2)
    return normalization*factor*exp_part*F0

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

def four_centers_integral_with_CGF(mu, nu, pi, ro, contraction_exponents_of_orbitals, centros, contraction_length,contraction_coefficients_of_orbitals):
    Fcenter = 0.0
    for p in range(0, contraction_length):
        for q in range(0, contraction_length):
            for r in range(0, contraction_length):
                for s in range(0, contraction_length):
                    Fcenter += contraction_coefficients_of_orbitals[mu][p]*contraction_coefficients_of_orbitals[nu][q]*contraction_coefficients_of_orbitals[pi][r]*contraction_coefficients_of_orbitals[ro][s]*four_centers_integral(contraction_exponents_of_orbitals[mu][p], centros[mu], contraction_exponents_of_orbitals[nu][q], centros[nu], contraction_exponents_of_orbitals[pi][r],centros[pi], contraction_exponents_of_orbitals[ro][s], centros[ro])
    return Fcenter


######################################################################################################
##    Functions for constructing the matrices that will change during SCG                           ##
##                                                                                                  ##
######################################################################################################

def Compute_P_matrix(C, basis_size, num_electrons):
    P = np.zeros((2,2))
    for mu in range (0, basis_size):
        for nu in range(0, basis_size):
            for a in range(0, int(num_electrons/2)):
                P[mu][nu] = 2.0*C[mu][a]*C[nu][a]
    return P

def Compute_G_Matrix(P_matrix, basis_size, contraction_exponents_of_orbitals, centros, contraction_length,contraction_coefficients_of_orbitals):
    G= np.zeros((2,2))
    for mu in range(0, basis_size):
        for nu in range(0, basis_size):
            for lambd in range(0, basis_size):
                for sigma in range(0, basis_size):
                    G[mu][nu] += P_matrix[lambd][sigma]*(four_centers_integral_with_CGF(mu, nu,sigma, lambd, contraction_exponents_of_orbitals, centros, contraction_length,contraction_coefficients_of_orbitals)- 0.5*four_centers_integral_with_CGF(mu, lambd, sigma, nu, contraction_exponents_of_orbitals, centros, contraction_length,contraction_coefficients_of_orbitals))
    return G


#def four_centers_integral_with_CGF(mu, nu, pi, ro, contraction_exponents_of_orbitals, centros, contraction_length,contraction_coefficients_of_orbitals):


def Compute_Fock_matrix(HCore, G_matrix, basis_size):
    F = np.zeros((2,2))
    for mu in range(0,basis_size):
        for nu in range(0,basis_size):
            F[mu][nu] = HCore[mu][nu] + G_matrix[mu][nu]
    return F



########################################################################################################

orbitals_coefficients = [[],[]]
orbital_exponents = [[],[]]
CENTERS = [[],[]]

R_A = 0.0
R_B = 1.4

scaled_exponents = scale_exponents(3, contraction_exponents_zeta_1, 1.69)
#scaled_exponents = scale_exponents(3, contraction_exponents_zeta_1, 1.69)



construct_initial_orbitals = CGF_of_each_orbital(orbitals_coefficients, contraction_coefficients[2])
construct_initial_orbitals_exponents = CGF_exponents_of_each_orbital(orbital_exponents, scaled_exponents)
array_of_centers = Centers(CENTERS, R_A, R_B)

#S_12 = overlap_integral_with_CGF(1,1,construct_initial_orbitals_exponents, array_of_centers, 3,construct_initial_orbitals)

#T_11 = kinetic_integral_with_CGF(1,1,construct_initial_orbitals_exponents, array_of_centers, 3,construct_initial_orbitals)
#V1_11 = nuclear_attraction_integral_with_CGF(1,1,construct_initial_orbitals_exponents, array_of_centers, 3,construct_initial_orbitals,1,0)

#FCENTER_11 = four_centers_integral_with_CGF(1,0,1,0, construct_initial_orbitals_exponents, array_of_centers,3,construct_initial_orbitals)

#print(FCENTER_11)

###################################################################
## We now construct the matrices that will not change upon SC    ##
###################################################################

T = np.zeros((2,2))
V_1 = np.zeros((2,2))
V_2 = np.zeros((2,2))


for i in range(0,2):
    for j in range(0,2):
        T[i][j] = kinetic_integral_with_CGF(i,j,construct_initial_orbitals_exponents, array_of_centers, 3,construct_initial_orbitals)
        V_1[i][j] = nuclear_attraction_integral_with_CGF(i,j,construct_initial_orbitals_exponents, array_of_centers, 3,construct_initial_orbitals,2,0)
        V_2[i][j] = nuclear_attraction_integral_with_CGF(i,j,construct_initial_orbitals_exponents, array_of_centers, 3,construct_initial_orbitals,1,1.4)

#####################################################################
## We now construct the matrices we will need Hcore matrix         ##
#####################################################################

Hcore = np.zeros((2,2))
S_munu = np.zeros((2,2))


for i in range(0,2):
    for j in range(0,2):
        Hcore[i][j] = T[i][j] + V_1[i][j] + V_2[i][j]



for i in range (0,2):
    for j in range(0,2):
        S_munu[i][j] = overlap_integral_with_CGF(i,j,construct_initial_orbitals_exponents, array_of_centers, 3,construct_initial_orbitals)


#############################################################################
## We know diagonalize the overlap matrix and obtain X through symmetrical ##
##diagonalization                                                          ##
#############################################################################

eigenvalues, U = np.linalg.eig(S_munu)
s = np.dot(U.T,np.dot(S_munu,U))
s_half_minus = np.diag(np.diagonal(s)**-0.5)
X = np.dot(U,np.dot(s_half_minus,U.T))


##############################################################################
## We know start the SCF procedure                                         ###
##############################################################################

P = Compute_P_matrix(np.zeros((2,2)), 2, 2)

#We set the initial P equal to HCore
#for i in range(0,2):
#    for j in range(0,2):
#       P[i][j] = Hcore[i][j]


##########################################
## From here should begin the loop      ##
##########################################

num_iter = 200
ii = 0
tolerance = 1e-11


while(ii < num_iter):
    ii +=1
    print("Iteration number: ", ii)
    G = Compute_G_Matrix(P, 2, construct_initial_orbitals_exponents, array_of_centers,3,construct_initial_orbitals)

    F = Compute_Fock_matrix(Hcore, G, 2)

    Energy = np.sum(0.5*P*(Hcore+F))


    F_prime = np.dot(X.T, np.dot(F,X))

    epsilon, C_prime = np.linalg.eig(F_prime)



    C = np.dot(X,C_prime)

    OldP = np.array(P)
    P = np.zeros([2,2])

    P = Compute_P_matrix(C, 2, 2)


    Delta = (P-OldP)
    Delta = np.sqrt(np.sum(Delta**2)/4.0)

    print("Delta",Delta)
    print("Energy: ", Energy)

    if (Delta<tolerance):
        break



#print(S_munu)
#print(Hcore)
#print(Hcore)
#print(T)
#print(V_1)
#print(V_2)

#print(construct_initial_orbitals_exponents)
#print(construct_initial_orbitals)
#print(array_of_centers)
#print(S_12)
#print(T_11)
#print(V1_11)
######################################################
## Now we write the code to compute the integrals ###
######################################################


####################################################
## We write the code to ensamble de matrices #######
###################################################
