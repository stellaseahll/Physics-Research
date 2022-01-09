# cd Documents/CQT/Quantum Thermodynamics/Rotor

# Imports all qutip functions
from qutip import *
# For sqrt, exp and other functions
from math import exp
from math import pi
from math import sin
from math import cos
from math import sqrt
import scipy.integrate as integrate
import numpy as np

#For tic-toc
import time

# For plots
print(time.ctime())
tic = time.clock()

# NOTE: hbar=1 throughout
#       Indices 0,1,2,...,m_min+m_max correspond to -m_min,-m_min+1,...,0,..., m_max-1,m_max

##############################
## Parse Arguments
#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-th", "--th", dest="th_value", type=float, metavar='<float>', default=2.0, help="TH value (default=5.0)")
parser.add_argument("-tc", "--tc", dest="tc_value", type=float, metavar='<float>', default=1.0, help="TC value (default=1.0)")
parser.add_argument("-tw", "--tw", dest="tw_value", type=float, metavar='<float>', default=8.0, help="TR value (default=2.0)")
parser.add_argument("-omegaC", "--omegaC", dest="omegaC_value", type=float, metavar='<float>', default=1.0, help="W1 value (default=1.0)")
parser.add_argument("-omegaW", "--omegaW", dest="omegaW_value", type=float, metavar='<float>', default=4.0, help="W2 value (default=5.0)")
parser.add_argument("-k", "--k", dest="k_value", type=float, metavar='<float>', default=0.001, help="k value (default=0.001)")
parser.add_argument("-g", "--g", dest="g_value", type=float, metavar='<float>', default=0.001, help="g value (default=0.001)")
parser.add_argument("-J", "--J", dest="J_value", type=int, metavar='<int>', default=1, help="J value (default=2.0)")

args = parser.parse_args()

def meanOccupation(E,Tr):
    if (Tr != 0):
        x = 1.0/(exp(E/Tr)-1)
        return x
    x = 0
    return x

#------PARAMETERS----------

# Dimensions for the Modes
dC = 4
dH = 4
dW = 12

# Operators for Harmonic Oscillators
NC = num(dC)
NH = num(dH)
NW = num(dW)
aC = destroy(dC)
aH = destroy(dH)
aW = destroy(dW)
adC = create(dC)
adH = create(dH)
adW = create(dW)

# Energy of the qubits
omegaC = args.omegaC_value
omegaW = args.omegaW_value
omegaH = omegaW+omegaC

TC = args.tc_value #Qubit 1
TH = args.th_value  #Qubit 2

# Coupling strength
k1 = args.k_value
k2 = args.k_value
k3 = args.k_value

# System Hamiltonian
hC = omegaC*tensor(NC,qeye(dH),qeye(dW))
hH = omegaH*tensor(qeye(dC),NH,qeye(dW)) 
hW = omegaW*tensor(qeye(dC),qeye(dH),NW)
h2 = 0.5*omegaR*tensor(qeye(d),sigmaz)
H0 = h1 + h2


# Interaction Hamiltonian
g = args.g_value
H0 = h1 + h2#+ Hint
# Total Hamiltonian
# H = H0 + Hint

TR = np.zeros(51)
E = np.zeros((51,9))

for j in range(51):
	TR[j] = j*0.25
# Temperature of baths


	# Mean Occupation of bath
	N1 = meanOccupation(omegaC/(d-1.0),TC)
	N2 = meanOccupation(omegaR/(d-1.0),TR[j])
	N3 = meanOccupation(omegaH/(d-1.0),TH)

	# # Thermal state of the qubits
	# tau1 = np.exp(-tmp*0.5*omegaC/(TC*1.0))
	# tau2 = np.exp(-tmp*0.5*omegaR/(TR*1.0))
	# tau3 = np.exp(-tmp*0.5*omegaH/(TH[j]*1.0))
	# tau1 = np.diag(tau1/np.sum(tau1))
	# tau2 = np.diag(tau2/np.sum(tau1))
	# tau3 = np.diag(tau3/np.sum(tau1))

	# Lindblad Dissipators

	C1=sqrt(k1*(N1+1))*tensor(sigmam,qeye(d))
	C1d=sqrt(k1*(N1))*tensor(sigmap,qeye(d))
	C2=sqrt(k2*(N2+1))*tensor(qeye(d),sigmam)
	C2d=sqrt(k2*(N2))*tensor(qeye(d),sigmap)
	C3=sqrt(k3*(N3+1))*tensor(sigmam,sigmam)
	C3d=sqrt(k3*(N3))*tensor(sigmap,sigmap)

	Q1 = C1.dag()*H0*C1 - 0.5*H0*C1.dag()*C1 - 0.5*C1.dag()*C1*H0 + C1d.dag()*H0*C1d - 0.5*H0*C1d.dag()*C1d - 0.5*C1d.dag()*C1d*H0
	Q2 = C2.dag()*H0*C2 - 0.5*H0*C2.dag()*C2 - 0.5*C2.dag()*C2*H0 + C2d.dag()*H0*C2d - 0.5*H0*C2d.dag()*C2d - 0.5*C2d.dag()*C2d*H0
	Q3 = C3.dag()*H0*C3 - 0.5*H0*C3.dag()*C3 - 0.5*C3.dag()*C3*H0 + C3d.dag()*H0*C3d - 0.5*H0*C3d.dag()*C3d - 0.5*C3d.dag()*C3d*H0

	# Solution to master equation
	rho = steadystate(H0, [C1,C1d,C2,C2d,C3,C3d], tol = 1e-15)
	# order = [1,0,2]
	
	E[j][0] = TR[j]
	E[j][1] = TH
	E[j][2] = TC
	E[j][3] = expect(Q1,rho)
	E[j][4] = expect(Q2,rho)
	E[j][5] = expect(Q3,rho)
	# E[j][6] = expect(h1,rho)
	# E[j][7] = expect(h2,rho)
	# E[j][8] = expect(h3,rho)


file_data_store('2Qubit_J%.1f_g%.3f_k%.3f_TC%.2f_TH%.2f_omegaC%.2f_omegaR%.2f.dat'%(J,g,k1,TC,TH,omegaC,omegaR),E,numtype="real")
# file_data_store('State_J%.1f_g%.3f_k%.3f_TC%.2f_TR%.2f_omegaC%.2f_omegaR%.2f.dat'%(J,g,k1,TC,TR,omegaC,omegaR),rho)
