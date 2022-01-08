# cd Documents/CQT/Quantum Thermodynamics/Rotor

# Test Qubit Thermalisation Model using Coarse Graining Master Equation https://arxiv.org/abs/1710.09939
# Model: Two qubits with same frequency interacting, One qubit in a thermal bath

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
parser.add_argument("-t", "--t", dest="t_value", type=int, metavar='<int>', default=1000, help="t value (default=5.0)")
parser.add_argument("-TH", "--TH", dest="TH_value", type=float, metavar='<float>', default=4.0, help="TH value (default=5.0)")
parser.add_argument("-TC", "--TC", dest="TC_value", type=float, metavar='<float>', default=1.0, help="TC value (default=1.0)")
parser.add_argument("-k", "--k", dest="k_value", type=float, metavar='<float>', default=0.01, help="K value (default=1.0)")
parser.add_argument("-g", "--g", dest="g_value", type=float, metavar='<float>', default=0.01, help="K value (default=1.0)")
parser.add_argument("-omegaH", "--omegaH", dest="omegaH_value", type=float, metavar='<float>', default=5.0, help="omega value (default=5.0)")
parser.add_argument("-omegaC", "--omegaC", dest="omegaC_value", type=float, metavar='<float>', default=1.0, help="omega value (default=1.0)")

# parser.add_argument("-delta", "--delta", dest="delta_value", type=float, metavar='<float>', default=0.0, help="DELTA value (default=0.0)")
args = parser.parse_args()

def thermalState(omega,T): #of qubit
    x = np.zeros((2,2),dtype=complex)
    if (T == 0):
    	x[0][0] = 1
    	x = Qobj(x)
    	return x 
    if (T == -1):
	x[1][1] = 1
    	x = Qobj(x)
    	return x 
    x[0][0] = 1.0/(1.0+exp(-omega/T))
    x[1][1] = 1-x[0][0]
    x = Qobj(x)
    return x

#---------------INITIALISATION OF QUBITS---------------
omegaH = args.omegaH_value
omegaC = args.omegaC_value
TC = args.TC_value
TH = args.TH_value
tauH = thermalState(omegaH,TH) 
tauC = thermalState(omegaC,TC)
# Coupling strength
k = args.k_value    
# Interaction Strength 
g = args.g_value
# Timesteps
t=args.t_value
NT=t+1
timesteps=np.linspace(0,t,NT)
options=Options()
options.store_states=True
options.tol= 1e-14

#---------------HAMILTONIANS----------------

# Creating kets and operators
sigmap  = np.zeros((2,2),dtype=complex) 
sigmap[1][0] = 1
sigmap = Qobj(sigmap)
sigmam  = np.zeros((2,2),dtype=complex) 
sigmam[0][1] = 1
sigmam = Qobj(sigmam)
sigmaz  = np.zeros((2,2),dtype=complex) 
sigmaz[0][0] = -1
sigmaz[1][1] = 1
sigmaz = Qobj(sigmaz)
sigmax  = np.zeros((2,2),dtype=complex) 
sigmax[0][1] = 1
sigmax[1][0] = 1
sigmax = Qobj(sigmax)

# Free Hamiltonian
# Hint = g*(tensor(sigmam,sigmap,sigmap)+tensor(sigmap,sigmam,sigmam))
H = 0.5*omegaH*tensor(sigmaz,qeye(2))+0.5*omegaC*tensor(qeye(2),sigmaz)

# Noise
H1 = sqrt(g)*(tensor(sigmax,qeye(2)))
H2 = sqrt(g)*tensor(qeye(2),sigmax)
# H2 = sqrt(g)*tensor(qeye(2),sigmax)
# Coefficients for different dissipators
nX = 1.0/(exp(omegaH*1.0/TH)-1.0)
nY = 1.0/(exp(omegaC*1.0/TC)-1.0)

# # Dissipators 
X = tensor(sigmap,qeye(2))
Y = tensor(qeye(2),sigmap)
CXM = sqrt(k)*sqrt(nX+1)*X.dag() 
CX = sqrt(k)*sqrt(nX)*X
CYM = sqrt(k)*sqrt(nY+1)*Y.dag() 
CY = sqrt(k)*sqrt(nY)*Y

#Initial State
rho = tensor(tauH,tauC)

# Expectation Values
E1 = tensor(sigmaz,qeye(2))
E2 = tensor(qeye(2),sigmaz)
QH = CXM.dag()*H*CXM - 0.5*H*CXM.dag()*CXM - 0.5*CXM.dag()*CXM*H + CX.dag()*H*CX - 0.5*H*CX.dag()*CX - 0.5*CX.dag()*CX*H 
QC = CYM.dag()*H*CYM - 0.5*H*CYM.dag()*CYM - 0.5*CYM.dag()*CYM*H + CY.dag()*H*CY - 0.5*H*CY.dag()*CY - 0.5*CY.dag()*CY*H 

#---------------EXPECTATION VALUES----------------
solution=mesolve(H, rho, timesteps, [CXM,CX,CYM,CY,H1,H2], [E1,E2,QH,QC],options=options)
a1 = solution.expect[0]
a2 = solution.expect[1]
QH = solution.expect[2]
QC = solution.expect[3]
E=np.vstack((timesteps,a1,a2,QH,QC))
file_data_store('Local_t%s_kappa%.3f_omegaH%.3f_omegaC%.3f_TH%.3f_TC%.3f_g%.3f.dat'%(t,k,omegaH,omegaC,TH,TC,g),E,numtype="real")