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
parser.add_argument("-t", "--t", dest="t_value", type=float, metavar='<float>', default=10000, help="t value (default=5.0)")
parser.add_argument("-th", "--th", dest="th_value", type=float, metavar='<float>', default=5.0, help="TH value (default=5.0)")
parser.add_argument("-tc", "--tc", dest="tc_value", type=float, metavar='<float>', default=1.0, help="TC value (default=1.0)")
parser.add_argument("-tr", "--tr", dest="tr_value", type=float, metavar='<float>', default=2.0, help="TR value (default=2.0)")
parser.add_argument("-omegaC", "--omegaC", dest="omegaC_value", type=float, metavar='<float>', default=1.0, help="W1 value (default=1.0)")
parser.add_argument("-omegaR", "--omegaR", dest="omegaR_value", type=float, metavar='<float>', default=5.0, help="W2 value (default=5.0)")
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
# Timesteps
t=args.t_value
NT=1001
timesteps=np.linspace(0,t,NT)

# Spin operators
d = args.J_value+1
J = args.J_value*0.5
sigmap = jmat(J,'+')
sigmam = jmat(J,'-')
sigmaz = jmat(J,'z')/J #normalise sigmaz to range +/- 1
sigmax = jmat(J,'x')
sigmay = jmat(J,'y')
tmp =  np.diag(sigmaz.full().real)
# Energy of the qubits
omegaC = args.omegaC_value
omegaR = args.omegaR_value
omegaH = omegaR-omegaC

TC = args.tc_value #Qubit 1
TR = args.tr_value  #Qubit 2
TH = args.th_value 
# Coupling strength
k1 = args.k_value
k2 = args.k_value
k3 = args.k_value

# System Hamiltonian
h1 = 0.5*omegaC*tensor(sigmaz,qeye(d),qeye(d)) 
h2 = 0.5*omegaR*tensor(qeye(d),sigmaz,qeye(d))
h3 = 0.5*omegaH*tensor(qeye(d),qeye(d),sigmaz) 
H0 = h1 + h2 + h3


# Interaction Hamiltonian
g = args.g_value
Hint = g*(tensor(sigmap,sigmam,sigmap)+tensor(sigmam,sigmap,sigmam))

# Total Hamiltonian
H = H0 + Hint

# Temperature of baths


# Mean Occupation of bath
N1 = meanOccupation(omegaC/(d-1.0),TC)
N2 = meanOccupation(omegaR/(d-1.0),TR)
N3 = meanOccupation(omegaH/(d-1.0),TH)
# tauH = thermalState(omegaH,TH) #Qubit 1 starts with ground state
# tauC = thermalState(omegaC,TC) #Qubit 2 is in an excited state 
# tauR = thermalState(omegaR,TR)

# Thermal state of the qubits
tau1 = np.exp(-tmp*0.5*omegaC/(TC*1.0))
tau2 = np.exp(-tmp*0.5*omegaR/(TR*1.0))
tau3 = np.exp(-tmp*0.5*omegaH/(TH*1.0))
tau1 = np.diag(tau1/np.sum(tau1))
tau2 = np.diag(tau2/np.sum(tau2))
tau3 = np.diag(tau3/np.sum(tau3))

tau1 = Qobj(tau1)
tau2 = Qobj(tau2)
tau3 = Qobj(tau3)
# Lindblad Dissipators

C1=sqrt(k1*(N1+1))*tensor(sigmam,qeye(d),qeye(d))
C1d=sqrt(k1*(N1))*tensor(sigmap,qeye(d),qeye(d))
C2=sqrt(k2*(N2+1))*tensor(qeye(d),sigmam,qeye(d))
C2d=sqrt(k2*(N2))*tensor(qeye(d),sigmap,qeye(d))
C3=sqrt(k3*(N3+1))*tensor(qeye(d),qeye(d),sigmam)
C3d=sqrt(k3*(N3))*tensor(qeye(d),qeye(d),sigmap)

Q1 = C1.dag()*h1*C1 - 0.5*h1*C1.dag()*C1 - 0.5*C1.dag()*C1*h1 + C1d.dag()*h1*C1d - 0.5*h1*C1d.dag()*C1d - 0.5*C1d.dag()*C1d*h1
Q2 = C2.dag()*h2*C2 - 0.5*h2*C2.dag()*C2 - 0.5*C2.dag()*C2*h2 + C2d.dag()*h2*C2d - 0.5*h2*C2d.dag()*C2d - 0.5*C2d.dag()*C2d*h2
Q3 = C3.dag()*h3*C3 - 0.5*h3*C3.dag()*C3 - 0.5*C3.dag()*C3*h3 + C3d.dag()*h3*C3d - 0.5*h3*C3d.dag()*C3d - 0.5*C3d.dag()*C3d*h3

rho = tensor(tau1,tau2,tau3)
# Frequency of Mode
options=Options()
options.store_states=True
options.tol= 1e-14
# Solution to master equation
solution=mesolve(H, rho, timesteps, [C1,C1d,C2,C2d,C3,C3d],[Q1,Q2,Q3,h1,h2,h3],options=options)

E1 = np.zeros((1,NT))
E2 = np.zeros((1,NT))
E3 = np.zeros((1,NT))
S = np.zeros((1,NT))
for j in range(0,NT):
    E1[0][j] = expect(solution.states[j].ptrace(0),solution.states[j].ptrace(0))
    E2[0][j] = expect(solution.states[j].ptrace(1),solution.states[j].ptrace(1))
    E3[0][j] = expect(solution.states[j].ptrace(2),solution.states[j].ptrace(2))
    S[0][j] = entropy_vn(solution.states[j])
QH = solution.expect[0]
QC = solution.expect[1]
QR = solution.expect[2]
H1 = solution.expect[3]
H2 = solution.expect[4]
H3 = solution.expect[5]
E=np.vstack((timesteps,H1,H2,H3,QH,QC,QR,S,QH/TH,QC/TC,QR/TR,E1,E2,E3))
file_data_store('Data_J%.1f_g%.3f_k%.3f_TC%.2f_TR%.2f_omegaC%.2f_omegaR%.2f.dat'%(J,g,k1,TC,TR,omegaC,omegaR),E,numtype="real")
# file_data_store('State_J%.1f_g%.3f_k%.3f_TC%.2f_TR%.2f_omegaC%.2f_omegaR%.2f.dat'%(J,g,k1,TC,TR,omegaC,omegaR),rho)
