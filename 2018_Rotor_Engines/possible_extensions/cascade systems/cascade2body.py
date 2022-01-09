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
parser.add_argument("-t", "--t", dest="t_value", type=int, metavar='<int>', default=10, help="T value (default=100)")
parser.add_argument("-nt", "--nt", dest="nt_value", type=int, metavar='<int>', default=1001, help="NT value (default=201)")
parser.add_argument("-th", "--th", dest="th_value", type=float, metavar='<float>', default=2.0, help="TH value (default=5.0)")
parser.add_argument("-tc", "--tc", dest="tc_value", type=float, metavar='<float>', default=1.0, help="TC value (default=1.0)")
parser.add_argument("-tw", "--tw", dest="tw_value", type=float, metavar='<float>', default=8.0, help="TR value (default=2.0)")
parser.add_argument("-omegaC", "--omegaC", dest="omegaC_value", type=float, metavar='<float>', default=0.5, help="W1 value (default=1.0)")
parser.add_argument("-omegaH", "--omegaH", dest="omegaH_value", type=float, metavar='<float>', default=1.0, help="W2 value (default=5.0)")
parser.add_argument("-k1", "--k1", dest="k1_value", type=float, metavar='<float>', default=1, help="k value (default=0.001)")
parser.add_argument("-k2", "--k2", dest="k2_value", type=float, metavar='<float>', default=1, help="k value (default=0.001)")
parser.add_argument("-J", "--J", dest="J_value", type=int, metavar='<int>', default=1, help="J value (default=2.0)")

args = parser.parse_args()
options=Options()
options.nsteps=1e6
options.atol=1e-12
options.rtol=1e-10
options.store_states=True
def meanOccupation(E,Tr):
    if (Tr != 0):
        x = 1.0/(exp(E/Tr)-1)
        return x
    x = 0
    return x

#------PARAMETERS----------

# Timesteps
T=args.t_value
NT=args.nt_value
dt = T*1.0/(NT*1.0)
timesteps=np.linspace(0,T,NT)

# Dimensions for the Modes
dC = 10
dH = 10

# Operators for Harmonic Oscillators
NC = tensor(num(dC),qeye(dH))
NH = tensor(qeye(dC),num(dH))
aC = tensor(destroy(dC),qeye(dH))
aH = tensor(qeye(dC),destroy(dH))
adC = tensor(create(dC),qeye(dH))
adH = tensor(qeye(dC),create(dH))

# Energy of the qubits
omegaC = args.omegaC_value
omegaH = args.omegaH_value

# Temperature and Occupation number
TC = args.tc_value #Qubit 1
TH = args.th_value  #Qubit 2
nc = meanOccupation(omegaC,TC)
nh = meanOccupation(omegaH,TH)
print(nh)
print(meanOccupation(omegaC,TH))
# Coupling strength
k1 = args.k1_value
k2 = args.k2_value

# System Hamiltonian
hC = omegaC*NC
hH = omegaH*NH
hInt = 1j*sqrt(k1*k2)*0.5*(adH*aC - adC*aH)
H = hC + hH + hInt

# Dissipators
C1 = sqrt(nh+1)*(sqrt(k1)*aC + sqrt(k2)*aH)
C2 = sqrt(nh)*(sqrt(k1)*adC + sqrt(k2)*adH)

rho = tensor(fock(dC,0),fock(dH,0))

# Q1 = C1.dag()*H0*C1 - 0.5*H0*C1.dag()*C1 - 0.5*C1.dag()*C1*H0 + C1d.dag()*H0*C1d - 0.5*H0*C1d.dag()*C1d - 0.5*C1d.dag()*C1d*H0
# Q2 = C2.dag()*H0*C2 - 0.5*H0*C2.dag()*C2 - 0.5*C2.dag()*C2*H0 + C2d.dag()*H0*C2d - 0.5*H0*C2d.dag()*C2d - 0.5*C2d.dag()*C2d*H0
# Q3 = C3.dag()*H0*C3 - 0.5*H0*C3.dag()*C3 - 0.5*C3.dag()*C3*H0 + C3d.dag()*H0*C3d - 0.5*H0*C3d.dag()*C3d - 0.5*C3d.dag()*C3d*H0

solution=mesolve(H, rho, timesteps, [C1,C2], [NC,NH],options=options)
probC = np.zeros((NT,dC))
probH = np.zeros((NT,dH))
for j in range(0,NT):
    # rho = 0.5*(solution.states[j] +solution.states[j].dag())
    # rho = rho.full()
    rhoC = solution.states[j].ptrace(0)
    rhoH = solution.states[j].ptrace(1)
    probC[j] = rhoC.diag()
    probH[j] = rhoH.diag()
E=np.vstack((timesteps,solution.expect[0],solution.expect[1]))

file_data_store('Cascade_2HO_kout%.3f_kin%.3f_TC%.2f_TH%.2f_omegaC%.2f_omegaH%.2f.dat'%(k1,k2,TC,TH,omegaC,omegaH),E,numtype="real")
file_data_store('Cascade_2HO_probC_kout%.3f_kin%.3f_TC%.2f_TH%.2f_omegaC%.2f_omegaH%.2f.dat'%(k1,k2,TC,TH,omegaC,omegaH),probC,numtype="real")
file_data_store('Cascade_2HO_probH_kout%.3f_kin%.3f_TC%.2f_TH%.2f_omegaC%.2f_omegaH%.2f.dat'%(k1,k2,TC,TH,omegaC,omegaH),probH,numtype="real")

# file_data_store('State_J%.1f_g%.3f_k%.3f_TC%.2f_TR%.2f_omegaC%.2f_omegaR%.2f.dat'%(J,g,k1,TC,TR,omegaC,omegaR),rho)
