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
parser.add_argument("-t", "--t", dest="t_value", type=int, metavar='<int>', default=50, help="T value (default=100)")
parser.add_argument("-nt", "--nt", dest="nt_value", type=int, metavar='<int>', default=1001, help="NT value (default=201)")
parser.add_argument("-th", "--th", dest="th_value", type=float, metavar='<float>', default=1, help="TH value (default=5.0)")
parser.add_argument("-tc", "--tc", dest="tc_value", type=float, metavar='<float>', default=0.5, help="TC value (default=1.0)")
parser.add_argument("-tw", "--tw", dest="tw_value", type=float, metavar='<float>', default=4, help="TR value (default=2.0)")
parser.add_argument("-omegaC", "--omegaC", dest="omegaC_value", type=float, metavar='<float>', default=1.0, help="W1 value (default=1.0)")
parser.add_argument("-omegaW", "--omegaW", dest="omegaW_value", type=float, metavar='<float>', default=4.0, help="W2 value (default=5.0)")
parser.add_argument("-g", "--g", dest="g_value", type=float, metavar='<float>', default=0.1, help="g value (default=0.001)")
parser.add_argument("-k1", "--k1", dest="k1_value", type=float, metavar='<float>', default=0.1, help="k1 value (default=0.001)")
parser.add_argument("-k2", "--k2", dest="k2_value", type=float, metavar='<float>', default=0.1, help="k2 value (default=0.001)")
parser.add_argument("-kc", "--kc", dest="kc_value", type=float, metavar='<float>', default=0.1, help="kc value (default=0.001)")
parser.add_argument("-kw", "--kw", dest="kw_value", type=float, metavar='<float>', default=0.1, help="kw value (default=0.001)")

args = parser.parse_args()

def meanOccupation(E,Tr):
    if (Tr != 0):
        x = 1.0/(exp(E/Tr)-1)
        return x
    x = 0
    return x
options=Options()
options.nsteps=1e6
options.atol=1e-12
options.rtol=1e-10
options.store_states=True
#------PARAMETERS----------
# Timesteps
T=args.t_value
NT=args.nt_value
dt = T*1.0/(NT*1.0)
timesteps=np.linspace(0,T,NT)
# Dimensions for the Modes
dC = 5
dH = 5
dW = 12
# Operators for Harmonic Oscillators
NC = tensor(num(dC),qeye(dH),qeye(dW))
NH = tensor(qeye(dC),num(dH),qeye(dW))
NW = tensor(qeye(dC),qeye(dH),num(dW))
aC = tensor(destroy(dC),qeye(dH),qeye(dW))
aH = tensor(qeye(dC),destroy(dH),qeye(dW))

aW = tensor(qeye(dC),qeye(dH),destroy(dW))
adC = tensor(create(dC),qeye(dH),qeye(dW))
adH = tensor(qeye(dC),create(dH),qeye(dW))
adW = tensor(qeye(dC),qeye(dH),create(dW))

# Energy of the qubits
omegaC = args.omegaC_value
omegaW = args.omegaW_value
omegaH = omegaW+omegaC

TC = args.tc_value #Qubit 1
TH = args.th_value  #Qubit 2
TW = args.tw_value  #Qubit 2
nc = meanOccupation(omegaC,TC)
nw = meanOccupation(omegaW,TW)
nh = meanOccupation(omegaH,TH)
print(nc)
print(nw)
print(nh)
# Coupling strength
k1 = args.k1_value
k2 = args.k2_value
kc = args.kc_value
kw = args.kw_value


# System Hamiltonian
hC = omegaC*NC
hH = omegaH*NH 
hW = omegaW*NW


# Interaction Hamiltonian
g = args.g_value
hInt = -g*(NC)*(aW+adW)

# Full Hamiltonian 
H = hInt + hC + hH + hW

# Dissipators
C1 = sqrt(nh+1)*(sqrt(k1)*aC + sqrt(k2)*aH)
C2 = sqrt(nh)*(sqrt(k1)*adC + sqrt(k2)*adH)
C3 = sqrt(kc*(nc+1))*aC
C4 = sqrt(kc*nc)*adC
C5 = sqrt(kw*(nw+1))*aW
C6 = sqrt(kw*nw)*adW
H = H + 1j*sqrt(k1*k2)*0.5*(adH*aC - adC*aH)
rho = tensor(thermal_dm(dC,nc),thermal_dm(dH,nh),thermal_dm(dW,nw))

QC = C3.dag()*hC*C3 - 0.5*hC*C3.dag()*C3 - 0.5*C3.dag()*C3*hC + C4.dag()*hC*C4 - 0.5*hC*C4.dag()*C4 - 0.5*C4.dag()*C4*hC

solution=mesolve(H, rho, timesteps, [C1,C2,C3,C4,C5,C6], [NC,NH,NW,QC],options=options)
probC = np.zeros((NT,dC))
probH = np.zeros((NT,dH))
probW = np.zeros((NT,dW))
for j in range(0,NT):
    rhoC = solution.states[j].ptrace(0)
    rhoH = solution.states[j].ptrace(1)
    rhoW = solution.states[j].ptrace(2)
    probC[j] = rhoC.diag()
    probH[j] = rhoH.diag()
    probW[j] = rhoW.diag()

E=np.vstack((timesteps,solution.expect[0],solution.expect[1],solution.expect[2],solution.expect[3]))

file_data_store('Cascade_3HO_kout%.3f_kin%.3f_TC%.2f_TH%.2f_TW%.2f_omegaC%.2f_omegaH%.2f.dat'%(k1,k2,TC,TH,TW,omegaC,omegaH),E,numtype="real")
file_data_store('Cascade_3HO_probC_kout%.3f_kin%.3f_TC%.2f_TH%.2f_TW%.2f_omegaC%.2f_omegaH%.2f.dat'%(k1,k2,TC,TH,TW,omegaC,omegaH),probC,numtype="real")
file_data_store('Cascade_3HO_probH_kout%.3f_kin%.3f_TC%.2f_TH%.2f_TW%.2f_omegaC%.2f_omegaH%.2f.dat'%(k1,k2,TC,TH,TW,omegaC,omegaH),probH,numtype="real")