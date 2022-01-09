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


# NOTE: hbar=1 throughout
#       Indices 0,1,2,...,m_min+m_max correspond to -m_min,-m_min+1,...,0,..., m_max-1,m_max

##############################
## Parse Arguments
#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-nh", "--nh", dest="nh_value", type=float, metavar='<float>', default=1.0, help="TH value (default=5.0)")
parser.add_argument("-nc", "--nc", dest="nc_value", type=float, metavar='<float>', default=0.0, help="TC value (default=1.0)")
parser.add_argument("-g", "--g", dest="g_value", type=float, metavar='<float>', default=3.0, help="G value (default=1.0)")
parser.add_argument("-k", "--k", dest="k_value", type=float, metavar='<float>', default=10.0, help="K value (default=10.0)")
parser.add_argument("-omega", "--omega", dest="omega_value", type=float, metavar='<float>', default=0.01, help="omega value (default=10.0)")

# parser.add_argument("-delta", "--delta", dest="delta_value", type=float, metavar='<float>', default=0.0, help="DELTA value (default=0.0)")
args = parser.parse_args()

def thermalState(N): #of qubit
    x = np.zeros((2,2),dtype=complex)
    x[1][1] = N/(2.0*N+1) 
    x[0][0] = 1-x[1][1]
    x = Qobj(x)
    return x

# For engine to work, TH/TC > wh/wc > 1
# Dimension of rotor

#---------------TWO QUBITS---------------
# Energy of the qubits
n1 = args.nh_value
n3 = args.nc_value

tau1 = thermalState(n1)
tau3 = thermalState(n3)

# Coupling strength
k1 = args.k_value    #effective thermalisation rate for qubit 1 (high)
k3 = args.k_value    #effective thermalisation rate for qubit 3 (high)
# k2 = args.k_value/10.0
# Interaction Strength 
g = sqrt(k1*0.8540927113)

# Driving Frequency
global omega 
omega = 10**(-4-args.omega_value/10.0)
# Timesteps
T=2*pi/omega
NT=101
timesteps=np.linspace(0,T,NT)


# Frequency of Mode
omega0 = 0
options=Options()
options.nsteps=1e8
options.atol=1e-12
options.method = 'bdf'
options.store_states=True
options.max_step = (T*1.0)/(NT*1.0)

#---------------HAMILTONIANS----------------

# Creating Sigma+ and Sigma-
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

# Interaction Hamiltonian
Hint1=g*tensor(sigmam,sigmap)
Hint2=g*tensor(sigmap,sigmam)
# Total Hamiltonian

def Hint1_coeff(t,args):
    global omega
    return cos(omega*t) + 1j *sin(omega*t)

def Hint2_coeff(t,args):
    global omega
    return cos(omega*t) - 1j *sin(omega*t)

H = [[Hint1,Hint1_coeff],[Hint2,Hint2_coeff]]

# Lindblad Dissipators
C_1=sqrt(k1*(n1+1))*tensor(sigmam,qeye(2))
C_1_d=sqrt(k1*n1)*tensor(sigmap,qeye(2))
C_3=sqrt(k3*(n3+1))*tensor(qeye(2),sigmam)
C_3_d=sqrt(k3*n3)*tensor(qeye(2),sigmap)


rho = tensor(tau1,tau3)


#---------------EXPECTATION VALUES----------------
a1 = tensor(sigmam,qeye(2))
ad1 = tensor(sigmap,qeye(2))
a3 = tensor(qeye(2),sigmam)
ad3 = tensor(qeye(2),sigmap)
Nh = tensor(sigmap*sigmam,qeye(2))
Nc = tensor(qeye(2),sigmap*sigmam)
W = 1j*g*omega*tensor(sigmam,sigmap)
K1 = tensor(sigmap*sigmam,qeye(2))
K3 = tensor(qeye(2),sigmap*sigmam)
EQh = k1*(n1+1)*(ad1*K1*a1 - 0.5*ad1*a1*K1 - 0.5*K1*ad1*a1) + k1*n1*(a1*K1*ad1 - 0.5*a1*ad1*K1 - 0.5*K1*a1*ad1)
EQc = k3*(n3+1)*(ad3*K3*a3 - 0.5*ad3*a3*K3 - 0.5*K3*ad3*a3) + k3*n3*(a3*K3*ad3 - 0.5*a3*ad3*K3 - 0.5*K3*a3*ad3)
solution=mesolve(H, rho, timesteps, [C_1,C_1_d,C_3,C_3_d],[W,EQh,EQc,Nh,Nc],options=options)

sinPhi = np.sin(timesteps*omega)
cosPhi = np.cos(timesteps*omega)
W = np.multiply(solution.expect[0],cosPhi + 1j*sinPhi)
W = 2.0*np.real(W)
Qh = solution.expect[1]
Qc = solution.expect[2]
Nh = solution.expect[3]
Nc = solution.expect[4]

dataThermo=np.vstack((timesteps,W,Qh,Qc,Nh,Nc)) 
file_data_store('Rotor_NT%s_kappa%s_nH%s_nC%s_omega%s.dat'%(NT,k1,n1,n3,0.4+args.omega_value/100.0),dataThermo,numtype="real")