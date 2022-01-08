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
parser.add_argument("-TH", "--TH", dest="TH_value", type=float, metavar='<float>', default=2.0, help="TH value (default=5.0)")
parser.add_argument("-TC", "--TC", dest="TC_value", type=float, metavar='<float>', default=1.0, help="TC value (default=1.0)")
parser.add_argument("-TR", "--TR", dest="TR_value", type=float, metavar='<float>', default=5.0, help="TR value (default=2.0)")
parser.add_argument("-g", "--g", dest="g_value", type=float, metavar='<float>', default=0.5, help="G value (default=1.0)")
parser.add_argument("-k", "--k", dest="k_value", type=float, metavar='<float>', default=0.01, help="K value (default=1.0)")
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
omegaR = omegaH - omegaC


TC = args.TC_value
TH = args.TH_value
TR = args.TR_value
tauH = thermalState(omegaH,TH) #Qubit 1 starts with ground state
tauC = thermalState(omegaC,TC) #Qubit 2 is in an excited state 
tauR = thermalState(omegaR,TR)

# Coupling strength
k = args.k_value#10**(-(m+1.0)/10.0)  

# Interaction Strength 
g = args.g_value#10**(-(mm+1.0)/10.0)


# Timesteps
t=args.t_value
NT=t+1
timesteps=np.linspace(0,t,NT)


# Frequency of Mode
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
ket0 = basis(2,0)
ket1 = basis(2,1)
ket000 = tensor(ket0,ket0,ket0)
ket001 = tensor(ket0,ket0,ket1)
ket010 = tensor(ket0,ket1,ket0)
ket011 = tensor(ket0,ket1,ket1)
ket100 = tensor(ket1,ket0,ket0)
ket101 = tensor(ket1,ket0,ket1)
ket110 = tensor(ket1,ket1,ket0)
ket111 = tensor(ket1,ket1,ket1)
ketplus = (ket100+ket011)/sqrt(2)
ketminus = (ket100-ket011)/sqrt(2)


# Hamiltonians
Hint = g*(tensor(sigmax,sigmax,sigmax))
H0 = omegaH*tensor(sigmap*sigmam,qeye(2),qeye(2))+omegaC*tensor(qeye(2),sigmap*sigmam,qeye(2))+omegaR*tensor(qeye(2),qeye(2),sigmap*sigmam)
H = H0 + Hint
tmpEig = H.eigenstates()
ketE0m = tmpEig[1][3]
ketE1m = tmpEig[1][2]
ketE2m = tmpEig[1][1]
ketE3m = tmpEig[1][0]
ketE0p = tmpEig[1][4]
ketE1p = tmpEig[1][5]
ketE2p = tmpEig[1][6]
ketE3p = tmpEig[1][7]
E0m = tmpEig[0][3]
E1m = tmpEig[0][2]
E2m = tmpEig[0][1]
E3m = tmpEig[0][0]
E0p = tmpEig[0][4]
E1p = tmpEig[0][5]
E2p = tmpEig[0][6]
E3p = tmpEig[0][7]

# Lindblad Dissipators Qubit 1
sigmax1 = tensor(sigmax,qeye(2),qeye(2))
X1 = ketE3p*ketE0m.dag()* (ketE3p.dag()*sigmax1*ketE0m) + ketE0p*ketE3m.dag()* (ketE0p.dag()*sigmax1*ketE3m)
omegaX1 = E3p-E0m
X2 = ketE3p*ketE0p.dag()* (ketE3p.dag()*sigmax1*ketE0p) + ketE0m*ketE3m.dag()* (ketE0m.dag()*sigmax1*ketE3m)
omegaX2 = E3p-E0p
X3 = ketE2p*ketE1m.dag()* (ketE2p.dag()*sigmax1*ketE1m) + ketE1p*ketE2m.dag()* (ketE1p.dag()*sigmax1*ketE2m)
omegaX3 = E2p-E1m
X4 = ketE2p*ketE1p.dag()* (ketE2p.dag()*sigmax1*ketE1p) + ketE1m*ketE2m.dag()* (ketE1m.dag()*sigmax1*ketE2m)
omegaX4 = E2p-E1p

# Lindblad Dissipators Qubit 2
sigmax2 = tensor(qeye(2),sigmax,qeye(2))
if (omegaC<omegaR):
    Y1 = ketE3p*ketE2m.dag()* (ketE3p.dag()*sigmax2*ketE2m) + ketE2p*ketE3m.dag()* (ketE2p.dag()*sigmax2*ketE3m)
    omegaY1 = E3p-E2m
    Y2 = ketE3p*ketE2p.dag()* (ketE3p.dag()*sigmax2*ketE2p) + ketE2m*ketE3m.dag()* (ketE2m.dag()*sigmax2*ketE3m)
    omegaY2 = E3p-E2p
    Y3 = ketE1p*ketE0m.dag()* (ketE1p.dag()*sigmax2*ketE0m) + ketE0p*ketE1m.dag()* (ketE0p.dag()*sigmax2*ketE1m)
    omegaY3 = E1p-E0m
    Y4 = ketE1p*ketE0p.dag()* (ketE1p.dag()*sigmax2*ketE0p) + ketE0m*ketE1m.dag()* (ketE0m.dag()*sigmax2*ketE1m)
    omegaY4 = E1p-E0p
else:
    Y1 = ketE3p*ketE1m.dag()* (ketE3p.dag()*sigmax2*ketE1m) + ketE1p*ketE3m.dag()* (ketE1p.dag()*sigmax2*ketE3m)
    omegaY1 = E3p-E1m
    Y2 = ketE3p*ketE1p.dag()* (ketE3p.dag()*sigmax2*ketE1p) + ketE1m*ketE3m.dag()* (ketE1m.dag()*sigmax2*ketE3m)
    omegaY2 = E3p-E1p
    Y3 = ketE2p*ketE0m.dag()* (ketE2p.dag()*sigmax2*ketE0m) + ketE0p*ketE2m.dag()* (ketE0p.dag()*sigmax2*ketE2m)
    omegaY3 = E2p-E0m
    Y4 = ketE2p*ketE0p.dag()* (ketE2p.dag()*sigmax2*ketE0p) + ketE0m*ketE2m.dag()* (ketE0m.dag()*sigmax2*ketE2m)
    omegaY4 = E2p-E0p

# Lindblad Dissipators Qubit 3
sigmax3 = tensor(qeye(2),qeye(2),sigmax)
if (omegaC<omegaR):
    Z1 = ketE3p*ketE1m.dag()* (ketE3p.dag()*sigmax3*ketE1m) + ketE1p*ketE3m.dag()* (ketE1p.dag()*sigmax3*ketE3m)
    omegaZ1 = E3p-E1m
    Z2 = ketE3p*ketE1p.dag()* (ketE3p.dag()*sigmax3*ketE1p) + ketE1m*ketE3m.dag()* (ketE1m.dag()*sigmax3*ketE3m)
    omegaZ2 = E3p-E1p
    Z3 = ketE2p*ketE0m.dag()* (ketE2p.dag()*sigmax3*ketE0m) + ketE0p*ketE2m.dag()* (ketE0p.dag()*sigmax3*ketE2m)
    omegaZ3 = E2p-E0m
    Z4 = ketE2p*ketE0p.dag()* (ketE2p.dag()*sigmax3*ketE0p) + ketE0m*ketE2m.dag()* (ketE0m.dag()*sigmax3*ketE2m)
    omegaZ4 = E2p-E0p
else:
    Z1 = ketE3p*ketE2m.dag()* (ketE3p.dag()*sigmax3*ketE2m) + ketE2p*ketE3m.dag()* (ketE2p.dag()*sigmax3*ketE3m)
    omegaZ1 = E3p-E2m
    Z2 = ketE3p*ketE2p.dag()* (ketE3p.dag()*sigmax3*ketE2p) + ketE2m*ketE3m.dag()* (ketE2m.dag()*sigmax3*ketE3m)
    omegaZ2 = E3p-E2p
    Z3 = ketE1p*ketE0m.dag()* (ketE1p.dag()*sigmax3*ketE0m) + ketE0p*ketE1m.dag()* (ketE0p.dag()*sigmax3*ketE1m)
    omegaZ3 = E1p-E0m
    Z4 = ketE1p*ketE0p.dag()* (ketE1p.dag()*sigmax3*ketE0p) + ketE0m*ketE1m.dag()* (ketE0m.dag()*sigmax3*ketE1m)
    omegaZ4 = E1p-E0p

# Coefficients for different dissipators
nX1 = 1.0/(exp(abs(omegaX1)*1.0/TH)-1.0)
nX2 = 1.0/(exp(abs(omegaX2)*1.0/TH)-1.0)
nX3 = 1.0/(exp(abs(omegaX3)*1.0/TH)-1.0)
nX4 = 1.0/(exp(abs(omegaX4)*1.0/TH)-1.0)
nY1 = 1.0/(exp(abs(omegaY1)*1.0/TC)-1.0)
nY2 = 1.0/(exp(abs(omegaY2)*1.0/TC)-1.0)
nY3 = 1.0/(exp(abs(omegaY3)*1.0/TC)-1.0)
nY4 = 1.0/(exp(abs(omegaY4)*1.0/TC)-1.0)
nZ1 = 1.0/(exp(abs(omegaZ1)*1.0/TR)-1.0)
nZ2 = 1.0/(exp(abs(omegaZ2)*1.0/TR)-1.0)
nZ3 = 1.0/(exp(abs(omegaZ3)*1.0/TR)-1.0)
nZ4 = 1.0/(exp(abs(omegaZ4)*1.0/TR)-1.0)
# # Dissipators 
# if (omega>g):
#     x1 = 0.5*(tensor(qeye(2),sigmap)+tensor(sigmap,sigmaz))
#     x2 = 0.5*(tensor(qeye(2),sigmap)-tensor(sigmap,sigmaz))
# else:
#     x1 = 0.5*(tensor(qeye(2),sigmam)+tensor(sigmam,sigmaz))
#     x2 = 0.5*(tensor(qeye(2),sigmap)-tensor(sigmap,sigmaz)) 


# Lindblad Dissipators
CXM4= sqrt(k)*sqrt(nX4+1) * X3.dag()
CXM3= sqrt(k)*sqrt(nX3+1) * X3.dag()
CXM2= sqrt(k)*sqrt(nX2+1) * X2.dag()
CXM1= sqrt(k)*sqrt(nX1+1) * X1.dag()
CX1 = sqrt(k)*sqrt(nX1) * X1
CX2 = sqrt(k)*sqrt(nX2) * X2
CX3 = sqrt(k)*sqrt(nX3) * X3
CX4 = sqrt(k)*sqrt(nX4) * X4

CYM4= sqrt(k)*sqrt(nY4+1) * Y4.dag()
CYM3= sqrt(k)*sqrt(nY3+1) * Y3.dag()
CYM2= sqrt(k)*sqrt(nY2+1) * Y2.dag()
CYM1= sqrt(k)*sqrt(nY1+1) * Y1.dag()
CY1 = sqrt(k)*sqrt(nY1) * Y1
CY2 = sqrt(k)*sqrt(nY2) * Y2
CY3 = sqrt(k)*sqrt(nY3) * Y3
CY4 = sqrt(k)*sqrt(nY4) * Y4

CZM4= sqrt(k)*sqrt(nZ4+1) * Z4.dag()
CZM3= sqrt(k)*sqrt(nZ3+1) * Z3.dag()
CZM2= sqrt(k)*sqrt(nZ2+1) * Z2.dag()
CZM1= sqrt(k)*sqrt(nZ1+1) * Z1.dag()
CZ1 = sqrt(k)*sqrt(nZ1) * Z1
CZ2 = sqrt(k)*sqrt(nZ2) * Z2
CZ3 = sqrt(k)*sqrt(nZ3) * Z3
CZ4 = sqrt(k)*sqrt(nZ4) * Z4

#Initial State
rho = tensor(tauH,tauC,tauR)

# Expectation Values
E1 = tensor(sigmaz,qeye(2),qeye(2))
E2 = tensor(qeye(2),sigmaz,qeye(2))
E3 = tensor(qeye(2),qeye(2),sigmaz)

# Probabilities
P1 = ket000*ket000.dag()
P2 = ket001*ket001.dag()
P3 = ket010*ket010.dag()
P4 = ket011*ket011.dag()
P5 = ket100*ket100.dag()
P6 = ket101*ket101.dag()
P7 = ket110*ket110.dag()
P8 = ket111*ket111.dag()
P9 = ketplus*ketplus.dag()
P10 = ketminus*ketminus.dag()

QH = CXM4.dag()*H*CXM4 - 0.5*H*CXM4.dag()*CXM4 - 0.5*CXM4.dag()*CXM4*H + CXM3.dag()*H*CXM3 - 0.5*H*CXM3.dag()*CXM3 - 0.5*CXM3.dag()*CXM3*H + CXM2.dag()*H*CXM2 - 0.5*H*CXM2.dag()*CXM2 - 0.5*CXM2.dag()*CXM2*H + CXM1.dag()*H*CXM1 - 0.5*H*CXM1.dag()*CXM1 - 0.5*CXM1.dag()*CXM1*H + CX4.dag()*H*CX3 - 0.5*H*CX4.dag()*CX3 - 0.5*CX4.dag()*CX3*H + CX3.dag()*H*CX3 - 0.5*H*CX3.dag()*CX3 - 0.5*CX3.dag()*CX3*H + CX2.dag()*H*CX2 - 0.5*H*CX2.dag()*CX2 - 0.5*CX2.dag()*CX2*H + CX1.dag()*H*CX1 - 0.5*H*CX1.dag()*CX1 - 0.5*CX1.dag()*CX1*H
QC = CYM4.dag()*H*CYM4 - 0.5*H*CYM4.dag()*CYM4 - 0.5*CYM4.dag()*CYM4*H + CYM3.dag()*H*CYM3 - 0.5*H*CYM3.dag()*CYM3 - 0.5*CYM3.dag()*CYM3*H + CYM2.dag()*H*CYM2 - 0.5*H*CYM2.dag()*CYM2 - 0.5*CYM2.dag()*CYM2*H + CYM1.dag()*H*CYM1 - 0.5*H*CYM1.dag()*CYM1 - 0.5*CYM1.dag()*CYM1*H + CY4.dag()*H*CY3 - 0.5*H*CY4.dag()*CY3 - 0.5*CY4.dag()*CY3*H + CY3.dag()*H*CY3 - 0.5*H*CY3.dag()*CY3 - 0.5*CY3.dag()*CY3*H + CY2.dag()*H*CY2 - 0.5*H*CY2.dag()*CY2 - 0.5*CY2.dag()*CY2*H + CY1.dag()*H*CY1 - 0.5*H*CY1.dag()*CY1 - 0.5*CY1.dag()*CY1*H
QR = CZM4.dag()*H*CYM4 - 0.5*H*CZM4.dag()*CZM4 - 0.5*CZM4.dag()*CZM4*H + CZM3.dag()*H*CZM3 - 0.5*H*CZM3.dag()*CZM3 - 0.5*CZM3.dag()*CZM3*H + CZM2.dag()*H*CZM2 - 0.5*H*CZM2.dag()*CZM2 - 0.5*CZM2.dag()*CZM2*H + CZM1.dag()*H*CZM1 - 0.5*H*CZM1.dag()*CZM1 - 0.5*CZM1.dag()*CZM1*H + CZ4.dag()*H*CZ3 - 0.5*H*CZ4.dag()*CZ3 - 0.5*CZ4.dag()*CZ3*H + CZ3.dag()*H*CZ3 - 0.5*H*CZ3.dag()*CZ3 - 0.5*CZ3.dag()*CZ3*H + CZ2.dag()*H*CZ2 - 0.5*H*CZ2.dag()*CZ2 - 0.5*CZ2.dag()*CZ2*H + CZ1.dag()*H*CZ1 - 0.5*H*CZ1.dag()*CZ1 - 0.5*CZ1.dag()*CZ1*H

#---------------EXPECTATION VALUES----------------
solution=mesolve(H, rho, timesteps, [CXM4,CXM3,CXM2,CXM1,CX1,CX2,CX3,CX4,CYM4,CYM3,CYM2,CYM1,CY1,CY2,CY3,CY4,CZM4,CZM3,CZM2,CZM1,CZ1,CZ2,CZ3], [E1,E2,E3,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,QH,QC,QR],options=options) #
S = np.zeros((1,NT))
# prob1 = np.zeros((NT,2))
# prob2 = np.zeros((NT,2))

for j in range(0,NT):
    S[0][j] = entropy_vn(solution.states[j])

a1 = solution.expect[0]
a2 = solution.expect[1]
a3 = solution.expect[2]
p1 = solution.expect[3]
p2 = solution.expect[4]
p3 = solution.expect[5]
p4 = solution.expect[6]
p5 = solution.expect[7]
p6 = solution.expect[8]
p7 = solution.expect[9]
p8 = solution.expect[10]
p9 = solution.expect[11]
p10 = solution.expect[12]
QH = solution.expect[13]
QC = solution.expect[14]
QR = solution.expect[15]
E=np.vstack((timesteps,a1,a2,a3,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,QH,QC,QR,S,QH/TH,QC/TC,QR/TR))
file_data_store('Global_t%s_kappa%.3f_omegaH%.3f_omegaC%.3f_TH%.3f_TC%.3f_TR%.3f_g%.3f.dat'%(t,k,omegaH,omegaC,TH,TC,TR,g),E,numtype="real")
