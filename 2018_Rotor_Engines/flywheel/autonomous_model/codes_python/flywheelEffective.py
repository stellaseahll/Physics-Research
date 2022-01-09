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
from extractWork import extractWork 
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
parser.add_argument("-T", "--T", dest="t_value", type=float, metavar='<float>', default=200.0, help="T value (default=10.0)")
parser.add_argument("-NT", "--NT", dest="nt_value", type=int, metavar='<int>', default=1001, help="NT value (default=51.0)")
parser.add_argument("-nh", "--nh", dest="nh_value", type=float, metavar='<float>', default=1.0, help="TH value (default=5.0)")
parser.add_argument("-nc", "--nc", dest="nc_value", type=float, metavar='<float>', default=0.0, help="TC value (default=1.0)")
parser.add_argument("-mmin", "--mmin", dest="mmin_value", type=int, metavar='<int>', default=-20, help="MMIN value (default=-20)")
parser.add_argument("-mmax", "--mmax", dest="mmax_value", type=int, metavar='<int>', default=100, help="MMAX value (default=100)")
parser.add_argument("-g", "--g", dest="g_value", type=float, metavar='<float>', default=1.0, help="G value (default=1.0)")
parser.add_argument("-k", "--k", dest="k_value", type=float, metavar='<float>', default=1.0, help="K value (default=10.0)")
parser.add_argument("-I", "--I", dest="I_value", type=float, metavar='<float>', default=10.0, help="I value (default=10.0)")

# parser.add_argument("-delta", "--delta", dest="delta_value", type=float, metavar='<float>', default=0.0, help="DELTA value (default=0.0)")
args = parser.parse_args()



#---------------TWO QUBITS---------------
# Energy of the qubits
n1 = args.nh_value
n3 = args.nc_value


#---------------ROTOR----------------
# Maximum angular momentum quantum number of the rotor
# 40 to 80 or 20 to 250
m_max=args.mmax_value
m_min=args.mmin_value

# Dimension of rotor Hilbert space
d2=(m_max-m_min)+1

# Momentum of inertia (If too small, dispersion occurs too quickly)
#10 or 100
MI=args.I_value

# Kappa parameter of von Mises distribution (The larger it is the more peaked the distribution will be) 
kappa=1

# Mean of von Mises distribution (Initial expectation value of the angle) 
mu=0.0

# von Mises function
def vonMises(x,mu,kappa):
    return exp(kappa*cos(x-mu))

# Normalization constant
norm_vonMises=integrate.quad(lambda x: (vonMises(x,mu,kappa))**2,0,2*pi)[0] 

# von Mises wavefunction
def int_psi(m):
    I1=integrate.quad(lambda x: sqrt(1/norm_vonMises)*(1/(sqrt(2*pi)))*cos(-(m_min+m)*x)*vonMises(x,mu,kappa),0,2*pi,limit=m_max*100)
    I2=integrate.quad(lambda x: sqrt(1/norm_vonMises)*(1/(sqrt(2*pi)))*sin(-(m_min+m)*x)*vonMises(x,mu,kappa),0,2*pi,limit=m_max*100)    
    return I1[0]+1j*I2[0]  

#Initial state of rotor
rot_ket=Qobj([[int_psi(m)] for m in range(d2)])
rot_ket=rot_ket.unit()
rot_bra=rot_ket.dag()

Rot=rot_ket*rot_bra
Rot=Rot/(Rot.tr())

# Initial state of rotor
rho = Rot


# Coupling strength
k1 = args.k_value    #thermalisation rate for qubit 1 (high)
k3 = args.k_value    #thermalisation rate for qubit 3 (high)
# k2 = args.k_value/10.0
# Interaction Strength 
g = 0.924


# Thermalisation rate
kh = (2.0*g)**2/k1*(n1+1)*n3/(2.0*n1+1+2.0*n3+1)/(2.0*n1+1)/(2.0*n3+1)
kc = (2.0*g)**2/k3*(n3+1)*n1/(2.0*n1+1+2.0*n3+1)/(2.0*n1+1)/(2.0*n3+1)
# n = ((wh*1.0/TH)+(wc*1.0/TC))/(wh-wc)
print('kh = %.8f s\n'%(kh))  
print('kc = %.8f s\n'%(kc))   
#-------OPERATORS DEFINITION---------
# Integral for the matrix elements of any function of the angle operator
# def int_function(m,n):
#     I1=integrate.quad(lambda x: (1/(2*pi))*cos(((-m_max+n)-(-m_max+m))*x)*function(x),0,2*pi,limit=m_max*100)
#     I2=integrate.quad(lambda x: (1/(2*pi))*sin(((-m_max+n)-(-m_max+m))*x)*function(x),0,2*pi,limit=m_max*100)
#     return I1[0]+1j*I2[0]

# Initialize matrix with zeros
lm=np.zeros((d2,d2)) # exp(-i*phi)
lp=np.zeros((d2,d2)) # exp(i*phi)
lz=np.zeros((d2,d2)) # angular momentum
lzSquare=np.zeros((d2,d2)) # angular momentum squared
# function_phi=np.zeros((d2,d2),dtype=complex) # uncomment to intialize any function

# Compute matrix elements
for m in range(d2):
    lz[m][m]=float(m_min+m)
    lzSquare[m][m]=float((m_min+m)*(m_min+m))
    if m!=0:
        lm[m-1][m]=1 #Decrease in momentum
        lp[m][m-1]=1 #Increase in momentum


# So qutip recognizes it     
Lp=Qobj(lp)      
Lm=Qobj(lm)
Lz=Qobj(lz)
LzSquare=Qobj(lzSquare)


#---------------HAMILTONIANS----------------
# System Hamiltonian
H = LzSquare/2.0/MI


# Lindblad Dissipators
C_1=sqrt(kh)*Lm
C_2=sqrt(kc)*Lp

#---------------TIMESTEPS AND INITIALISATION----------------
T=args.t_value
NT=args.nt_value
timesteps=np.linspace(0,T,NT)

options=Options()
options.nsteps=10e5
options.store_states=True
Qh = kh*(Lm*H*Lp - 0.5*Lm*Lp*H - 0.5*H*Lm*Lp)
Qc = kc*(Lp*H*Lm - 0.5*Lp*Lm*H - 0.5*H*Lp*Lm)

# Solution to master equation
solution=mesolve(H, rho, timesteps, [C_1,C_2],[Lz,LzSquare,Qh,Qc],options=options)
ELz = solution.expect[0]
ELz2 = solution.expect[1]
QH = solution.expect[2]
QC = solution.expect[3]
W = np.zeros(NT)
dataThermo=np.vstack((timesteps,ELz,ELz2)) 
p1 = np.zeros(NT)
p2 = np.zeros(NT)
KE = lz*lz/2.0/MI
for j in range(0,NT):
    rho = solution.states[j].full()
    W[j] = extractWork(rho,KE)
    p1[j] = rho[0][0]
    p2[j] = rho[d2-1][d2-1]
dataThermo=np.vstack((timesteps,ELz,ELz2,W,QH,QC,p1,p2)) 
file_data_store('Approximate_T%s_NT%s_mmin%s_mmax%s_I%.3f_g%.3f_kappa%.3f_nH%.2f_nC%.2f.dat'%(T,NT,m_min,m_max,MI,g,k1,n1,n3),dataThermo,numtype="real")

