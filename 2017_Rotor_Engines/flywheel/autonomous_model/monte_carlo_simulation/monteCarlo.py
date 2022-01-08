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
parser.add_argument("-T", "--T", dest="t_value", type=float, metavar='<float>', default=500.0, help="T value (default=10.0)")
parser.add_argument("-NT", "--NT", dest="nt_value", type=int, metavar='<int>', default=1001, help="NT value (default=51.0)")
parser.add_argument("-th", "--th", dest="th_value", type=float, metavar='<float>', default=5.0, help="TH value (default=5.0)")
parser.add_argument("-tc", "--tc", dest="tc_value", type=float, metavar='<float>', default=0.0, help="TC value (default=1.0)")
parser.add_argument("-wh", "--wh", dest="wh_value", type=float, metavar='<float>', default=5.0, help="wh value (default=5.0)")
parser.add_argument("-wc", "--wc", dest="wc_value", type=float, metavar='<float>', default=3.0, help="wc value (default=3.0)")
parser.add_argument("-mmin", "--mmin", dest="mmin_value", type=int, metavar='<int>', default=-10, help="MMIN value (default=-10)")
parser.add_argument("-mmax", "--mmax", dest="mmax_value", type=int, metavar='<int>', default=30, help="MMAX value (default=30)")
parser.add_argument("-g", "--g", dest="g_value", type=float, metavar='<float>', default=1.0, help="G value (default=5.0)")
parser.add_argument("-k", "--k", dest="k_value", type=float, metavar='<float>', default=1.0, help="K value (default=100.0)")
parser.add_argument("-I", "--I", dest="I_value", type=float, metavar='<float>', default=1.0, help="I value (default=10.0)")
parser.add_argument("-run", "--run", dest="run_value", type=float, metavar='<int>', default=0, help="run value (default=0)")

# parser.add_argument("-delta", "--delta", dest="delta_value", type=float, metavar='<float>', default=0.0, help="DELTA value (default=0.0)")
args = parser.parse_args()

# For engine to work, TH/TC > wh/wc > 1
# Dimension of rotor

#---------------TWO QUBITS---------------
# Energy of the qubits
wh = args.wh_value
wc = args.wc_value
# delta = args.delta_value
# T2 = wh-wc -delta

# Temperature of baths
TH = args.th_value #Qubit 1 : Get from argument
TC = args.tc_value #Qubit 3

# Mean Occupation of the qubits
if (TH != 0):
    n1 = exp(-wh/TH)  #Ratio of finding qubit1 in excited state : ground state
    ph = n1*1.0/(1.0-n1)
else: 
    n1 = 0.0
    ph = 0.0
if (TC != 0):
    n3 = exp(-wc/TC)  #Ratio of finding qubit3 in excited state : ground state
    pc = n3*1.0/(1.0-n3)
else: 
    n3 = 0.0
    pc = 0.0

ground = np.zeros((2,1),dtype=complex)
ground[0] = 1
ground = Qobj(ground)
excited = np.zeros((2,1),dtype=complex)
excited[1] = 1
excited = Qobj(excited)

#---------------ROTOR----------------
# Maximum angular momentum quantum number of the rotor
# 40 to 80 or 20 to 250
m_max=args.mmax_value
m_min=args.mmin_value

# Dimension of rotor Hilbert space
d2=(m_max-m_min)+1

# Number of trajectories
N = 2
run = args.run_value

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
tau2 = rot_ket


# Coupling strength
k1 = args.k_value    #effective thermalisation rate for qubit 1 (high)
k3 = args.k_value    #effective thermalisation rate for qubit 3 (high)
# k2 = args.k_value/10.0
# Interaction Strength 
g = args.g_value


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
h1 = np.zeros((2,2),dtype=complex) 
h1[1][1] = wh
h2 = LzSquare/2.0/MI
h3 = np.zeros((2,2),dtype=complex)
h3[1][1] = wc
h1 = Qobj(h1)
h1 = tensor(h1,qeye(d2),qeye(2)) 
h2 = tensor(qeye(2),h2,qeye(2))
h3 = Qobj(h3)
h3 = tensor(qeye(2),qeye(d2),h3) 
H0 = h1 + h2 + h3

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
Hint = g*(tensor(sigmam,Lp,sigmap)+tensor(sigmap,Lm,sigmam))

# Total Hamiltonian
H = H0 + Hint

# Lindblad Dissipators
C_1=sqrt(k1)*tensor(sigmam,qeye(d2),qeye(2))
C_1_d=sqrt(k1*n1)*tensor(sigmap,qeye(d2),qeye(2))
# C_2=sqrt(k2*(n2+1))*tensor(qeye(2),Lm,qeye(2))
# C_2_d=sqrt(k2*(n2))*tensor(qeye(2),Lp,qeye(2))
C_3=sqrt(k3)*tensor(qeye(2),qeye(d2),sigmam)
C_3_d=sqrt(k3*n3)*tensor(qeye(2),qeye(d2),sigmap)

#---------------TIMESTEPS AND INITIALISATION----------------
T=args.t_value
NT=args.nt_value
timesteps=np.linspace(0,T,NT)

options=Options()
options.nsteps=1e5
options.ntraj= N
options.store_states=True


#---------------EXPECTATION VALUES----------------
a1 = tensor(sigmam,qeye(d2),qeye(2))
ad1 = tensor(sigmap,qeye(d2),qeye(2))
a3 = tensor(qeye(2),qeye(d2),sigmam)
ad3 = tensor(qeye(2),qeye(d2),sigmap)
Qh = k1*(ad1*H*a1 - 0.5*ad1*a1*H - 0.5*H*ad1*a1) + k1*n1*(a1*H*ad1 - 0.5*a1*ad1*H - 0.5*H*a1*ad1)
Qc = k3*(ad3*H*a3 - 0.5*ad3*a3*H - 0.5*H*ad3*a3) + k3*n3*(a3*H*ad3 - 0.5*a3*ad3*H - 0.5*H*a3*ad3)
ELz = tensor(qeye(2),Lz,qeye(2))
ELz2 = tensor(qeye(2),LzSquare,qeye(2))
W = np.zeros(NT)
W1 = np.zeros(NT)
p1 = np.zeros(NT)
p2 = np.zeros(NT)
# Solution to master equation
if (TC == 0):
    ket1 = tensor(ground,tau2,ground)
    bra1 = ket1.dag()
    ket2 = tensor(excited,tau2,ground)
    solution1=mcsolve(H, ket1, timesteps, [C_1,C_1_d,C_3,C_3_d],[Qh,Qc,ELz,ELz2],options=options)
    solution2=mcsolve(H, ket2, timesteps, [C_1,C_1_d,C_3,C_3_d],[Qh,Qc,ELz,ELz2],options=options)
    Qh = solution1.expect[0]*(1-ph) + solution2.expect[0]*(ph) 
    Qc = solution1.expect[1]*(1-ph) + solution2.expect[1]*(ph) 
    ELz = solution1.expect[2]*(1-ph) + solution2.expect[2]*(ph) 
    ELz2 = solution1.expect[3]*(1-ph) + solution2.expect[3]*(ph) 
    KE = lz*lz/2.0/MI
    for j in range(0,NT):
        state1 = np.zeros((2*d2*2,2*d2*2),dtype=complex)
        state1 = Qobj(state1)
        state1.dims = [[2,d2,2],[2,d2,2]]
        print(state1)
        state2 = state1
        for k in range(0,N):
            r1 = solution1.states[k][j]
            print(r1)
            r2 = solution2.states[k][j]
            state1 += r1*r1.dag() 
            state2 += r2*r2.dag()
        state = state1*(1-ph)/(N*1.0) + state2*ph/(N*1.0)
        rhor = state.ptrace(1).full()
        rho = state.full()
        W[j] = getErgotropy(rhor,KE)
        W1[j] = getErgotropy(rho,H)
        p1[j] = rhor[0][0]
        p2[j] = rhor[d2-1][d2-1]
else:   
    ket1 = tensor(ground,tau2,ground)
    ket2 = tensor(excited,tau2,ground)
    ket3 = tensor(ground,tau2,excited)
    ket4 = tensor(excited,tau2,excited)
    solution1=mcsolve(H, ket1, timesteps, [C_1,C_1_d,C_3,C_3_d],[Qh,Qc,ELz,ELz2],options=options)
    solution2=mcsolve(H, ket2, timesteps, [C_1,C_1_d,C_3,C_3_d],[Qh,Qc,ELz,ELz2],options=options)
    solution3=mcsolve(H, ket3, timesteps, [C_1,C_1_d,C_3,C_3_d],[Qh,Qc,ELz,ELz2],options=options)
    solution4=mcsolve(H, ket4, timesteps, [C_1,C_1_d,C_3,C_3_d],[Qh,Qc,ELz,ELz2],options=options)
    Qh = solution1.expect[0]*(1-ph)*(1-pc) + solution2.expect[0]*(ph)*(1-pc) + solution3.expect[0]*(1-ph)*(pc) + solution4.expect[0]*(ph)*(pc) 
    Qc = solution1.expect[1]*(1-ph)*(1-pc) + solution2.expect[1]*(ph)*(1-pc) + solution3.expect[1]*(1-ph)*(pc) + solution4.expect[1]*(ph)*(pc) 
    ELz = solution1.expect[2]*(1-ph)*(1-pc) + solution2.expect[2]*(ph)*(1-pc) + solution3.expect[2]*(1-ph)*(pc) + solution4.expect[2]*(ph)*(pc) 
    ELz2 = solution1.expect[3]*(1-ph)*(1-pc) + solution2.expect[3]*(ph)*(1-pc) + solution3.expect[3]*(1-ph)*(pc) + solution4.expect[3]*(ph)*(pc) 
    KE = lz*lz/2.0/MI
    for j in range(0,NT):
        ket1 = solution1.states[j]
        bra1 = ket1[0].dag()
        ket2 = solution2.states[j]
        bra2 = ket2[0].dag()
        ket3 = solution3.states[j]
        bra3 = ket3[0].dag()
        ket4 = solution4.states[j]
        bra4 = ket4[0].dag()
        state = ket1*bra1*(1-ph)*(1-pc) + ket2*bra2*ph*(1-pc) +ket3*bra3*(1-ph)*(pc) + ket4*bra2*ph*pc
        rhor = state[0].ptrace(1).full()
        rho = state[0].full()
        W[j] = extractWork(rhor,KE)
        W1[j] = extractWork(rho,H)
        p1[j] = rhor[0][0]
        p2[j] = rhor[d2-1][d2-1]


dataThermo=np.vstack((timesteps,ELz,ELz2,W,W1,Qh,Qc,p1,p2)) 
# file_data_store('MCT%s_NT%s_mmin%s_mmax%s_I%.3f_g%.3f_kappa%.3f_TH%.2f_TC%.2f_wh%.2f_wc%.2f.dat'%(T,NT,m_min,m_max,MI,g,k1,TH,TC,wh,wc),dataThermo,numtype="real")
file_data_store('MCrun%s_T%s_NT%s_mmin%s_mmax%s_I%.3f_g%.3f_kappa%.3f_TH%.2f_TC%.2f_wh%.2f_wc%.2f.dat'%(run,T,NT,m_min,m_max,MI,g,k1,TH,TC,wh,wc),dataThermo,numtype="real")

file_data_store('MCStaterun%sT%s_NT%s_mmin%s_mmax%s_I%.3f_g%.3f_kappa%.3f_TH%.2f_TC%.2f_wh%.2f_wc%.2f.dat'%(run,T,NT,m_min,m_max,MI,g,k1,TH,TC,wh,wc),rho)
file_data_store('MCRotorStaterun%sT%s_NT%s_mmin%s_mmax%s_I%.3f_g%.3f_kappa%.3f_TH%.2f_TC%.2f_wh%.2f_wc%.2f.dat'%(run,T,NT,m_min,m_max,MI,g,k1,TH,TC,wh,wc),rhor)
