# Imports all qutip functions
from qutip import *
# For sqrt, exp and other functions
from math import exp
from math import pi
from math import sin
from math import cos
from math import sqrt
from extractWork import extractWork
import scipy.integrate as integrate
import numpy as np

#For tic-toc
import time


print(time.ctime())
tic = time.clock()

# NOTE: hbar=1 throughout
#       Indices 0,1,2,...,m_min+m_max correspond to -m_min,-m_min+1,...,0,..., m_max-1,m_max

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--t", dest="t_value", type=int, metavar='<int>', default=100, help="T value (default=100)")
parser.add_argument("-nt", "--nt", dest="nt_value", type=int, metavar='<int>', default=501, help="NT value (default=201)")
parser.add_argument("-mmin", "--mmin", dest="mmin_value", type=int, metavar='<int>', default=-40, help="MMIN value (default=-20)")
parser.add_argument("-mmax", "--mmax", dest="mmax_value", type=int, metavar='<int>', default=100, help="MMAX value (default=60)")
parser.add_argument("-n1", "--n1", dest="n1_value", type=float, metavar='<float>', default=1.0, help="N1 value (default=0.33)")
parser.add_argument("-n2", "--n2", dest="n2_value", type=float, metavar='<float>', default=0.0, help="N2 value (default=0.0)")
parser.add_argument("-g", "--g", dest="g_value", type=float, metavar='<float>', default=1.0, help="G value (default=3.0)")
parser.add_argument("-k", "--k", dest="k_value", type=float, metavar='<float>', default=1.0, help="K value (default=10.0)")
parser.add_argument("-I", "--I", dest="I_value", type=float, metavar='<float>', default=10.0, help="I value (default=10.0)")
parser.add_argument("-run", "--run", dest="run_value", type=int, metavar='<int>', default=1, help="run value (default=1)")

args = parser.parse_args()
#------PARAMETERS----------
run = args.run_value
# Dimension of the driving mode
d1=2

# Mean occupation number of hot bath
n1=args.n1_value

# Mean occupation number of cold bath
n2=args.n2_value

# Initial state of the mode (starts off cold)
s1 = np.zeros((2,2),dtype=complex) 
s1[0][0] = 1
s1[1][1] = 0 
s1 = Qobj(s1)
sigmap  = np.zeros((2,2),dtype=complex) 
sigmap[1][0] = 1
sigmap = Qobj(sigmap)
sigmam  = np.zeros((2,2),dtype=complex) 
sigmam[0][1] = 1
sigmam = Qobj(sigmam)

# Maximum angular momentum quantum number of the rotor
# 40 to 80 or 20 to 250
m_max=args.mmax_value + (run-1)*10
m_min=args.mmin_value + (run-1)*10

# Dimension of rotor Hilbert space
d2=(m_max-m_min)+1

# Kappa parameter of von Mises distribution (The larger it is the more peaked the distribution will be) 
kappa=10

# Mean of von Mises distribution (Initial expectation value of the angle) 
mu=pi/2

# Momentum of inertia (If too small, dispersion occurs too quickly)
#10 or 100
MI=args.I_value

# Interaction strength
# 1 or 10
g=args.g_value

# Dissipation strengths of collapse operators
k1=args.k_value #hot bath
k2=args.k_value #cold bath

# Timesteps
T=args.t_value
NT=args.nt_value
dt = T*1.0/(NT*1.0)
timesteps=np.linspace(0,T,NT)

options=Options()
options.nsteps=1e7
options.atol=1e-12
options.rtol=1e-10
options.store_states=True


#-----ROTOR---------

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

def shiftMomentum(x,d1,d2,step):
    totalDim = d1*d2
    # Do rows first
    for i in range(0,totalDim):
        currRemainder = i%d2
        currQuotientIndex = (i/d2)*d2
        if (currRemainder >= step):
            x[currRemainder-step + currQuotientIndex] = x[currRemainder + currQuotientIndex]
        if (d2-currRemainder <= step):
            x[currRemainder + currQuotientIndex] = 0
            
    # Now do column
    for i in range(0,totalDim):
        currRemainder = i%d2
        currQuotientIndex = (i/d2)*d2
        if (currRemainder >= step):
            x[:,currRemainder-step + currQuotientIndex] = x[:,currRemainder + currQuotientIndex]
        if (d2-currRemainder <= step):
            x[:,currRemainder + currQuotientIndex] = 0
    return x
            

# Initial state
if (run == 1):
    rho = tensor(s1,Rot)
else:
    tmp = np.loadtxt('kstate%s.dat'%(run-1),dtype=complex, delimiter = ',')
    tmp = shiftMomentum(tmp,d1,d2,10)
    rho = Qobj(tmp)
    rho.dims = [[d1,d2],[d1,d2]]


#-------OPERATORS DEFINITION---------
# Integral for the matrix elements of any function of the angle operator
# def int_function(m,n):
#     I1=integrate.quad(lambda x: (1/(2*pi))*cos(((-m_max+n)-(-m_max+m))*x)*function(x),0,2*pi,limit=m_max*100)
#     I2=integrate.quad(lambda x: (1/(2*pi))*sin(((-m_max+n)-(-m_max+m))*x)*function(x),0,2*pi,limit=m_max*100)
#     return I1[0]+1j*I2[0]

# Initialize matrix with zeros
cos_phi=np.zeros((d2,d2),dtype=complex) # cos
sin_phi=np.zeros((d2,d2),dtype=complex) # sin
sinSquare_phi=np.zeros((d2,d2),dtype=complex) # sin^2
lz=np.zeros((d2,d2)) # angular momentum
lzSquare=np.zeros((d2,d2)) # angular momentum squared
# function_phi=np.zeros((d2,d2),dtype=complex) # uncomment to intialize any function

# Compute matrix elements
for m in range(d2):
    lz[m][m]=float(m_min+m)
    lzSquare[m][m]=float((m_min+m)*(m_min+m))
    if m!=0:
        cos_phi[m-1][m]=0.5+1j*0.0
        cos_phi[m][m-1]=0.5+1j*0.0
        sin_phi[m-1][m]=0.0+1j*0.5
        sin_phi[m][m-1]=0.0-1j*0.5
        sinSquare_phi[m][m]=0.5+1j*0.0
        sinSquare_phi[m-2][m]=-0.25+1j*0.0
        sinSquare_phi[m][m-2]=-0.25+1j*0.0

# So qutip recognizes it     
cos_Phi=Qobj(cos_phi)      
sin_Phi=Qobj(sin_phi)
sinSquare_Phi=Qobj(sinSquare_phi)
Lz=Qobj(lz)
LzSquare=Qobj(lzSquare)


#-------HAMILTONIANS---------

# Rotor Hamiltonian
H_r=tensor(qeye(d1),LzSquare/(2.0*MI))
ELz = tensor(qeye(d1),Lz)
ELz2 = tensor(qeye(d1),LzSquare)
EN = tensor(sigmap*sigmam,qeye(d2))
# Interaction Hamiltonian
Hint=g*tensor(sigmap*sigmam,cos_Phi)

# Total Hamiltonian
H=H_r+Hint


#-------COLLAPSE OPERATORS----------

# Using ((1+sin(x))/2) for modulating the hot bath (approx=ON for phi < pi, approx=OFF otherwise)
fHot_Phi=((qeye(d2)+sin_Phi)/2)

# Using ((1-sin(x))/2) for modulating the cold bath (approx=OFF for phi < pi, approx=ON otherwise)
fCold_Phi=((qeye(d2)-sin_Phi)/2)

# Lindblad collapse operators
C_hot=sqrt(k1*(1+n1))*tensor(sigmam,fHot_Phi)
C_hot_d=sqrt(k1*(n1))*tensor(sigmap,fHot_Phi)
C_cold=sqrt(k2*(1+n2))*tensor(sigmam,fCold_Phi)
C_cold_d=sqrt(k2*(n2))*tensor(sigmap,fCold_Phi)

#------- THERMODYNAMIC OPERATORS ----------

# Work rate (power)
P=(1.0/2.0)*(g/MI)*(tensor(sigmap*sigmam,sin_Phi*Lz+Lz*sin_Phi))

# Heat rate (in units of mode frequency)
fHotSquare_Phi=(qeye(d2)+2*sin_Phi+sinSquare_Phi)/4
QH=k1*(n1*tensor(qeye(d1),qeye(d2))-tensor(sigmap*sigmam,qeye(d2)))*tensor(qeye(d1),fHotSquare_Phi)
fColdSquare_Phi=(qeye(d2)-2*sin_Phi+sinSquare_Phi)/4
QC=k2*(n2*tensor(qeye(d1),qeye(d2))-tensor(sigmap*sigmam,qeye(d2)))*tensor(qeye(d1),fColdSquare_Phi)

#--------SIMULATION----------
# Inform that the real computation is starting
tocBeforeSolving = time.clock()
print('Initialization time = %.2f s'%(tocBeforeSolving-tic))   
print("Solving master equation...")

# Solution to master equation
solution=mesolve(H, rho, timesteps, [C_hot,C_hot_d,C_cold,C_cold_d], [P,QH,QC,ELz,ELz2, EN],options=options)

tocFinal = time.clock()
print('Total elapsed time = %.2f s'%(tocFinal-tic))


#--------SIMULATION OUTPUT----------

#Expectation values work rate
# lz=solution.expect[0]
# for j in range(0,NT):
#     file_data_store('T%s_NT%s_mMin%s_mMax%s_%.2fpi_k%s_I%s_g%s_Gamma%s_nH%s_nC%sts%s.dat'%(T,NT,m_min,m_max,mu/pi,kappa,MI,g,k1,n1,n2,j),solution.states[j])

# file_data_store('RotorThermo_T%s_NT%s_mMin%s_mMax%s_%.2fpi_k%s_I%s_g%s_Gamma%s_nH%s_nC%sLz.dat'%(T,NT,m_min,m_max,mu/pi,kappa,MI,g,k1,n1,n2),dataThermo,numtype="real")

rhor = solution.states[0].ptrace(1)
lzdist = np.hstack(rhor.diag())
for j in range(1,NT):
    rhor = solution.states[j].ptrace(1)
    tmpdist = rhor.diag()
    lzdist = np.vstack((lzdist,tmpdist))
    
W = np.zeros(NT)
W1 = np.zeros(NT)
W2 = solution.expect[0]
QH = solution.expect[1]
QC = solution.expect[2]
ELz = solution.expect[3]
ELz2 = solution.expect[4]
EN = solution.expect[5]
Hke = lz*lz/2/MI
p1 = np.zeros(NT)
p2 = np.zeros(NT)
for j in range(0,NT):
    rho = solution.states[j].full()
    rhor = solution.states[j].ptrace(1).full()
    W1[j] = extractWork(rhor,Hke)
    p1[j] = rho[0][0]
    p2[j] = rho[d2-1][d2-1]
    for k in range(0,d1):
        rhotmp = [[rho[i][m] for m in range(k*d2,k*d2+d2)] for i in range(k*d2,k*d2+d2)]
        Htmp = Hke + g*cos_phi*k
        W[j] = W[j] + extractWork(rhotmp,Htmp)
dataThermo=np.vstack((timesteps,W,W1,W2,QH,QC,ELz,ELz2,EN,p1,p2))
#     print(j)
#     x = solution.states[j]
#     y = x.ptrace(1)
#     W1[j] = extractWork(x,H)
#     W2[j] = extractWork(y,LzSquare/(2*MI))

# dataThermo=np.vstack((timesteps,W1,W2))
file_data_store('run%s_kappa%s_g%s_I%s_mmin%s_mmax%s.dat'%(run,k1,g,MI,m_min,m_max),dataThermo,numtype="real")
file_data_store('run%s_kappa%s_g%s_I%s_mmin%s_mmax%sState.dat'%(run,k1,g,MI,m_min,m_max),solution.states[NT-1],numtype="complex")
file_data_store('run%s_kappa%s_g%s_I%s_mmin%s_mmax%sState.dat'%(run,k1,g,MI,m_min,m_max),lzdist,numtype="real")

