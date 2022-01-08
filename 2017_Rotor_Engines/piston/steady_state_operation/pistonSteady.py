# cd Documents/CQT/Quantum Thermodynamics/Rotor

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
parser.add_argument("-nt", "--nt", dest="nt_value", type=int, metavar='<int>', default=401, help="NT value (default=201)")
parser.add_argument("-mmin", "--mmin", dest="mmin_value", type=int, metavar='<int>', default=-150, help="MMIN value (default=-20)")
parser.add_argument("-mmax", "--mmax", dest="mmax_value", type=int, metavar='<int>', default=400, help="MMAX value (default=60)")
parser.add_argument("-n1", "--n1", dest="n1_value", type=float, metavar='<float>', default=1.0, help="N1 value (default=0.33)")
parser.add_argument("-n2", "--n2", dest="n2_value", type=float, metavar='<float>', default=0.0, help="N2 value (default=0.0)")
parser.add_argument("-g", "--g", dest="g_value", type=float, metavar='<float>', default=3.0, help="G value (default=0.0)")
parser.add_argument("-k", "--k", dest="k_value", type=float, metavar='<float>', default=1.0, help="K value (default=10.0)")
parser.add_argument("-I", "--I", dest="I_value", type=float, metavar='<float>', default=10.0, help="I value (default=10.0)")
parser.add_argument("-gamma", "--gamma", dest="gamma_value", type=float, metavar='<float>', default=1.0, help="gamma value (default=10.0)")
parser.add_argument("-temp", "--temp", dest="temp_value", type=float, metavar='<float>', default=1.0, help="temp value (default=10.0)")


args = parser.parse_args()
#------PARAMETERS----------

# Dimension of the driving mode
d1=2

# Mean occupation number of hot bath
n1=args.n1_value

# Mean occupation number of cold bath
n2=args.n2_value

# Initial state of the mode (starts off cold)
s1 = np.zeros((2,2),dtype=complex) 
s1[0][0] = (n1+1.0)*1.0/(2.0*n1+1.0) 
s1[1][1] = n1*1.0/(2.0*n1+1.0) 
s1 = Qobj(s1)
sigmap  = np.zeros((2,2),dtype=complex) 
sigmap[1][0] = 1
sigmap = Qobj(sigmap)
sigmam  = np.zeros((2,2),dtype=complex) 
sigmam[0][1] = 1
sigmam = Qobj(sigmam)

# Maximum angular momentum quantum number of the rotor
# 40 to 80 or 20 to 250
m_max=args.mmax_value
m_min=args.mmin_value

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

# dissipation rate
# 1 or 10
gamma= args.gamma_value
temp = args.temp_value
# Timesteps
T=args.t_value
NT=args.nt_value
timesteps=np.linspace(0,T,NT)

options=Options()
options.nsteps=10e5
options.store_states=True
options.tol = 1e-16

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
rho=tensor(s1,Rot)

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
CosPhi = tensor(qeye(d1),cos_Phi)
SinPhi = tensor(qeye(d1),sin_Phi)


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
ch = tensor(sigmam,fHot_Phi)
chd = tensor(sigmap,fHot_Phi)
cc = tensor(sigmam,fCold_Phi)
ccd = tensor(sigmap,fCold_Phi)
C_hot=sqrt(k1*(1+n1))*tensor(sigmam,fHot_Phi)
C_hot_d=sqrt(k1*(n1))*tensor(sigmap,fHot_Phi)
C_cold=sqrt(k2*(1+n2))*tensor(sigmam,fCold_Phi)
C_cold_d=sqrt(k2*(n2))*tensor(sigmap,fCold_Phi)
ax = (cos_Phi - 1j*1.0/4.0/temp/MI * sin_Phi *Lz)
ay = (sin_Phi + 1j*1.0/4.0/temp/MI * cos_Phi *Lz)
r1 = tensor(qeye(2),ax)
r1d = r1.dag()
r2 = tensor(qeye(2),ay)
r2d = r2.dag()
#------- THERMOqDYNAMIC OPERATORS ----------

# Work rate (power)
P=(1.0/2.0)*(g/MI)*(tensor(sigmap*sigmam,sin_Phi*Lz+Lz*sin_Phi))

# Heat rate (in units of mode frequency)
fHotSquare_Phi=(qeye(d2)+2*sin_Phi+sinSquare_Phi)/4
QH=k1*(n1*tensor(qeye(d1),qeye(d2))-tensor(sigmap*sigmam,qeye(d2)))*tensor(qeye(d1),fHotSquare_Phi)
fColdSquare_Phi=(qeye(d2)-2*sin_Phi+sinSquare_Phi)/4
QC=k2*(n2*tensor(qeye(d1),qeye(d2))-tensor(sigmap*sigmam,qeye(d2)))*tensor(qeye(d1),fColdSquare_Phi)
Wload = r1d*H*r1 - 0.5*r1d*r1*H - 0.5*H*r1d*r1 + r2d*H*r2 - 0.5*r2d*r2*H - 0.5*H*r2d*r2
QHload = k1*(n1+1)*(chd*H*ch - 0.5*H*chd*ch - 0.5*chd*ch*H) +k1*n1*(ch*H*chd - 0.5*H*ch*chd - 0.5*ch*chd*H) 
QCload = k2*(n2+1)*(ccd*H*cc - 0.5*H*ccd*cc - 0.5*ccd*cc*H) +k2*n2*(cc*H*ccd - 0.5*H*cc*ccd - 0.5*cc*ccd*H) 
T = tensor(sigmap*sigmam,cos_Phi)
#gamma array
NT = 1
gamma = np.zeros(NT)
for m in range(0,NT):
    gamma[m]=10.0**(args.gamma_value/10.0)
#--------SIMULATION----------
tmpW = np.zeros(NT,dtype=complex)
tmpW1 = np.zeros(NT,dtype=complex)
tmpW2 = np.zeros(NT,dtype=complex)
tmpWload = np.zeros(NT,dtype=complex)
tmpQHload = np.zeros(NT,dtype=complex)
tmpQCload = np.zeros(NT,dtype=complex)
tmpQH = np.zeros(NT,dtype=complex)
tmpQC = np.zeros(NT,dtype=complex)
tmpELz = np.zeros(NT,dtype=complex)
tmpELz2 = np.zeros(NT,dtype=complex)
tmpEN = np.zeros(NT,dtype=complex)
tmpCos = np.zeros(NT,dtype=complex)
tmpSin = np.zeros(NT,dtype=complex)
tmp1 = np.zeros(NT,dtype=complex)
p1 = np.zeros(NT,dtype=complex)
p2 = np.zeros(NT,dtype=complex)
for kk in range(0,NT):
    A1 = sqrt(2*temp*gamma[kk]*MI)*r1
    A2 = sqrt(2*temp*gamma[kk]*MI)*r2
    state=steadystate(H, [C_hot,C_hot_d,C_cold,C_cold_d,A1,A2],tol=1e-12)
    tmpW2[kk] = expect(P,state)*1000000
    tmpQH[kk] = expect(QH,state)*1000000
    tmpQC[kk] = expect(QC,state)*1000000
    tmpELz[kk] = expect(ELz,state)*1000000
    tmpELz2[kk] = expect(ELz2,state)*1000000
    tmpEN[kk] = expect(EN,state)*1000000
    tmpCos[kk] = expect(CosPhi,state)*1000000
    tmpSin[kk] = expect(SinPhi,state)*1000000
    tmpQCload[kk] = expect(QCload,state)*1000000
    tmpQHload[kk] = expect(QHload,state)*1000000
    tmpWload[kk] = expect(Wload,state)*(2*temp*gamma[kk]*MI)*1000000
    tmp1[kk] = expect(T,state)*g/16/temp/MI*1000000
    Hke = lz*lz/2/MI
    rho = state.full()
    rhor = state.ptrace(1).full()
    tmpW1[kk] = extractWork(rhor,Hke)
    p1[kk] = rho[0][0]
    p2[kk] = rho[d2-1][d2-1]
    for k in range(0,d1):
        rhotmp = [[rho[i][m] for m in range(k*d2,k*d2+d2)] for i in range(k*d2,k*d2+d2)]
        Htmp = Hke + g*cos_phi*k
        tmpW[kk] = tmpW[kk] + extractWork(rhotmp,Htmp)




dataThermo=np.vstack((gamma,tmpW,tmpW1,tmpW2,tmpQH,tmpQC,tmpWload,tmpQCload,tmpQHload,tmpELz,tmpELz2,tmpEN,tmpCos,tmpSin,tmp1,p1,p2))
file_data_store('I%s_g%s_kappa%s_Gamma%s_nH%s_nC%s_T%s.dat'%(m_min,m_max,MI,g,k1,args.gamma_value/10.0,n1,n2,temp),dataThermo,numtype="real")
