# Imports all qutip functions
from qutip import *
from qutip.piqs import *
# For sqrt, exp and other functions
from math import exp
from math import pi
from math import sin
from math import cos
from math import sqrt

import numpy as np
from scipy.linalg import *

# NOTE: hbar=1 throughout
#	   Indices 0,1,2,...,m_min+m_max correspond to -m_min,-m_min+1,...,0,..., m_max-1,m_max

##############################
## Parse Arguments
#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-th", "--th", dest="th_value", type=float, metavar='<float>', default=1.0, help="TH value (default=5.0)")
parser.add_argument("-tc", "--tc", dest="tc_value", type=float, metavar='<float>', default=0.1, help="TC value (default=1.0)")
parser.add_argument("-kh", "--kh", dest="kh_value", type=float, metavar='<float>', default=0.001, help="kh value (default=0.001)")
parser.add_argument("-kc", "--kc", dest="kc_value", type=float, metavar='<float>', default=0.001, help="kc value (default=0.001)")
parser.add_argument("-g", "--g", dest="g_value", type=float, metavar='<float>', default=0.01, help="g value (default=0.001)")
parser.add_argument("-N1", "--N1", dest="N1_value", type=int, metavar='<int>', default=3, help="N1 value (default=1.0)")
parser.add_argument("-N2", "--N2", dest="N2_value", type=int, metavar='<int>', default=3, help="N2 value (default=1.0)")

args = parser.parse_args()

def meanOccupation(E,Tr):
	if (Tr != 0):
		n = 1.0/(exp(E/Tr)-1)
		return n
	n = 0
	return n

def leftMultiply(A):
	d = A.dims[1]
	return tensor(qeye(d),A)

def rightMultiply(A):
	d = A.dims[1]
	return tensor(A.trans(),qeye(d))

def lrMultiply(A):
	return tensor(A.conj(),A)

def liouvilleH(H):
	return 1j*(rightMultiply(H)-leftMultiply(H))

def liouvilleDiss(A):
	B = A.dag()*A
	return lrMultiply(A)-0.5*leftMultiply(B)-0.5*rightMultiply(B)

def createLiouville(H,Cops):
	L = liouvilleH(H)
	for j in range(0,len(Cops)):
		L = L + liouvilleDiss(Cops[j])
	return L

def vec2state(rho,d):
	tmp = np.reshape(rho,(d,d))
	tmp = 0.5*(tmp+tmp.conj().T)
	return tmp/np.trace(tmp)

def getQOp(H,Cops):
	L = H*0
	for j in range(0,len(Cops)):
		b = Cops[j].dag()*Cops[j]
		L = L + Cops[j].dag()*H*Cops[j]- 0.5*(H*b + b*H)
	return L

#------PARAMETERS----------

# Spin operators 1st chain
N1 = args.N1_value
sp1 = jspin(N1,'+')
sm1 = jspin(N1,'-')
sz1 = jspin(N1,'z')
sx1 = jspin(N1,'x')
sy1 = jspin(N1,'y')
d1 = sp1.dims[0][0]

# Spin operators 2nd chain
N2 = args.N2_value
sp2 = jspin(N2,'+')
sm2 = jspin(N2,'-')
sz2 = jspin(N2,'z')
sx2 = jspin(N2,'x')
sy2 = jspin(N2,'y')
d2 = sp2.dims[0][0]

# Total dimension
dtot = d1*d2

# Hamiltonians
h1 = tensor(sz1,qeye(d2))
h2 = tensor(qeye(d1),sz2)
g = 10**(args.g_value*0.01-4)
Hint = g*(tensor(sp1,sm2)+tensor(sm1,sp2))
H = h1 + h2 + Hint

# Dissipator coefficients
Th = args.th_value
Tc = args.tc_value
kh = args.kh_value
kc = args.kc_value
nh = meanOccupation(1,Th)
nc = meanOccupation(1,Tc)

# Lindblad Dissipators
chm = sqrt(kh*(nh+1))*tensor(sm1,qeye(d2))
chp = sqrt(kh*nh)*tensor(sp1,qeye(d2))
ccm = sqrt(kc*(nc+1))*tensor(qeye(d1),sm2)
ccp = sqrt(kc*nc)*tensor(qeye(d1),sp2)

# Liouville Superoperator
cL = createLiouville(H,[chm,chp,ccm,ccp])

# Solve for null vectors
M = null_space(cL.full())
row,col = M.shape
cP = np.zeros((dtot**2,dtot**2))
for j in range(0,col):
	v = M[:,j:j+1]
	cP = cP + np.dot(v,(v.conj().T))

cQHop = getQOp(H,[chm,chp])
cQCop = getQOp(H,[ccm,ccp])
nSS1 = int(np.ceil((N1+1)/2))
nSS2 = int(np.ceil((N2+1)/2))
cQh = np.zeros(nSS1*nSS2+1)
cQc = np.zeros(nSS1*nSS2+1)
cJ1 = np.zeros(nSS1*nSS2+1)
cJ2 = np.zeros(nSS1*nSS2+1)
# for j in range(0,dtot):
# 	rho = np.zeros((dtot**2,1))
# 	rho[j*(dtot+1)][0] = 1
# 	rhoSS = Qobj(vec2state(np.dot(P,rho),dtot))
# 	rhoSS.dims = H.dims
# 	Qh[j] = expect(QHop,rhoSS)
# 	Qc[j] = expect(QCop,rhoSS)
idx = 0
tol = 1e-12
for j1 in np.arange(0.5*N1,-0.5,-1):
	for j2 in np.arange(0.5*N2,-0.5,-1):
		# 	for m2 in np.arange(-j2,j2+0.5,1):
		rho = tensor(dicke(N1,j1,j1),dicke(N2,j2,j2)).full()
		rho = np.reshape(rho,(dtot**2,1))
		rhoSS = Qobj(vec2state(np.dot(cP,rho),dtot))
		rhoSS.dims = H.dims
		cJ1[idx] = j1
		cJ2[idx] = j2
		cQh[idx] = expect(cQHop,rhoSS)
		cQc[idx] = expect(cQCop,rhoSS)
		if (abs(cQh[idx])<tol or abs(cQc[idx])<tol):
			cQh[idx] = 0
			cQc[idx] = 0
		idx = idx+1

# Spin operators qubit
sp = jmat(0.5,'+')
sm = jmat(0.5,'-')
sz = jmat(0.5,'z')
sx = jmat(0.5,'x')
sy = jmat(0.5,'y')

# Hamiltonians
h1 = tensor(sz,qeye(2))
h2 = tensor(qeye(2),sz)
Hint = g*(tensor(sp,sm)+tensor(sm,sp))
H = h1 + h2 + Hint

# Lindblad Dissipators
sphm = sqrt(kh*(nh+1))*tensor(sm,qeye(2))
sphp = sqrt(kh*nh)*tensor(sp,qeye(2))
spcm = sqrt(kc*(nc+1))*tensor(qeye(2),sm)
spcp = sqrt(kc*nc)*tensor(qeye(2),sp)

# # Liouville Superoperator
# spL = createLiouville(H,spops)
# Solve for null vectors
# M = null_space(spL.full())
# rhoSS = Qobj(vec2state(M[:,0],4))
# rhoSS.dims = H.dims
rhoSS = steadystate(H,[sphm,sphp,spcm,spcp])
spQHop = getQOp(H,[sphm,sphp])
spQCop = getQOp(H,[spcm,spcp])
cQh[nSS1*nSS2] = expect(spQHop,rhoSS)
cQc[nSS1*nSS2] = expect(spQCop,rhoSS)
data = np.vstack((cJ1,cJ2,cQh,cQc))
file_data_store('Nh%s_Nc%s_kh%.3f_kc%.3f_g%.4f_Th%.3f_Tc%.3f.dat'%(N1,N2,kh,kc,g,Th,Tc),data,numtype="real")
