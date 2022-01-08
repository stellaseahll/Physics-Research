# Imports all qutip functions
from qutip import *
from qutip.piqs import *
# For sqrt, exp and other functions
from math import exp
from math import pi
from math import sin
from math import cos
from math import sqrt
from math import log10

import numpy as np
from scipy.linalg import *

# NOTE: hbar=1 throughout
#	   Indices 0,1,2,...,m_min+m_max correspond to -m_min,-m_min+1,...,0,..., m_max-1,m_max


#N1 = 3
#N2 = 3
Nlist = ( (3,3), (3,1), (1,3), (1,1) )
gmax = 0.1
gmin = 0.0001
Ng = 101
kh = 0.001
kc = 0.001
Th = 1.0
Tc = 0.1

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

nh = meanOccupation(1,Th)
nc = meanOccupation(1,Tc)

cQh = np.zeros([Ng,len(Nlist)])
cQc = np.zeros([Ng,len(Nlist)])
gvals = np.logspace(log10(gmin),log10(gmax),Ng)

for j in range(0,len(Nlist)):
    N1 = Nlist[j][0]
    N2 = Nlist[j][1]
    # Spin operators
    id1 = qeye(N1+1)
    id2 = qeye(N2+1)
    sp1 = jmat(N1/2,'+')
    sm1 = jmat(N1/2,'-')
    sz1 = jmat(N1/2,'z')
    sx1 = jmat(N1/2,'x')
    sy1 = jmat(N1/2,'y')
    sp2 = jmat(N2/2,'+')
    sm2 = jmat(N2/2,'-')
    sz2 = jmat(N2/2,'z')
    sx2 = jmat(N2/2,'x')
    sy2 = jmat(N2/2,'y')
    # Hamiltonians
    h1 = tensor(sz1,id2)
    h2 = tensor(id1,sz2)
    for k in range(0,Ng):
        g = gvals[k]
        Hint = g*(tensor(sp1,sm2)+tensor(sm1,sp2))
        H = h1 + h2 + Hint
        # Lindblad Dissipators
        sphm = sqrt(kh*(nh+1))*tensor(sm1,id2)
        sphp = sqrt(kh*nh)*tensor(sp1,id2)
        spcm = sqrt(kc*(nc+1))*tensor(id1,sm2)
        spcp = sqrt(kc*nc)*tensor(id1,sp2)
    
        rhoSS = steadystate(H,[sphm,sphp,spcm,spcp])
        spQHop = getQOp(H,[sphm,sphp])
        spQCop = getQOp(H,[spcm,spcp])
        
        cQh[k,j] = expect(spQHop,rhoSS)
        cQc[k,j] = expect(spQCop,rhoSS)
        

data = np.hstack((cQh,cQc))
