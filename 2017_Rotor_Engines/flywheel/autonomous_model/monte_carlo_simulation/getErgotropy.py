from qutip import *

def getErgotropy(rho,H):
	rhotmp = Qobj(rho)
	Htmp = Qobj(H)
	r = rhotmp.eigenstates(sort='high')
	e = Htmp.eigenstates()
	r = r[0]
	e = e[0]
	s = sum(sum((rho*H))) - sum(r*e)
	return s