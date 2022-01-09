from qutip import *

def extractWork(rho,H):
	rhotmp = Qobj(rho)
	Htmp = Qobj(H)
	r = rhotmp.eigenstates(sort='high')
	e = Htmp.eigenstates()
	r = r[0]
	e = e[0]
	s = sum(sum((rho*H))) - sum(r*e)
	return s


	# r = rho.eigenstates()
	# e = H.eigenstates()
	# d = len(r[0])
	# s = 0
	# for j in range(0,d-1):
	# 	for k in range(0,d-1):
	# 		m = (d-1)-j
	# 		if (j==k):
	# 			s = s + r[0][m]*e[0][k] * (((e[1][k].dag()*r[1][m]).norm())**2 - 1.0)
	# 		else:
	# 			s = s + r[0][m]*e[0][k] * (((e[1][k].dag()*r[1][m]).norm())**2)
	# return s
