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
parser.add_argument("-g", "--g", dest="g_value", type=float, metavar='<float>', default=-10, help="G value (default=1.0)")
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
entropy = np.zeros((30,30))
heatflow = np.zeros((30,30))
for m in range(30):
    for mm in range(30):
# Energy of the qubits 
        omegaH = 2 + m/2.0
        omegaC = args.omegaC_value
        omegaR = omegaH - omegaC

        
        TC = args.TC_value
        TH = args.TH_value
        TR = TH + mm/2.0
        tauH = thermalState(omegaH,TH) #Qubit 1 starts with ground state
        tauC = thermalState(omegaC,TC) #Qubit 2 is in an excited state 
        tauR = thermalState(omegaR,TR)

        # Coupling strength
        k = args.k_value    

        # Interaction Strength 
        g = 10**(args.g_value/10.0)

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
        Hint = g*(tensor(sigmam,sigmap,sigmap)+tensor(sigmap,sigmam,sigmam))
        H0 = 0.5*omegaH*tensor(sigmaz,qeye(2),qeye(2))+0.5*omegaC*tensor(qeye(2),sigmaz,qeye(2))+0.5*omegaR*tensor(qeye(2),qeye(2),sigmaz)
        H = H0 + Hint

        # Coefficients for different dissipators
        nX = 1.0/(exp(omegaH*1.0/TH)-1.0)
        nY = 1.0/(exp(omegaC*1.0/TC)-1.0)
        nZ = 1.0/(exp(omegaR*1.0/TR)-1.0)
        # # Dissipators 
        # if (omega>g):
        #     x1 = 0.5*(tensor(qeye(2),sigmap)+tensor(sigmap,sigmaz))
        #     x2 = 0.5*(tensor(qeye(2),sigmap)-tensor(sigmap,sigmaz))
        # else:
        #     x1 = 0.5*(tensor(qeye(2),sigmam)+tensor(sigmam,sigmaz))
        #     x2 = 0.5*(tensor(qeye(2),sigmap)-tensor(sigmap,sigmaz)) 


        # Lindblad Dissipators
        X = tensor(sigmap,qeye(2),qeye(2))
        Y = tensor(qeye(2),sigmap,qeye(2))
        Z = tensor(qeye(2),qeye(2),sigmap)

        CXM = sqrt(k)*sqrt(nX+1)*X.dag() 
        CX = sqrt(k)*sqrt(nX)*X
        CYM = sqrt(k)*sqrt(nY+1)*Y.dag() 
        CY = sqrt(k)*sqrt(nY)*Y
        CZM = sqrt(k)*sqrt(nZ+1)*Z.dag() 
        CZ = sqrt(k)*sqrt(nZ)*Z

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

        QH = CXM.dag()*H*CXM - 0.5*H*CXM.dag()*CXM - 0.5*CXM.dag()*CXM*H + CX.dag()*H*CX - 0.5*H*CX.dag()*CX - 0.5*CX.dag()*CX*H 
        QC = CYM.dag()*H*CYM - 0.5*H*CYM.dag()*CYM - 0.5*CYM.dag()*CYM*H + CY.dag()*H*CY - 0.5*H*CY.dag()*CY - 0.5*CY.dag()*CY*H 
        QR = CZM.dag()*H*CZM - 0.5*H*CZM.dag()*CZM - 0.5*CZM.dag()*CZM*H + CZ.dag()*H*CZ - 0.5*H*CZ.dag()*CZ - 0.5*CZ.dag()*CZ*H 

        #---------------EXPECTATION VALUES----------------
        solution=mesolve(H, rho, timesteps, [CXM,CX,CYM,CY,CZM,CZ], [E1,E2,E3,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,QH,QC,QR],options=options)
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
        entropy[m][mm] = -QH[NT-1]/TH-QC[NT-1]/TC-QR[NT-1]/TR
        heatflow[m][mm] = QC[NT-1]
        E=np.vstack((timesteps,a1,a2,a3,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,QH,QC,QR,S,QH/TH,QC/TC,QR/TR))
        file_data_store('Local_t%s_kappa%.3f_omegaH%.3f_omegaC%.3f_TH%.3f_TC%.3f_TR%.3f_g%.3f.dat'%(t,k,omegaH,omegaC,TH,TC,TR,args.g_value/10.0),E,numtype="real")
        file_data_store('FinalStateLocal_t%s_kappa%.3f_omegaH%.3f_omegaC%.3f_TH%.3f_TC%.3f_TR%.3f_g%.3f.dat'%(t,k,omegaH,omegaC,TH,TC,TR,args.g_value/10.0),solution.states[NT-1],numtype="complex")

file_data_store('Local_Entropy_omegaH(2,17)_TH(5,20)_omegaC%.3f_TC%.3f_TR%.3f_g%.3f_kappa%.3f.dat'%(args.omegaC_value,args.TC_value,args.TR_value,args.g_value/10.0,args.k_value),entropy,numtype="real")
file_data_store('Local_HeatFlow_omegaH(2,17)_TH(5,20)_omegaC%.3f_TC%.3f_TR%.3f_g%.3f_kappa%.3f.dat'%(args.omegaC_value,args.TC_value,args.TR_value,args.g_value/10.0,args.k_value),heatflow,numtype="real")


