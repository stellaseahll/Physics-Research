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
        k = 10**(-(m+1.0)/10.0)  

        # Interaction Strength 
        g = 10**(-(mm+1.0)/10.0)


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
        H0 = omegaH*tensor(sigmap*sigmam,qeye(2),qeye(2))+omegaC*tensor(qeye(2),sigmap*sigmam,qeye(2))+omegaR*tensor(qeye(2),qeye(2),sigmap*sigmam)
        H = H0 + Hint

        # Coefficients for different dissipators
        nX1 = 1.0/(exp(abs(omegaH-g)*1.0/TH)-1.0)
        nX2 = 1.0/(exp(omegaH*1.0/TH)-1.0)
        nX3 = 1.0/(exp(abs(omegaH+g)*1.0/TH)-1.0) 
        nY1 = 1.0/(exp(abs(omegaC-g)*1.0/TC)-1.0)
        nY2 = 1.0/(exp(omegaC*1.0/TC)-1.0)
        nY3 = 1.0/(exp(abs(omegaC+g)*1.0/TC)-1.0) 
        nZ1 = 1.0/(exp(abs(omegaR-g)*1.0/TR)-1.0)
        nZ2 = 1.0/(exp(omegaR*1.0/TR)-1.0)
        nZ3 = 1.0/(exp(abs(omegaR+g)*1.0/TR)-1.0) 

        # # Dissipators 
        # if (omega>g):
        #     x1 = 0.5*(tensor(qeye(2),sigmap)+tensor(sigmap,sigmaz))
        #     x2 = 0.5*(tensor(qeye(2),sigmap)-tensor(sigmap,sigmaz))
        # else:
        #     x1 = 0.5*(tensor(qeye(2),sigmam)+tensor(sigmam,sigmaz))
        #     x2 = 0.5*(tensor(qeye(2),sigmap)-tensor(sigmap,sigmaz)) 


        # Lindblad Dissipators
        X1 = (ketminus*ket000.dag()+ket111*ketplus.dag())/sqrt(2)
        X2 = ket101*ket001.dag() + ket110*ket010.dag()
        X3 = (ketplus*ket000.dag()- ket111*ketminus.dag())/sqrt(2)
        Y1 = (ket110*ketplus.dag()-ketminus*ket001.dag())/sqrt(2)
        Y2 = ket010*ket000.dag() + ket111*ket101.dag()
        Y3 = (ketplus*ket001.dag()+ ket110*ketminus.dag())/sqrt(2)
        Z1 = (ket101*ketplus.dag()-ketminus*ket010.dag())/sqrt(2)
        Z2 = ket001*ket000.dag() + ket111*ket110.dag()
        Z3 = (ketplus*ket010.dag()+ ket101*ketminus.dag())/sqrt(2)

        CXM3= sqrt(k)*sqrt(nX3+1) * X3.dag()
        CXM2= sqrt(k)*sqrt(nX2+1) * X2.dag()
        CXM1= sqrt(k)*sqrt(nX1+1) * X1.dag()
        CX1 = sqrt(k)*sqrt(nX1) * X1
        CX2 = sqrt(k)*sqrt(nX2) * X2
        CX3 = sqrt(k)*sqrt(nX3) * X3

        CYM3= sqrt(k)*sqrt(nY3+1) * Y3.dag()
        CYM2= sqrt(k)*sqrt(nY2+1) * Y2.dag()
        CYM1= sqrt(k)*sqrt(nY1+1) * Y1.dag()
        CY1 = sqrt(k)*sqrt(nY1) * Y1
        CY2 = sqrt(k)*sqrt(nY2) * Y2
        CY3 = sqrt(k)*sqrt(nY3) * Y3

        CZM3= sqrt(k)*sqrt(nZ3+1) * Z3.dag()
        CZM2= sqrt(k)*sqrt(nZ2+1) * Z2.dag()
        CZM1= sqrt(k)*sqrt(nZ1+1) * Z1.dag()
        CZ1 = sqrt(k)*sqrt(nZ1) * Z1
        CZ2 = sqrt(k)*sqrt(nZ2) * Z2
        CZ3 = sqrt(k)*sqrt(nZ3) * Z3

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

        QH = CXM3.dag()*H*CXM3 - 0.5*H*CXM3.dag()*CXM3 - 0.5*CXM3.dag()*CXM3*H + CXM2.dag()*H*CXM2 - 0.5*H*CXM2.dag()*CXM2 - 0.5*CXM2.dag()*CXM2*H + CXM1.dag()*H*CXM1 - 0.5*H*CXM1.dag()*CXM1 - 0.5*CXM1.dag()*CXM1*H + CX3.dag()*H*CX3 - 0.5*H*CX3.dag()*CX3 - 0.5*CX3.dag()*CX3*H + CX2.dag()*H*CX2 - 0.5*H*CX2.dag()*CX2 - 0.5*CX2.dag()*CX2*H + CX1.dag()*H*CX1 - 0.5*H*CX1.dag()*CX1 - 0.5*CX1.dag()*CX1*H
        QC = CYM3.dag()*H*CYM3 - 0.5*H*CYM3.dag()*CYM3 - 0.5*CYM3.dag()*CYM3*H + CYM2.dag()*H*CYM2 - 0.5*H*CYM2.dag()*CYM2 - 0.5*CYM2.dag()*CYM2*H + CYM1.dag()*H*CYM1 - 0.5*H*CYM1.dag()*CYM1 - 0.5*CYM1.dag()*CYM1*H + CY3.dag()*H*CY3 - 0.5*H*CY3.dag()*CY3 - 0.5*CY3.dag()*CY3*H + CY2.dag()*H*CY2 - 0.5*H*CY2.dag()*CY2 - 0.5*CY2.dag()*CY2*H + CY1.dag()*H*CY1 - 0.5*H*CY1.dag()*CY1 - 0.5*CY1.dag()*CY1*H
        QR = CZM3.dag()*H*CZM3 - 0.5*H*CZM3.dag()*CZM3 - 0.5*CZM3.dag()*CZM3*H + CZM2.dag()*H*CZM2 - 0.5*H*CZM2.dag()*CZM2 - 0.5*CZM2.dag()*CZM2*H + CZM1.dag()*H*CZM1 - 0.5*H*CZM1.dag()*CZM1 - 0.5*CZM1.dag()*CZM1*H + CZ3.dag()*H*CZ3 - 0.5*H*CZ3.dag()*CZ3 - 0.5*CZ3.dag()*CZ3*H + CZ2.dag()*H*CZ2 - 0.5*H*CZ2.dag()*CZ2 - 0.5*CZ2.dag()*CZ2*H + CZ1.dag()*H*CZ1 - 0.5*H*CZ1.dag()*CZ1 - 0.5*CZ1.dag()*CZ1*H

        #---------------EXPECTATION VALUES----------------
        solution=mesolve(H, rho, timesteps, [CXM3,CXM2,CXM1,CX1,CX2,CX3,CYM3,CYM2,CYM1,CY1,CY2,CY3,CZM3,CZM2,CZM1,CZ1,CZ2,CZ3], [E1,E2,E3,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,QH,QC,QR],options=options) #
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
        file_data_store('Global_t%s_kappa%.3f_omegaH%.3f_omegaC%.3f_TH%.3f_TC%.3f_TR%.3f_g%.3f.dat'%(t,k,omegaH,omegaC,TH,TC,TR,g),E,numtype="real")
 
file_data_store('Global_Entropy_omegaH%.3f_TH%.3f_omegaC%.3f_TC%.3f_TR%.3f_g(0.001,0)_kappa(0.001,0).dat'%(args.omegaH_value,args.TH_value,args.omegaC_value,args.TC_value,args.TR_value),entropy,numtype="real")
file_data_store('Global_HeatFlow_omegaH%.3f_TH%.3f_omegaC%.3f_TC%.3f_TR%.3f_g(0.001,0)_kappa(0.001,0).dat'%(args.omegaH_value,args.TH_value,args.omegaC_value,args.TC_value,args.TR_value),heatflow,numtype="real")


