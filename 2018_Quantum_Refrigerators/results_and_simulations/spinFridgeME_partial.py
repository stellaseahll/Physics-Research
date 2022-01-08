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
parser.add_argument("-TH", "--TH", dest="TH_value", type=float, metavar='<float>', default=2.0, help="TH value (default=2.0)")
parser.add_argument("-TC", "--TC", dest="TC_value", type=float, metavar='<float>', default=1.0, help="TC value (default=1.0)")
parser.add_argument("-TR", "--TR", dest="TR_value", type=float, metavar='<float>', default=5.0, help="TR value (default=5.0)")
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
entropy = np.zeros((29,29))
heatflow = np.zeros((29,29))
for m in range(29):
    for mm in range(29):
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
        H0 = 0.5*omegaH*tensor(sigmaz,qeye(2),qeye(2))+0.5*omegaC*tensor(qeye(2),sigmaz,qeye(2))+0.5*omegaR*tensor(qeye(2),qeye(2),sigmaz)
        H = H0 + Hint

        # Coefficients for different dissipators
        nX1 = 1.0/(exp(abs(omegaH-g)*1.0/TH)-1.0)
        nX2 = 1.0/(exp(omegaH*1.0/TH)-1.0)
        nX3 = 1.0/(exp(abs(omegaH+g)*1.0/TH)-1.0) 
        nX4 = 1.0/(exp(abs(omegaH+0.5*g)*1.0/TH)-1.0) 
        nX5 = 1.0/(exp(abs(omegaH-0.5*g)*1.0/TH)-1.0) 
        nY1 = 1.0/(exp(abs(omegaC-g)*1.0/TC)-1.0)
        nY2 = 1.0/(exp(omegaC*1.0/TC)-1.0)
        nY3 = 1.0/(exp(abs(omegaC+g)*1.0/TC)-1.0)
        nY4 = 1.0/(exp(abs(omegaC+0.5*g)*1.0/TC)-1.0) 
        nY5 = 1.0/(exp(abs(omegaC-0.5*g)*1.0/TC)-1.0)  
        nZ1 = 1.0/(exp(abs(omegaR-g)*1.0/TR)-1.0)
        nZ2 = 1.0/(exp(omegaR*1.0/TR)-1.0)
        nZ3 = 1.0/(exp(abs(omegaR+g)*1.0/TR)-1.0) 
        nZ4 = 1.0/(exp(abs(omegaR+0.5*g)*1.0/TR)-1.0) 
        nZ5 = 1.0/(exp(abs(omegaR-0.5*g)*1.0/TR)-1.0) 

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

        #TO DO: Write function to diagonalise the dissipators.

        # Diagonalise Dissipators of Hot Mode
        A = np.zeros((3,3))
        A[0][0] = nX3+1
        A[1][1] = nX2+1
        A[2][2] = nX1+1
        A[0][1] = nX4+1
        A[1][0] = nX4+1
        A[0][2] = nX2+1
        A[2][0] = nX2+1
        A[1][2] = nX5+1
        A[2][1] = nX5+1

        r = np.linalg.eig(A)
        eigval = r[0]
        eigvec = r[1]
        if (eigval[0]<0):
           eigval[0]=0
        if (eigval[1]<0):
           eigval[1]=0
        if (eigval[2]<0):
           eigval[2]=0
        LX1 = sqrt(eigval[0]*k)*(eigvec[0][0]*X3.dag() + eigvec[1][0]*X2.dag() + eigvec[2][0]*X1.dag())
        LX2 = sqrt(eigval[1]*k)*(eigvec[0][1]*X3.dag() + eigvec[1][1]*X2.dag() + eigvec[2][1]*X1.dag())
        LX3 = sqrt(eigval[2]*k)*(eigvec[0][2]*X3.dag() + eigvec[1][2]*X2.dag() + eigvec[2][2]*X1.dag())

        A = np.zeros((3,3))
        A[0][0] = nX1
        A[1][1] = nX2
        A[2][2] = nX3
        A[0][1] = nX5
        A[1][0] = nX5
        A[0][2] = nX2
        A[2][0] = nX2
        A[1][2] = nX4
        A[2][1] = nX4

        r = np.linalg.eig(A)
        eigval = r[0]
        eigvec = r[1]
        if (eigval[0]<0):
           eigval[0]=0
        if (eigval[1]<0):
           eigval[1]=0
        if (eigval[2]<0):
           eigval[2]=0
        LX4 = sqrt(eigval[0]*k)*(eigvec[0][0]*X1 + eigvec[1][0]*X2 + eigvec[2][0]*X3)
        LX5 = sqrt(eigval[1]*k)*(eigvec[0][1]*X1 + eigvec[1][1]*X2 + eigvec[2][1]*X3)
        LX6 = sqrt(eigval[2]*k)*(eigvec[0][2]*X1 + eigvec[1][2]*X2 + eigvec[2][2]*X3)

        # Diagonalise Dissipators of Cold Mode
        A = np.zeros((3,3))
        A[0][0] = nY3+1
        A[1][1] = nY2+1
        A[2][2] = nY1+1
        A[0][1] = nY4+1
        A[1][0] = nY4+1
        A[0][2] = nY2+1
        A[2][0] = nY2+1
        A[1][2] = nY5+1
        A[2][1] = nY5+1

        r = np.linalg.eig(A)
        eigval = r[0]
        eigvec = r[1]
        if (eigval[0]<0):
           eigval[0]=0
        if (eigval[1]<0):
           eigval[1]=0
        if (eigval[2]<0):
           eigval[2]=0
        LY1 = sqrt(eigval[0]*k)*(eigvec[0][0]*Y3.dag() + eigvec[1][0]*Y2.dag() + eigvec[2][0]*Y1.dag())
        LY2 = sqrt(eigval[1]*k)*(eigvec[0][1]*Y3.dag() + eigvec[1][1]*Y2.dag() + eigvec[2][1]*Y1.dag())
        LY3 = sqrt(eigval[2]*k)*(eigvec[0][2]*Y3.dag() + eigvec[1][2]*Y2.dag() + eigvec[2][2]*Y1.dag())

        A = np.zeros((3,3))
        A[0][0] = nY1
        A[1][1] = nY2
        A[2][2] = nY3
        A[0][1] = nY5
        A[1][0] = nY5
        A[0][2] = nY2
        A[2][0] = nY2
        A[1][2] = nY4
        A[2][1] = nY4

        r = np.linalg.eig(A)
        eigval = r[0]
        eigvec = r[1]
        if (eigval[0]<0):
           eigval[0]=0
        if (eigval[1]<0):
           eigval[1]=0
        if (eigval[2]<0):
           eigval[2]=0
        LY4 = sqrt(eigval[0]*k)*(eigvec[0][0]*Y1 + eigvec[1][0]*Y2 + eigvec[2][0]*Y3)
        LY5 = sqrt(eigval[1]*k)*(eigvec[0][1]*Y1 + eigvec[1][1]*Y2 + eigvec[2][1]*Y3)
        LY6 = sqrt(eigval[2]*k)*(eigvec[0][2]*Y1 + eigvec[1][2]*Y2 + eigvec[2][2]*Y3)

        # Diagonalise Dissipators of Work Mode
        A = np.zeros((3,3))
        A[0][0] = nZ3+1
        A[1][1] = nZ2+1
        A[2][2] = nZ1+1
        A[0][1] = nZ4+1
        A[1][0] = nZ4+1
        A[0][2] = nZ2+1
        A[2][0] = nZ2+1
        A[1][2] = nZ5+1
        A[2][1] = nZ5+1

        r = np.linalg.eig(A)
        eigval = r[0]
        eigvec = r[1]
        if (eigval[0]<0):
           eigval[0]=0
        if (eigval[1]<0):
           eigval[1]=0
        if (eigval[2]<0):
           eigval[2]=0
        LZ1 = sqrt(eigval[0]*k)*(eigvec[0][0]*Z3.dag() + eigvec[1][0]*Z2.dag() + eigvec[2][0]*Z1.dag())
        LZ2 = sqrt(eigval[1]*k)*(eigvec[0][1]*Z3.dag() + eigvec[1][1]*Z2.dag() + eigvec[2][1]*Z1.dag())
        LZ3 = sqrt(eigval[2]*k)*(eigvec[0][2]*Z3.dag() + eigvec[1][2]*Z2.dag() + eigvec[2][2]*Z1.dag())

        A = np.zeros((3,3))
        A[0][0] = nZ1
        A[1][1] = nZ2
        A[2][2] = nZ3
        A[0][1] = nZ5
        A[1][0] = nZ5
        A[0][2] = nZ2
        A[2][0] = nZ2
        A[1][2] = nZ4
        A[2][1] = nZ4

        r = np.linalg.eig(A)
        eigval = r[0]
        eigvec = r[1]
        if (eigval[0]<0):
           eigval[0]=0
        if (eigval[1]<0):
           eigval[1]=0
        if (eigval[2]<0):
           eigval[2]=0
        LZ4 = sqrt(eigval[0]*k)*(eigvec[0][0]*Z1 + eigvec[1][0]*Z2 + eigvec[2][0]*Z3)
        LZ5 = sqrt(eigval[1]*k)*(eigvec[0][1]*Z1 + eigvec[1][1]*Z2 + eigvec[2][1]*Z3)
        LZ6 = sqrt(eigval[2]*k)*(eigvec[0][2]*Z1 + eigvec[1][2]*Z2 + eigvec[2][2]*Z3)

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

        QH = LX6.dag()*H*LX6 - 0.5*H*LX6.dag()*LX6 - 0.5*LX6.dag()*LX6*H + LX5.dag()*H*LX5 - 0.5*H*LX5.dag()*LX5 - 0.5*LX5.dag()*LX5*H + LX4.dag()*H*LX4 - 0.5*H*LX4.dag()*LX4 - 0.5*LX4.dag()*LX4*H + LX3.dag()*H*LX3 - 0.5*H*LX3.dag()*LX3 - 0.5*LX3.dag()*LX3*H + LX2.dag()*H*LX2 - 0.5*H*LX2.dag()*LX2 - 0.5*LX2.dag()*LX2*H + LX1.dag()*H*LX1 - 0.5*H*LX1.dag()*LX1 - 0.5*LX1.dag()*LX1*H
        QC = LY6.dag()*H*LY6 - 0.5*H*LY6.dag()*LY6 - 0.5*LY6.dag()*LY6*H + LY5.dag()*H*LY5 - 0.5*H*LY5.dag()*LY5 - 0.5*LY5.dag()*LY5*H + LY4.dag()*H*LY4 - 0.5*H*LY4.dag()*LY4 - 0.5*LY4.dag()*LY4*H + LY3.dag()*H*LY3 - 0.5*H*LY3.dag()*LY3 - 0.5*LY3.dag()*LY3*H + LY2.dag()*H*LY2 - 0.5*H*LY2.dag()*LY2 - 0.5*LY2.dag()*LY2*H + LY1.dag()*H*LY1 - 0.5*H*LY1.dag()*LY1 - 0.5*LY1.dag()*LY1*H
        QR = LZ6.dag()*H*LZ6 - 0.5*H*LZ6.dag()*LZ6 - 0.5*LZ6.dag()*LZ6*H + LZ5.dag()*H*LZ5 - 0.5*H*LZ5.dag()*LZ5 - 0.5*LZ5.dag()*LZ5*H + LZ4.dag()*H*LZ4 - 0.5*H*LZ4.dag()*LZ4 - 0.5*LZ4.dag()*LZ4*H + LZ3.dag()*H*LZ3 - 0.5*H*LZ3.dag()*LZ3 - 0.5*LZ3.dag()*LZ3*H + LZ2.dag()*H*LZ2 - 0.5*H*LZ2.dag()*LZ2 - 0.5*LZ2.dag()*LZ2*H + LZ1.dag()*H*LZ1 - 0.5*H*LZ1.dag()*LZ1 - 0.5*LZ1.dag()*LZ1*H


        #---------------EXPECTATION VALUES----------------
        solution=mesolve(H, rho, timesteps, [LX1,LX2,LX3,LX4,LX5,LX6,LY1,LY2,LY3,LY4,LY5,LY6,LZ1,LZ2,LZ3,LZ4,LZ5,LZ6],  [E1,E2,E3,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,QH,QC,QR],options=options)
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
        E=np.vstack((timesteps,a1,a2,a3,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,QH,QC,QR,S,QH/TH,QC/TC,QR/TR))
        entropy[m][mm] = -QH[NT-1]/TH-QC[NT-1]/TC-QR[NT-1]/TR
        heatflow[m][mm] = QC[NT-1]
        file_data_store('Partial_t%s_kappa%.3f_omegaH%.3f_omegaC%.3f_TH%.3f_TC%.3f_TR%.3f_g%.3f.dat'%(t,k,omegaH,omegaC,TH,TC,TR,g),E,numtype="real")

file_data_store('Partial_Entropy_omegaH%.3f_TH%.3f_omegaC%.3f_TC%.3f_TR%.3f_g(0.001,0)_kappa(0.001,0).dat'%(args.omegaH_value,args.TH_value,args.omegaC_value,args.TC_value,args.TR_value),entropy,numtype="real")
file_data_store('Partial_HeatFlow_omegaH%.3f_TH%.3f_omegaC%.3f_TC%.3f_TR%.3f_g(0.001,0)_kappa(0.001,0).dat'%(args.omegaH_value,args.TH_value,args.omegaC_value,args.TC_value,args.TR_value),heatflow,numtype="real")


