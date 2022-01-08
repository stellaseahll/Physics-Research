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
import numpy as np
import matplotlib.pyplot as plt


# plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')

plt.figure(figsize=(8,6),frameon=False)

# plt.gcf().subplots_adjust(left=0.12)
# plt.gcf().subplots_adjust(right=0.95)
# plt.gcf().subplots_adjust(top=0.95)
# plt.gcf().subplots_adjust(bottom=0.12)
# x = np.loadtxt("ergotropy.dat");
# line0 = plt.plot(x[0]*(sqrt(1/10.0)),x[1]/(sqrt(1/10.0)),linewidth=4.5,color='k');
# x = np.loadtxt("qubitho_hospread_data.dat");
# line1 = plt.fill(x[0]/10.0,x[1], color='#ffb1b1', lw=0, alpha=0.75,edgecolor='none');
# x = np.loadtxt("qubitho_qubitspread_data.dat");
# line2 = plt.fill(x[0]/10.0,x[1], color='#cae6ff', lw=0, alpha=0.75,edgecolor='none');
# x = np.loadtxt("qubitho_qubitmean_data.dat");
# line3 = plt.plot(x[0]/10.0,x[1],linewidth=2,color='#4c8bff');

# x = np.loadtxt("qubitho_homean_data.dat");
# line3 = plt.plot(x[0]/10.0,x[1],linewidth=2,linestyle=':',color='#ff674c');

x = np.loadtxt("hotreg.dat");
line1 = plt.fill(x[0]-1,x[1]-1, color='#ffb3ba', lw=0, edgecolor='none');
x = np.loadtxt("coldreg.dat");
line2 = plt.fill(x[0]-1,x[1]-1, color='#bae1ff', lw=0, edgecolor='none');
x = np.loadtxt("THdata0.01.dat");
line3 = plt.plot(x[0],x[1], color='#1a3547', linewidth=2);
line4 = plt.plot(x[0],x[2], color='#1a3547', linewidth=2);
line5 = plt.plot(x[0],x[3], color='#1a3547', linewidth=2);
line6 = plt.plot(x[0],x[4], color='#1a3547', linewidth=2);
line7 = plt.plot(x[0],x[5], color='#1a3547', linewidth=2);
line8 = plt.plot(x[0],x[6], color='#1a3547', linewidth=2);
line9 = plt.plot(x[0],x[7], color='#1a3547', linewidth=2);
line8 = plt.plot(x[0],x[8], color='#1a3547', linewidth=2);
line9 = plt.plot(x[0],x[9], color='#1a3547', linewidth=2);
# line6 = plt.plot(x[0],x[4], color='#1a3547', linewidth=2);

# x = np.loadtxt("qubitho_qubitmean_data.dat");
# line3 = plt.plot(x[0]/10.0,x[1],linewidth=2,color='#1c2d83');

# x = np.loadtxt("qubitho_homean_data.dat");
# line3 = plt.plot(x[0]/10.0,x[1],linewidth=3.5,linestyle=':',color='#bd0687');
plt.xlabel(r'$\mathrm{T_w-T_c}$',fontsize=18)
plt.ylabel(r'$\mathrm{T_h-T_c}$',fontsize=18)
plt.gca().autoscale(enable=True, axis='both', tight=None)


# plt.grid(color='0.7', linestyle='--', linewidth=1.0)
# plt.gca().spines['top'].set_visible(False)
# plt.gca().spines['right'].set_visible(False)
plt.tick_params(top='off', bottom='on', left='on', right='off', labelleft='on', labelbottom='on',direction='in',length=6,width=1)
plt.tick_params(axis='both', which='major', labelsize=20)
ax1 = plt.subplot()
ax1.yaxis.set_label_coords(-0.07,0.5)
# ax.spines['bottom'].set_visible(False)
# ax.spines['left'].set_visible(False)
# # plt.legend()
# plt.axvline(x=0,linewidth=1.5, color='k')
# plt.axhline(y=-5,linewidth=1.5, color='k')
plt.axis([0,20, 0 ,3])
plt.yticks([0,1,2,3],fontsize=18)
plt.xticks([0,5,10,15,20],fontsize=18)
# plt.legend(frameon='False',edgecolor='w')
plt.savefig("test.pdf",dpi=400)


# x= np.loadtxt("test.dat",dtype=float)
# I = 10.0;
# dt = x[0][1]-x[0][0];
# Wrate0 = np.zeros(20000);
# Wrate1 = np.zeros(20000);
# t = np.zeros(20000);
# for j in range(0,20000):
#     Wrate0[j] = (x[3][j+1]-x[3][j])*x[3][j]/dt/I/(sqrt(1/10.0));
#     Wrate1[j] = (x[4][j+1]-x[4][j])/2/I/dt/(sqrt(1/10.0));
#     t[j] = x[0][j]*sqrt(1/10.0)
# plt.close()
# line1 = plt.plot(t,Wrate0,linewidth=3.5,color='#ECA845');
# line2 = plt.plot(x[0]*(sqrt(1/10.0)),x[8]/(sqrt(1/10.0)),'--',linewidth=3.5,color='k');
# line3 = plt.plot(t,Wrate1,':',linewidth=4.5,color='#52812C');
# # plt.legend()
# plt.grid(color='0.8', linestyle=':', linewidth=1.5)
# # plt.axis([20,100,0.2,1.4])
# plt.xticks(fontsize=12,stretch='extra-condensed',weight='bold')
# plt.yticks(fontsize=12,stretch='extra-condensed',weight='bold')
# # plt.legend(frameon='False',edgecolor='w')
# plt.savefig("test.png",dpi=1000)

