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
from collections import OrderedDict
from matplotlib.transforms import blended_transform_factory
import matplotlib.ticker as ticker
from matplotlib.ticker import FuncFormatter
# linestyles = OrderedDict(
#     [('solid',               (0, ())),
#      ('loosely dotted',      (0, (1, 10))),
#      ('dotted',              (0, (1, 5))),
#      ('densely dotted',      (0, (1, 1))),

#      ('loosely dashed',      (0, (5, 10))),
#      ('dashed',              (0, (5, 5))),
#      ('densely dashed',      (0, (5, 1))),

#      ('loosely dashdotted',  (0, (3, 10, 1, 10))),
#      ('dashdotted',          (0, (3, 5, 1, 5))),
#      ('densely dashdotted',  (0, (3, 1, 1, 1))),

#      ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
#      ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
#      ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])




# def millions(x, pos):
#     return '%1.6f' % (x*1e-4)

# formatter = FuncFormatter(millions)
fig, ax = plt.subplots()
# plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.figure(figsize=(8,6),frameon=False)
plt.grid(color='0.9', alpha=0.8,linestyle='--', linewidth=1.0)
# ax.yaxis.set_major_formatter(formatter)
plt.gcf().subplots_adjust(left=0.14)
plt.gcf().subplots_adjust(right=0.95)
plt.gcf().subplots_adjust(top=0.89)
plt.gcf().subplots_adjust(bottom=0.12)





x = np.loadtxt("compareModels.dat");
# line1 = plt.semilogx(x[0],x[1], color='#000000', linestyle='-',linewidth=2);#,marker='o',markevery=12,markersize=5.5); #partialxxx
# line4 = plt.semilogx(x[0],x[4], color='#5d5d5d', linestyle='-',linewidth=2);#,marker='s',markevery=12,markersize=5.5); #partialres
# line2 = plt.semilogx(x[0],x[2], color='#0000e6', linestyle=(0, (1, 1)),linewidth=2); #localxxx
# line3 = plt.semilogx(x[0],x[3], color='#e30074', linestyle=(0, (1, 1)),linewidth=2); #globalxxx
# line5 = plt.semilogx(x[0],x[5], color='#1d76fc', linestyle='--',linewidth=2); #localres
# line6 = plt.semilogx(x[0],x[6], color='#ff80a0', linestyle='--',linewidth=2); #globalres
line1 = plt.semilogx(x[0],x[1]/0.01, color='#000000', linestyle='-',linewidth=2);#,marker='o',markevery=12,markersize=5.5); #partialxxx
line4 = plt.semilogx(x[0],x[4]/0.01, color='#e30074', linestyle='-',linewidth=2);#,marker='s',markevery=12,markersize=5.5); #partialres
line2 = plt.semilogx(x[0],x[2]/0.01, color='#000000', linestyle=(0, (1, 1)),linewidth=3); #localxxx
line3 = plt.semilogx(x[0],x[3]/0.01, color='#000000', linestyle='--',linewidth=3); #globalxxx
line5 = plt.semilogx(x[0],x[5]/0.01, color='#e30074', linestyle=(0, (1, 1)),linewidth=3); #localres
line6 = plt.semilogx(x[0],x[6]/0.01, color='#e30074', linestyle='--',linewidth=3); #globalres
# line1 = plt.semilogx(x[0],x[1], color='#000000', linestyle='-',linewidth=0,marker='o',markevery=12,markersize=5.5,alpha = 0.2); #partialxxx
# line4 = plt.semilogx(x[0],x[4], color='#5d5d5d', linestyle='-',linewidth=0,marker='s',markevery=12,markersize=5.5,alpha = 0.2); #partialres

# line6 = plt.semilogx(x[0],x[4], color='#1a3547', linewidth=2);

# x = np.loadtxt("qubitho_qubitmean_data.dat");
# line3 = plt.semilogx(x[0]/10.0,x[1],linewidth=2,color='#1c2d83');

# x = np.loadtxt("qubitho_homean_data.dat");
# line3 = plt.semilogx(x[0]/10.0,x[1],linewidth=3.5,linestyle=':',color='#bd0687');
plt.xlabel(r'$\mathrm{Coupling\,Strength\,}\mathrm{g/\omega_c}$',fontsize=18)
plt.ylabel(r'$\mathrm{Heat\,Flow\,}\mathcal{P}_c/\hbar\omega_c\kappa_c$',fontsize=18)
plt.gca().autoscale(enable=True, axis='both', tight=None)
# ax.yaxis.set_major_formatter(ticker.ScalarFormatter("%.10f",useMathText=True))
# ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.3f"))

# plt.grid(color='0.7', linestyle='--', linewidth=1.0)
# plt.gca().spines['top'].set_visible(False)
# plt.gca().spines['right'].set_visible(False)
plt.tick_params(top='off', bottom='on', left='on', right='off', labelleft='on', labelbottom='on',direction='in',length=6,width=1)
plt.tick_params(axis='both', which='major', labelsize=18)
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_powerlimits((-3,2))
# ax.yaxis.set_major_formatter(formatter)
# y_labels = ax.get_yticks()
plt.gca().yaxis.set_major_formatter(formatter)
t = plt.gca().yaxis.get_offset_text()
t.set_size(16)
# ax1 = plt.subplot()
# ax1.yaxis.set_label_coords(-0.07,0.5)
# # ax.spines['bottom'].set_visible(False)
# # ax.spines['left'].set_visible(False)
# # # plt.legend()
# # plt.axvline(x=0,linewidth=1.5, color='k')
# # plt.axhline(y=-5,linewidth=1.5, color='k')
# plt.yticks([1,2,3,4],fontsize=12)
# plt.xticks([5,10,15,20],fontsize=12)
# plt.legend(frameon='False',edgecolor='w')
plt.savefig("compareModels.pdf",dpi=400)


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
# line1 = plt.semilogx(t,Wrate0,linewidth=3.5,color='#ECA845');
# line2 = plt.semilogx(x[0]*(sqrt(1/10.0)),x[8]/(sqrt(1/10.0)),'--',linewidth=3.5,color='k');
# line3 = plt.semilogx(t,Wrate1,':',linewidth=4.5,color='#52812C');
# # plt.legend()
# plt.grid(color='0.8', linestyle=':', linewidth=1.5)
# # plt.axis([20,100,0.2,1.4])
# plt.xticks(fontsize=12,stretch='extra-condensed',weight='bold')
# plt.yticks(fontsize=12,stretch='extra-condensed',weight='bold')
# # plt.legend(frameon='False',edgecolor='w')
# plt.savefig("test.png",dpi=1000)

