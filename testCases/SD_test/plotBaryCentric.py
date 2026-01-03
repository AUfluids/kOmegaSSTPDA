import numpy as np
import pandas as pd
import os
import Ofpp
print(os.getcwd())
import matplotlib.pyplot as plt

import matplotlib as mpl
from matplotlib.lines import Line2D
from scipy.interpolate import griddata



def baryCentricCalculation (Rij):
    # calculation Of baryCentric Points with Rij

    NumData = Rij.shape[0]

    tke = 0.5*(Rij[:,0] + Rij[:,3] + Rij[:,5])
    print(tke)
    zeta = np.zeros(NumData)
    eta = np.zeros(NumData)
    for i in range(NumData):
        b = np.zeros((3,3))
        b[0,0] =  Rij[i,0]/tke[i]/2-1/3
        b[0,1] =  Rij[i,1]/tke[i]/2
        b[0,2] =  Rij[i,2]/tke[i]/2
        b[1,0] =  Rij[i,1]/tke[i]/2
        b[1,1] =  Rij[i,3]/tke[i]/2-1/3
        b[1,2] =  Rij[i,4]/tke[i]/2
        b[2,0] =  Rij[i,2]/tke[i]/2
        b[2,1] =  Rij[i,4]/tke[i]/2
        b[2,2] =  Rij[i,5]/tke[i]/2-1/3
        w, v = np.linalg.eig(b) # eigen value calculation
        # sorting the eigen values
        idx = (-w).argsort()
        w = w[idx]
        v = v[:,idx]
        lamda1 = w[0]
        lamda2 = w[1]
        lamda3 = w[2]
        c1 = lamda1 - lamda2
        c2 = 2*(lamda2-lamda3)
        c3 = 3*lamda3 + 1
        zeta[i] = c1 - c2
        eta[i] = np.sqrt(3.0)*c3
    return eta, zeta



csvData = []
plt.rcParams.update({'font.family': 'serif',
                'font.weight': 'normal',
                'font.size': 24,
                })

fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.subplots_adjust(hspace=0.3)

fig.set_figwidth(21)
fig.set_figheight(7)


# plot frames

ax1.plot(np.array([-1,1,0,-1]), np.array([0,0,3**0.5,0]), 'k-')
ax1.set_title("KOSST")
ax2.plot(np.array([-1,1,0,-1]), np.array([0,0,3**0.5,0]), 'k-')
ax2.set_title("This Model")
ax3.plot(np.array([-1,1,0,-1]), np.array([0,0,3**0.5,0]), 'k-')
ax3.set_title("High Fidelity")

# KOSST

fname = './!RefData/Rij_KOSST'
Rij_KOSST = Ofpp.parse_internal_field(fname)
eta, zeta = baryCentricCalculation (Rij_KOSST)
ax1.plot(zeta, eta, 'g*')


# High Fidelity Data

fname = './!RefData/Aij_HF'
Aij_HF = Ofpp.parse_internal_field(fname)

fname = './!RefData/k_HF'
k_HF = Ofpp.parse_internal_field(fname)

temp = np.array([[1,0,0,1,0,1],])*np.reshape(k_HF,[k_HF.shape[0],1])

Rij_HF = Aij_HF + 2/3*np.array([[1,0,0,1,0,1],])*np.reshape(k_HF,[k_HF.shape[0],1])
eta, zeta = baryCentricCalculation (Rij_HF)
ax3.plot(zeta, eta, 'r*')


# This case

list_subfolders = [f.name for f in os.scandir(os.getcwd()) if f.is_dir()]
list_subfoders_num = []
for folder in list_subfolders:
    if folder != '0':
        if (folder.isnumeric()):
            list_subfoders_num.append(folder)
list_subfoders_num.sort()
res = [eval(i) for i in list_subfoders_num]
selectedFolder = str(np.max(res))


fname = os.getcwd()+ '/' + selectedFolder +  '/turbulenceProperties:R'
RijLinear_internal = Ofpp.parse_internal_field(fname)

fname = os.getcwd()+ '/' + selectedFolder +  '/dBij'
dBij_internal = Ofpp.parse_internal_field(fname)

fname = os.getcwd()+ '/' + selectedFolder +  '/k'
tke_internal = Ofpp.parse_internal_field(fname)

Rij_ = RijLinear_internal + dBij_internal*np.reshape(tke_internal,[tke_internal.shape[0],1])

eta, zeta = baryCentricCalculation (Rij_)
ax2.plot(zeta, eta, 'b*')


nameTosave = 'baryCentric' + selectedFolder + '.png'
plt.savefig(nameTosave)

