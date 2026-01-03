import numpy as np
import pandas as pd
import os
import Ofpp
print(os.getcwd())
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"],
    'font.size': 7,
    'axes.linewidth': 2,
    "grid.color": "#cccccc",
    "axes.grid": True,
    "grid.alpha": 0.8,
    "grid.linewidth": 0.5,
    "grid.linestyle": '-',
    "axes.grid.which": 'both',
    "axes.spines.right": False,
    "axes.spines.top": False,
    'axes.axisbelow': True})

cm = 1 / 2.54  # inches to cm
textWidth = 9  # cm
figSize = (textWidth * cm, textWidth * cm * 3 / 4)

figSizeMedium = (14 * cm, 14 * cm * 3 / 4)
figSizeFull = (19 * cm, 19 * cm)

Ny = 48
Nz = 48
coordinates = np.array([0.0124562000000000, 0.0368192000000000, 0.0601075000000000, 0.0823686000000000,0.103648000000000,0.123988000000000,0.143432000000000,0.162017000000000,0.179783000000000,0.196766000000000,0.212999000000000,0.228516000000000,0.243348000000000,0.257527000000000,0.271080000000000,0.284035000000000,0.296419000000000,0.308256000000000,0.319572000000000,0.330388000000000,0.340727000000000,0.350610000000000,0.360057000000000,0.369088000000000,0.377720000000000,0.385971000000000,0.393859000000000,0.401398000000000,0.408605000000000,0.415494000000000,0.422079000000000,0.428374000000000,0.434391000000000,0.440142000000000,0.445640000000000,0.450896000000000,0.455919000000000,0.460721000000000,0.465312000000000,0.469699000000000,0.473893000000000,0.477903000000000,0.481735000000000,0.485398000000000,0.488900000000000,0.492247000000000,0.495447000000000,0.498505000000000])
coordinates = (coordinates-coordinates[0])/(coordinates[-1]-coordinates[0]) * 0.5
# print(coordinates)
y = 2 * (coordinates-coordinates[0])
z = 2 * (0.5 - np.flip(coordinates))
# y = np.linspace(0,1,Ny)
# z = np.linspace(0,1,Nz)
# fRANS - NN(LES input)
folder = 'SavedPic/'

list_subfolders = [f.name for f in os.scandir(os.getcwd()) if f.is_dir()]
# print(list_subfolders)

list_subfoders_num = []
for folder in list_subfolders:
    if folder != '0':
        if (folder.isnumeric()):
            list_subfoders_num.append(folder)
list_subfoders_num.sort()
res = [eval(i) for i in list_subfoders_num]
selectedFolder = str(np.max(res))




plt.clf()
fname = os.getcwd()+ '/' + selectedFolder +  '/U'

U_internal = Ofpp.parse_internal_field(fname)
# U_internal = U_internal[:2304,:]
# np.savetxt('check',U_internal)
Ux_grid = np.reshape(U_internal[:,0],(Nz,Ny))
Uy_grid = np.reshape(U_internal[:,1],(Nz,Ny))
Uz_grid = np.reshape(U_internal[:,2],(Nz,Ny))
Ux_meanY = np.mean(Ux_grid,axis=1)
Ux_pp = Ux_grid - np.reshape(Ux_meanY,(Nz,1))
Ux_meanTop = np.mean(Ux_grid[:,-1])



################ EXTRA STUFF BY MARIO #################

# READ THETA VALUES #
dfTheta = pd.read_csv(os.getcwd() + '/constant/theta_values.csv')
dfTheta.drop(dfTheta.columns[0], axis=1, inplace=True)  # drop first column which is dummy
dfTheta.drop(index=dfTheta.index[0], axis=0, inplace=True)  # drop first row which contains NaNs


def VolumetricAverage(P_HF,Volume):
    return np.sum((np.multiply(np.squeeze(P_HF),Volume)))/np.sum(Volume)


V_name = os.getcwd() + '/V_SD'
V_ph = Ofpp.parse_internal_field(V_name)


fname = os.getcwd() + '/U_HF'
U_internal_LES = Ofpp.parse_internal_field(fname)
U_x_LES, U_y_LES, U_z_LES = np.split(U_internal_LES,3,axis=1)


fname = os.getcwd() + '/vort_HF'
Vort_internal_LES = Ofpp.parse_internal_field(fname)
Vort_x_LES, Vort_y_LES, Vort_z_LES = np.split(Vort_internal_LES,3,axis=1)
Diff_Vort_KE = np.abs(Vort_x_LES)
AV_Diff_Vort_KE = VolumetricAverage(Diff_Vort_KE, V_ph)

fname = os.getcwd() + '/U_KE'
U_internal_KE = Ofpp.parse_internal_field(fname)
U_x_KE, U_y_KE, U_z_KE = np.split(U_internal_KE,3,axis=1)


DiffUSt_KE = np.abs(U_x_KE - U_x_LES)
AV_DiffUSt_KE = VolumetricAverage(DiffUSt_KE, V_ph)

list_subfolders = [f.name for f in os.scandir(os.getcwd()) if f.is_dir()]
# print(list_subfolders)

list_subfoders_num = []
for folder in list_subfolders:
    if folder != '0':
        if (folder.isnumeric()):
            list_subfoders_num.append(folder)
list_subfoders_num.sort()
res = [eval(i) for i in list_subfoders_num]


if res != []:

    selectedFolder = str(np.max(res))
    fname = selectedFolder + '/U'
    U_internal_ = Ofpp.parse_internal_field(fname)
    U_x_, U_y_, U_z_ = np.split(U_internal_,3,axis=1)

    DiffUSt = np.abs(U_x_ - U_x_LES)
    VA_DiffUSt = VolumetricAverage(DiffUSt, V_ph)
    NormErrorUSt = VA_DiffUSt/AV_DiffUSt_KE

    fname = selectedFolder + '/vorticity'
    Vort_internal_ = Ofpp.parse_internal_field(fname)
    Vort_x_, Vort_y_, Vort_z_ = np.split(Vort_internal_,3,axis=1)
    
    Diff_Vorticity = np.abs(Vort_x_ - Vort_x_LES)
    VA_Diff_Vorticity = VolumetricAverage(Diff_Vorticity, V_ph)
    NormErrorUSF = VA_Diff_Vorticity/AV_Diff_Vort_KE

else:
    NormErrorUSt = 1000
    NormErrorUSF = 1000

NormErrorUav = (1 * NormErrorUSt + NormErrorUSF) / 2

differenceU = U_x_LES - U_x_
Ux_grid_2 = np.reshape(differenceU,(Nz,Ny))

plotSecondaryFlow = Diff_Vorticity

################################################

fig, ax = plt.subplots(figsize=figSize)
Ub = 0.4820072072378321

normalisedError = (U_x_ - U_x_LES) / VolumetricAverage((U_x_KE - U_x_LES), V_ph)
normalisedErrorGrid = np.reshape(normalisedError, (Nz,Ny))

vmin = -1
vmax = 1
delta = 0.25

cs = ax.contourf(y, z, normalisedErrorGrid, 35, cmap='RdBu_r', vmin=vmin, vmax=vmax)

# print(normalisedErrorGrid)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)


nNy = 11
nNz = 11
nnNy = np.floor(Ny / nNy)
nnNz = np.floor(Nz / nNz)

yQuiver = np.zeros(nNy)
zQuiver = np.zeros(nNz)
UyQuiver= np.zeros((nNz, nNy))
UzQuiver = np.zeros((nNz, nNy))

for i in range(nNy):
    for j in range(nNz):
        iy = int(np.floor(nnNy / 2 + 1 + nnNy * i))
        jz = int(np.floor(nnNz / 2 + 5 + nnNz * j))
        yQuiver[i] = y[iy]
        zQuiver[j] = z[jz]
        UyQuiver[j,i] = Uy_grid[jz, iy]
        UzQuiver[j,i] = Uz_grid[jz, iy]
        # UyQuiver[j,i] = 0.5
        # UzQuiver[j,i] = 0

### FOR COLOURBAR ###
norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
norm.autoscale(plotSecondaryFlow)
colormap = matplotlib.cm.turbo  # RdBu_r
# sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)   
# sm.set_array([])

ax.quiver(yQuiver,zQuiver, UyQuiver, UzQuiver,
           scale=0.07, 
           scale_units='y', 
           # width=0.006, 
           # headwidth=8, 
           # headlength=6,
           # headaxislength=4,
           edgecolor='k',
           linewidth=0.25,
           color=colormap(norm(plotSecondaryFlow)))

# plt.axis('equal')
ax.set_aspect('equal', 'box')

ax.set_xticks([0, 0.5, 1], ['0', '0.5', '1'])
ax.set_yticks([0, 0.5, 1], ['0', '0.5', '1'])
# ax.set_yticks([0, 0.5, 1], [])

ax.set_xlabel(r"$y/h$")
ax.set_ylabel(r"$z/h$")

ax.set_title(r'$\theta_{1} = ' + str(np.round(dfTheta.theta_1.iloc[0], 3)) +
 r',\theta_{2} = ' + str(np.round(dfTheta.theta_2.iloc[0], 3)) +
  r',\theta_{3} = ' + str(np.round(dfTheta.theta_3.iloc[0], 3)) + 
  '$\n$j_{1} = ' + str(np.round(NormErrorUSt, 3)) +
 r',j_{2} = ' + str(np.round(NormErrorUSF, 3)) +
  r',J = ' + str(np.round(NormErrorUav, 3)) + '$',
    loc='left')

# plt.title(sub +title_plot)
# plt.title(sub + 'Standard ' +  r'$k-\epsilon$')
# # plt.xlim([0,1]) 
# # plt.ylim([0,1]) 

# vmin= -1  #minimum value to show on colobar
# vmax = 1 #maximum value to show on colobar
norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
cax, _  = matplotlib.colorbar.make_axes(fig.gca())
cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='RdBu_r', norm=norm, label=r'$\frac{\langle u \rangle - \langle u \rangle_{LES}}{U_b}$')


ticks = np.arange(vmin, vmax + delta, delta)
ticklabels = np.round(np.arange(vmin, vmax+delta, delta), 2)

cbar.set_ticks(ticks)
cbar.set_ticklabels(ticklabels)

# cbar.set_ticks([-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0]) #[0, 0.01, 0.03, 0.05] [0, 0.5, 1, 1.5] [0, 5, 10, 15] [0, 0.3, 0.7, 1]
# cbar.set_ticklabels([-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0])

cbar.set_label(r'$\frac{\langle u \rangle - \langle u \rangle_{LES}}{\epsilon_{{u}_{k-\epsilon - LES}}}$', rotation=0, loc='bottom')
# cbar.set_label(r'$\frac{\langle v \rangle - \langle v \rangle_{LES}}{U_b}$', rotation=0, loc='top')

# plt.plot([np.pi/2-0.48,np.pi/2],[0,0], 'k-',linewidth=10)
nameTosave = selectedFolder + '.png'
fig.savefig(nameTosave, dpi=300, transparent=False, bbox_inches='tight')




