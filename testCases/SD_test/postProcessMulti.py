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

# coordinates = np.linspace(0, 0.5, Ny)
coordinates = (coordinates - coordinates.min()) / (coordinates.max() - coordinates.min()) * 0.5
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




# plt.clf()
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
# dfTheta = pd.read_csv(os.getcwd() + '/constant/theta_values.csv')
# dfTheta.drop(dfTheta.columns[0], axis=1, inplace=True)  # drop first column which is dummy
# dfTheta.drop(index=dfTheta.index[0], axis=0, inplace=True)  # drop first row which contains NaNs


def VolumetricAverage(P_HF,Volume):
    return np.sum((np.multiply(np.squeeze(P_HF),Volume)))/np.sum(Volume)


fname = os.getcwd() + '/!RefData/C_SD'
C_internal = Ofpp.parse_internal_field(fname)
C_x, C_y, C_z = np.split(C_internal, 3, axis=1)


V_name = os.getcwd() + '/!RefData/V_SD'
V_ph = Ofpp.parse_internal_field(V_name)


fname = os.getcwd() + '/!RefData/U_HF'
U_internal_LES = Ofpp.parse_internal_field(fname)
U_x_LES, U_y_LES, U_z_LES = np.split(U_internal_LES, 3, axis=1)


fname = os.getcwd() + '/!RefData/vort_HF'
Vort_internal_LES = Ofpp.parse_internal_field(fname)
Vort_x_LES, Vort_y_LES, Vort_z_LES = np.split(Vort_internal_LES, 3, axis=1)
# Diff_Vort_LES = np.abs(Vort_x_LES)
# AV_Diff_Vort_LES = VolumetricAverage(Diff_Vort_LES, V_ph)

fname = os.getcwd() + '/!RefData/vort_KOSST'
Vort_internal_KOSST = Ofpp.parse_internal_field(fname)
Vort_x_KOSST, Vort_y_KOSST, Vort_z_KOSST = np.split(Vort_internal_KOSST, 3, axis=1)
Diff_Vort_KOSST = np.abs(Vort_x_KOSST - Vort_x_LES)
AV_Diff_Vort_KOSST = VolumetricAverage(Diff_Vort_KOSST, V_ph)

# fname = os.getcwd() + '/U_KE'
# U_internal_KE = Ofpp.parse_internal_field(fname)
# U_x_KE, U_y_KE, U_z_KE = np.split(U_internal_KE, 3,axis=1)

fname = os.getcwd() + '/!RefData/U_KOSST'
U_internal_KOSST = Ofpp.parse_internal_field(fname)
U_x_KOSST, U_y_KOSST, U_z_KOSST = np.split(U_internal_KOSST, 3,axis=1)

DiffUSt_KOSST = np.abs(U_x_KOSST - U_x_LES)
AV_DiffUSt_KOSST = VolumetricAverage(DiffUSt_KOSST, V_ph)

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
    U_x_, U_y_, U_z_ = np.split(U_internal_, 3, axis=1)

    DiffUSt = np.abs(U_x_ - U_x_LES)
    VA_DiffUSt = VolumetricAverage(DiffUSt, V_ph)
    NormErrorUSt = VA_DiffUSt/AV_DiffUSt_KOSST

    fname = selectedFolder + '/vorticity'
    Vort_internal_ = Ofpp.parse_internal_field(fname)
    Vort_x_, Vort_y_, Vort_z_ = np.split(Vort_internal_,3,axis=1)
    
    Diff_Vorticity = np.abs(Vort_x_ - Vort_x_LES)
    VA_Diff_Vorticity = VolumetricAverage(Diff_Vorticity, V_ph)
    NormErrorUSF = VA_Diff_Vorticity / AV_Diff_Vort_KOSST
    # NormErrorUSF = VA_Diff_Vorticity / AV_Diff_Vort_LES

    # print(AV_Diff_Vort_KOSST, AV_Diff_Vort_LES, VA_Diff_Vorticity)

else:
    NormErrorUSt = 1000
    NormErrorUSF = 1000

NormErrorUav = (1 * NormErrorUSt + NormErrorUSF) / 2

fileName = os.getcwd() +'/objFunResults.txt'
with open(fileName, 'w') as f:
        f.write('j1 = ')
        f.write('%.8f' %(NormErrorUSt))
        f.write('\nj2 = ')
        f.write('%.8f' %(NormErrorUSF))
        f.write('\nJ = ')
        f.write('%.8f' %(NormErrorUav))

################################################

fig, ax = plt.subplots(figsize=figSizeMedium, ncols=3, nrows=1, sharey=True)
Ub = 0.4820072072378321

# normalisedErrorUSt = (U_x_ - U_x_LES) / VolumetricAverage((U_x_KOSST - U_x_LES), V_ph)
# normalisedErrorUSt = (U_x_ - U_x_LES) / (U_x_KOSST - U_x_LES)

# print( Vort_x_LES.shape[0])

normalisedVelocityField = np.reshape(U_x_ / Ub, (Nz, Ny))
normalisedErrorUStGrid = np.reshape((U_x_ - U_x_LES) / np.mean(np.abs((U_x_KOSST - U_x_LES))), (Nz, Ny))
normalisedErrorVortGrid = np.reshape((Vort_x_ - Vort_x_LES) / AV_Diff_Vort_KOSST, (Nz, Ny))
# normalisedErrorVortGrid = np.reshape((Vort_x_ - Vort_x_LES) / (Vort_x_KOSST - Vort_x_LES), (Nz, Ny))

vmin = -1
vmax = 1
delta = 0.5

norm1 = matplotlib.colors.Normalize(vmin=0, vmax=1.5)
norm2 = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
norm3 = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

cs = ax[0].contourf(-y, z, normalisedVelocityField, levels=np.linspace(0, 1.5, 35), cmap='turbo', vmin=0, vmax=1.5, norm=norm1, extend='both')
csUst = ax[1].contourf(-y, z, normalisedErrorUStGrid, levels=np.linspace(vmin, vmax, 35), cmap='RdBu_r', vmin=vmin, vmax=vmax, norm=norm2, extend='both')
csVort = ax[2].contourf(-y, z, normalisedErrorVortGrid, levels=np.linspace(vmin, vmax, 35), cmap='PuOr_r', vmin=vmin, vmax=vmax, norm=norm3, extend='both')

cbar1 = fig.colorbar(cs, cmap='turbo', pad=0.15, orientation='horizontal', norm=norm1)
cbar2 = fig.colorbar(csUst, cmap='RdBu_r', pad=0.15, orientation='horizontal', norm=norm2)
cbar3 = fig.colorbar(csVort, cmap='PuOr_r', pad=0.15, orientation='horizontal', norm=norm3)

ticks = np.arange(vmin, vmax + delta, delta)
ticklabels = np.round(np.arange(vmin, vmax + delta, delta), 2)

cbar1.set_ticks([0, 0.5, 1.0, 1.5])
cbar1.set_ticklabels([0, 0.5, 1.0, 1.5])
cbar1.set_label(r'$\frac{\langle u \rangle}{U_b}$', rotation=0, loc='right')

cbar2.set_ticks(ticks)
cbar2.set_ticklabels(ticklabels)
cbar2.set_label(r'$\frac{\langle u \rangle - \langle u \rangle_{LES}}{\epsilon_{u}}$', rotation=0, loc='right')

cbar3.set_ticks(ticks)
cbar3.set_ticklabels(ticklabels)
cbar3.set_label(r'$\frac{\langle \omega_{x} \rangle - \langle \omega_{x} \rangle_{LES}}{\epsilon_{\omega_{x}}}$', rotation=0, loc='right')

#### STREAMLINE PLOT ####
from scipy.interpolate import griddata

surfaceVelocity = np.sqrt(U_y_ ** 2 + U_z_ ** 2)

resolution = 100

H = 0.5

xx = -C_y / H
yy = C_z / H + 2 * H

x = np.linspace(-2 * H, H, resolution)
yS = np.linspace(0, 2 * H, resolution)

xi, yi = np.meshgrid(x, yS)

points = np.concatenate((xx, yy), axis=1)

gu = griddata(points, -U_y_, (xi, yi), method='linear')
gv = griddata(points, U_z_, (xi, yi), method='linear')
surfaceVelocityPlot = griddata(points, surfaceVelocity, (xi, yi), method='nearest')

gu = gu.reshape(resolution, resolution)
gv = gv.reshape(resolution, resolution)
surfaceVelocityPlot = surfaceVelocityPlot.reshape(resolution, resolution)

lw = surfaceVelocityPlot / surfaceVelocityPlot.max()

stream = ax[0].streamplot(xi, yi, gu, gv, density=2, linewidth=lw, color='k', arrowsize=0.5, arrowstyle='fancy')

# nNy = 48
# nNz = 48
# nnNy = int(np.floor(Ny / nNy))
# nnNz = int(np.floor(Nz / nNz))

# yQuiver = C_y[::nnNy]
# zQuiver = C_z[::nnNz]

# UyQuiver = U_y_[::nnNy]
# UzQuiver = U_z_[::nnNz]

# ax[0].quiver(-(yQuiver * 2), (zQuiver * 2) + 1, -UyQuiver, UzQuiver,
#            scale=0.08, 
#            scale_units='y', 
#            # width=0.006, 
#            # headwidth=8, 
#            # headlength=6,
#            # headaxislength=4,
#            edgecolor='k',
#            linewidth=0.25,
#            color='w')  # colormap(norm(normalisedErrorVortGrid)))

# plt.axis('equal')
ax[0].set_aspect('equal', 'box')
ax[1].set_aspect('equal', 'box')
ax[2].set_aspect('equal', 'box')

ax[0].set_xticks([-1, -0.5, 0], ['0', '0.5', '1'])
ax[1].set_xticks([-1, -0.5, 0], ['0', '0.5', '1'])
ax[2].set_xticks([-1, -0.5, 0], ['0', '0.5', '1'])

ax[0].set_yticks([0, 0.5, 1], ['0', '0.5', '1'])
# ax.set_yticks([0, 0.5, 1], [])

ax[0].set_xlabel(r"$y/h$")
ax[1].set_xlabel(r"$y/h$")
ax[2].set_xlabel(r"$y/h$")

ax[0].set_ylabel(r"$z/h$")

## to automate the plot titples
# for i, v in enumerate(dfTheta):
#     if i ==0:
#         supTitle1 = r'$\theta_{' + str(i) + '} = ' + str(np.round(eval(str('dfTheta.theta_' + str(i) + '.iloc[0]', 3)))
#     else:
#         supTitle = np.concatenate(supTitle1)
# fig.suptitle(supTitle + '$', y=0.72, fontsize=7)


# fig.suptitle(r'$C_{1} = ' + str(np.round(dfTheta.theta_1.iloc[0], 3)) +
#  r',C_{2} = ' + str(np.round(dfTheta.theta_2.iloc[0], 3)) +
#  r',C_{3} = ' + str(np.round(dfTheta.theta_3.iloc[0], 3)) + 
#  '$',
#   y=0.72, fontsize=7)

ax[0].set_title(r'$J = ' + str(np.round(NormErrorUav, 3)) + '$', loc='left', fontsize=7)
ax[1].set_title(r'$j_{1} = ' + str(np.round(NormErrorUSt, 3)) + '$', loc='left', fontsize=7)
ax[2].set_title(r'$j_{2} = ' + str(np.round(NormErrorUSF, 3)) + '$', loc='left', fontsize=7)

nameTosave = str(os.path.basename(os.getcwd())) + '_multi.png'
fig.savefig(nameTosave, dpi=300, transparent=False, bbox_inches='tight')