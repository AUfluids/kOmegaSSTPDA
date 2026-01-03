#!/usr/bin/env python3
# Mario J. Rincon Ph.D.
###############################
# PLOTTING OPENFOAM RESIDUALS #
###############################
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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

path = os.getcwd()

variableName = np.array(['Ux_0', 'Uy_0', 'Uz_0', 'p_0', 'k_0', 'omega_0'])
legendName = np.array(['$u$', '$v$', '$w$', '$p$', '$k$', r'$\omega$'])

print('Plotting linear residuals')
fig, ax = plt.subplots(figsize=figSize, ncols=1, nrows=1, sharey=True)
for i in range(len(variableName)):
    headers = ['Time', variableName[i]]
    #print(variableName[i])
    df = pd.read_csv(str(path) + '/logs/' + str(variableName[i]), delimiter=r'\s+', names=headers)

    print('Plotting ' + str(variableName[i]) + ' residuals')
    ax.plot(df['Time'], df[variableName[i]], label=legendName[i], lw=0.75)

ax.set_ylabel("Initial residual")
ax.set_xlabel("Iteration")
ax.set_yscale('log')
# plt.grid(True, which="both", ls=':', lw=1)
ax.legend(loc='best')
# fig.savefig('linear.png', dpi=300, transparent=False)
# fig.clf()

# ########################
# # CONTINUITY RESIDUALS #
# ########################

# variableName = ['contGlobal_0', 'contCumulative_0', 'contLocal_0']
# legendName = ['Global', 'Cumulative', 'Local']

# print('Plotting continuity residuals')
# for i in range(len(variableName)):
#     headers = ['Time', variableName[i]]
#     df = pd.read_csv(str(path) + '/logs/' + str(variableName[i]), delimiter=r'\s+', names=headers)

#     if i == 0:
#         df1 = df
#     elif i == 1:
#         df2 = df
#     else:
#         df3 = df

# # fig, ax1 = plt.subplots()

# color = 'tab:red'
# ax[1].set_xlabel('Iteration')
# ax[1].set_ylabel(legendName[0], color=color)
# ax[1].plot(df['Time'], df2['contCumulative_0'], color=color, label=legendName[0], lw=0.75)
# #ax1.hlines(y=0, xmin=0, xmax=df['Time'][len(df['Time']) - 1], color='tab:green', ls='dotted', lw=1.0, zorder=10)
# ax[1].tick_params(axis='y', labelcolor=color)
# # ax1.grid(True, which='both', ls=':')

# ax2 = ax[1].twinx()  # instantiate a second axes that shares the same x-axis

# color = 'tab:green'
# ax2.set_ylabel(legendName[1], color=color)  # we already handled the x-label with ax1
# ax2.plot(df['Time'], df1['contGlobal_0'], color=color, label=legendName[1], lw=0.5)
# ax2.hlines(y=0, xmin=0, xmax=df['Time'][len(df['Time']) - 1], color='tab:red', ls='dotted', lw=1.0, zorder=10)
# ax2.tick_params(axis='y', labelcolor=color)

# fig.tight_layout()  # otherwise the right y-label is slightly clipped
# # plt.savefig('cont.png', dpi=300, transparent=False)
# # plt.clf()

# print('Plotting local residuals')
# # PLOT LOCAL CONTINUITY RESIDUALS
# data = ax[2].plot(df['Time'], df3['contLocal_0'], label=legendName[2], lw=0.75)
# ax[2].set_ylabel(legendName[2])
# ax[2].set_xlabel("Iteration")
# # plt.grid(True, ls=':')
# # fig.savefig('local.png', dpi=300, transparent=False)
fig.savefig('residuals.png', dpi=300, transparent=False)
# plt.clf()