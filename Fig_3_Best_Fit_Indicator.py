# Project: Regional Variation in the Role of Humidity on City-level Heat-Related Mortality
# Author: Dr.GUO Qiang, The University of Tokyo
# Contact: qiang@rainbow.iis.u-tokyo.ac.jp
# This script is developed and tested in Python 3.7.12 and Linux platform, and can be run also in Windows with
# proper environment settings.
# Description:
# This script is the to plot the best-fit indicators (Figure 3 in the manuscript)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import colorbar
from matplotlib import gridspec
from matplotlib import cm
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cfeature


mpl.use('Qt5Agg')
plt.close('all')
one_clm = 3.46
two_clm = 7.08

hsi_tick_humi_weight_no = np.array([r'T$_{\mathrm{air}}$ (0)',
                                    'HI (1)', 'UTCI (2)', 'APT (3)', r'H$_{\mathrm{x}}$ (4)',
                                    r'T$_{\mathrm{sWBG}}$ (5)', r'T$_{\mathrm{WBG}}$ (6)',
                                    r'T$_{\mathrm{s}}$ (7)', r'T$_{\mathrm{w}}$ (8)'])
hsi_tick_humi_weight = np.array([r'T$_{\mathrm{air}}$',
                                 'HI', 'UTCI', 'APT', r'H$_{\mathrm{x}}$',
                                 r'T$_{\mathrm{sWBG}}$', r'T$_{\mathrm{WBG}}$',
                                 r'T$_{\mathrm{s}}$', r'T$_{\mathrm{w}}$'])

rg_uqe = ['North Europe', 'Central Europe', 'South Europe',
          'East Asia', 'South-East Asia', 'Middle-East Asia',
          'North America', 'Central America', 'South America',
          'South Africa', 'Australia']

# Load some pre-calculated data for plotting the figure
with np.load('fig3_BFI.npz', allow_pickle=True) as file:
    mcc_lon = file['mcc_lon']
    mcc_lat = file['mcc_lat']
    city_region = file['city_region']

# load the qAIC information in Table S3
qaic_all = np.zeros((9, 739))
qaic_csv = pd.read_csv('Regional Variation in the Role of Humidity on City-level Heat-Related Mortality_TableS3_QG.csv',
                       header=0)
qaic_all[0, :] = np.array(qaic_csv['Tair'])
qaic_all[1, :] = np.array(qaic_csv['HI'])
qaic_all[2, :] = np.array(qaic_csv['UTCI'])
qaic_all[3, :] = np.array(qaic_csv['APT'])
qaic_all[4, :] = np.array(qaic_csv['Hx'])
qaic_all[5, :] = np.array(qaic_csv['TsWBG'])
qaic_all[6, :] = np.array(qaic_csv['TWBG'])
qaic_all[7, :] = np.array(qaic_csv['Ts'])
qaic_all[8, :] = np.array(qaic_csv['Tw'])

# obtain the best-fit indicator, which is the indicators with smallest qAIC
rp_hsi_weight = np.argmin(qaic_all, axis=0)

# calculate the number of best performing cities in each region for each HSI
rp_hsi_num = np.zeros((9, 11), dtype=int)
for rg in range(11):
    rg_city_no = np.where(city_region == rg_uqe[rg])[0]
    for hs in range(9):
        rp_hsi_num[hs, rg] = np.nansum(rp_hsi_weight[rg_city_no] == hs)

######################################################################################################################
# Fig. 1: Global map of representative HSI with humidity weight
# colormap for subplot a
cmap_hsi = plt.cm.get_cmap('autumn_r', 9)
bounds_hsi = np.linspace(0, 9, 10)
norm_hsi = mpl.colors.BoundaryNorm(bounds_hsi, cmap_hsi.N)
cl_bar = cm.get_cmap('tab20c')(np.linspace(0, 1, 15))

fig1 = plt.figure(1, figsize=(two_clm, 7.5), constrained_layout=True)
gs = gridspec.GridSpec(figure=fig1,
                       nrows=3,
                       ncols=1,
                       height_ratios=[1,
                                      0.02,
                                      0.7],
                       width_ratios=[1])

# subplot a
ax = fig1.add_subplot(gs[0, :], projection=ccrs.Robinson())
ax.set_global()
ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor=[0.8, 0.8, 0.8], edgecolor='white', zorder=0)
ax.add_feature(cfeature.BORDERS.with_scale('50m'), edgecolor='white', linewidth=0.3, zorder=1)
plt.scatter(x=mcc_lon,
            y=mcc_lat,
            c=rp_hsi_weight + 0.5,
            s=8,
            norm=norm_hsi,
            cmap=cmap_hsi,
            alpha=0.6,
            zorder=3,
            edgecolors='k',
            linewidths=0.1,
            transform=ccrs.PlateCarree())
plt.text(-1.65e+07, 8.1e+06, 'a', fontsize=10)
plt.title('Indicators with the minimum qAIC', fontsize=10)

# colorbar
ax = fig1.add_subplot(gs[1, :])
cb = colorbar.ColorbarBase(ax,
                           cmap=cmap_hsi,
                           norm=norm_hsi,
                           orientation='horizontal',
                           label='BFI',
                           ticks=list(np.linspace(0.5, 8.5, 9)))
cb.ax.set_xticklabels(list(hsi_tick_humi_weight_no[:]), fontsize=9)
cb.ax.invert_yaxis()

# subplot b
ax = fig1.add_subplot(gs[2, :])
bottom_sum = np.zeros(9)

for rg in range(11):
    pp = ax.bar(np.linspace(1, 9, 9),
                rp_hsi_num[:, rg],  # (2, 9, 11)
                label=rg_uqe[rg],
                bottom=bottom_sum,
                facecolor=cl_bar[rg])
    bottom_sum = bottom_sum + rp_hsi_num[:, rg]
    if rg == 10:
        plt.bar_label(pp, label_type='edge', fontsize=8)

# add legend
ax.legend(ncol=4, fontsize='x-small', loc='upper right')
# add labels, etc.
ax.set(xlim=[0.5, 9.5],
       ylim=[0, 250])
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
ax.set_xlabel('BFI', fontsize=10)
ax.set_ylabel('City numbers', fontsize=9)
plt.xticks(np.linspace(1, 9, 9), hsi_tick_humi_weight, rotation=0, fontsize=9)
plt.text(0.5, 255, 'b', fontsize=10)

plt.show()

# fig1.savefig('Fig_3.jpg',
#              format="jpg",
#              dpi=1200)

print('All Finished.')
