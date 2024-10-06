# Project: Regional Variation in the Role of Humidity on City-level Heat-Related Mortality
# Author: Dr.GUO Qiang, The University of Tokyo
# Contact: qiang@rainbow.iis.u-tokyo.ac.jp
# This script is developed and tested in Python 3.7.12 and Linux platform, and can be run also in Windows with
# proper environment settings.
# Description:
# This script is the to plot the feature importance (Figure 4 in the manuscript)


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import colorbar
from matplotlib import gridspec
import matplotlib as mpl
from scipy.stats import gaussian_kde
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import statsmodels.api as sm
from matplotlib.transforms import offset_copy

# Self Defined Packages
import newcolors as ncl

mpl.use('Qt5Agg')
plt.close('all')
one_clm = 3.46
two_clm = 7.08

# Load some pre-calculated data for plotting the figure
with np.load('fig45_feature_importance.npz') as file:
    imp_500 = file['imp_500']   # feature importance obtained from Random Forest
    cf_mtx_500 = file['cf_mtx_500']  # confusion matrix obtained from Random Forest
    mcc_lon = file['mcc_lon']  # lon for each cities
    mcc_lat = file['mcc_lat']  # lat for each cities

imp_mean = np.mean(imp_500, 1)   # average for 500 implementations
imp_rank = np.argsort(imp_mean)
imp_show = np.sort(imp_mean)
err = np.std(imp_500, axis=1)
err_show = err[imp_rank]

# load the input of 12 features for random forest
rf_csv = pd.read_csv('Regional Variation in the Role of Humidity on City-level Heat-Related Mortality_TableS4_QG.csv',
                     header=0)
ctrh = np.array(rf_csv['corr_t_rh'])

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

row_name = np.array([r'T$_{\mathrm{mean}}$', r'T$_{\mathrm{std}}$', r'RH$_{\mathrm{mean}}$',
                     r'RH$_{\mathrm{std}}$', r'C$_{\mathrm{T-RH}}$', 'Elv',
                     r'D$_{\mathrm{coast}}$', 'Lat', 'GDP',
                     'KG', 'Ctry', 'Rgn'])

hsi_tick_humi_weight = [r'T$_{\mathrm{air}}$',
                        'HI', 'UTCI', 'APT', r'H$_{\mathrm{x}}$',
                        r'T$_{\mathrm{sWBG}}$', r'T$_{\mathrm{WBG}}$',
                        r'T$_{\mathrm{s}}$', r'T$_{\mathrm{w}}$']

#######################################################################################################################
# Plot the figure 4
fig1 = plt.figure(1, figsize=(two_clm, 6.1))
gs = gridspec.GridSpec(figure=fig1,
                       nrows=3,
                       ncols=3,
                       height_ratios=[1,
                                      0.05,
                                      0.6],
                       width_ratios=[1, 1, 0.05])
od_no = ['a', 'b', 'c', 'd', 'e']

set_pos = [[0.070, 0.47, 0.22, 0.49],  # a: feature importance
           [0.294, 0.48, 0.703, 0.48],  # b: map
           [0.294, 0.47, 0.703, 0.015],  # colorbar 1
           [0.10, 0.08, 0.80, 0.3],  # c: hh vs corr
           [0.91, 0.08, 0.01, 0.3]]  # colorbar 2

# feature importance, subplot a
ax = fig1.add_subplot(gs[0: 2, 0])
ax.set_position(set_pos[0])

b1 = ax.barh(np.linspace(1, 12, 12), imp_show,
             align='center', edgecolor='k', xerr=err_show, linewidth=0.2,
             color=['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown',
                    'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'lightcoral', 'teal'][::-1])
plt.yticks(np.linspace(1, 12, 12), labels=row_name[imp_rank], fontsize=8)
plt.xticks([0, 5, 10, 15], labels=[0, 5, 10, 15], fontsize=8)
ax.set_xlabel('Feature importance (%)', fontsize=9)
# ax.set_ylabel('Feature')
ax.bar_label(b1, label_type='edge', fontsize=8, fmt='%.1f')
plt.axis([0, 17, 0.5, 12.5])
plt.text(0, 12.6, 'a', fontsize=9)


# global map of corr_TRH, subplot b
ax = fig1.add_subplot(gs[0, 1], projection=ccrs.Robinson())
ax.set_position(set_pos[1])

cmap_corr = ncl.newcmp('t8_r')
# cmap_corr = plt.cm.get_cmap('coolwarm', 8)
bounds_corr = np.linspace(-0.6, 0.6, 9)
norm_corr = mpl.colors.BoundaryNorm(bounds_corr, cmap_corr.N)

ax.set_global()
ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor=[0.8, 0.8, 0.8], edgecolor='white', zorder=0)
ax.add_feature(cfeature.BORDERS.with_scale('50m'), edgecolor='white', linewidth=0.3, zorder=1)
plt.scatter(mcc_lon,
            mcc_lat,
            c=ctrh,
            s=8,
            norm=norm_corr,
            cmap=cmap_corr,
            alpha=0.8,
            zorder=3,
            edgecolors='k',
            linewidths=0.1,
            transform=ccrs.PlateCarree())
plt.text(-1.65e+07, 8.1e+06, 'b', fontsize=9)

# colorbar, subplot b
ax = fig1.add_subplot(gs[1, 1])
ax.set_position(set_pos[2])
cb = colorbar.ColorbarBase(ax,
                           cmap=cmap_corr,
                           norm=norm_corr,
                           orientation='horizontal',
                           extend='both',
                           ticks=[-0.6, -0.3, 0, 0.3, 0.6],
                           label=r'C$_{\mathrm{T-RH}}$')
cb.ax.set_xticklabels([-0.6, -0.3, 0, 0.3, 0.6], fontsize=8)

# density plot, subplot c
ax = fig1.add_subplot(gs[2, :2])
ax.set_position(set_pos[3])
cmap_ovlp = plt.cm.get_cmap('Spectral_r')
norm_ovlp = colors.Normalize(vmin=0,
                             vmax=1)
plt.axvline(0, color='k', linewidth=1.1, linestyle='--')

for i in range(9):
    i_index = np.where(rp_hsi_weight == i)[0]
    x = np.array(ctrh[i_index])
    y = rp_hsi_weight[i_index]
    z = gaussian_kde(x)(x)  # the second (x) is use the previous (gaussian_kde(x)) as a function
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    plt.scatter(x, y, c=z, s=30, cmap=cmap_ovlp, norm=norm_ovlp)

plt.yticks(np.linspace(0, 8, 9), labels=hsi_tick_humi_weight, fontsize=9)
plt.xticks(np.linspace(-0.8, 0.8, 9), fontsize=8)

plt.axis([-0.9, 0.9, -0.5, 8.5])
plt.ylabel('BFI')
plt.xlabel(r'C$_{\mathrm{T-RH}}$')
plt.text(-0.9, 8.7, 'c', fontsize=9)

# colorbar, subplot c
ax = fig1.add_subplot(gs[2, 2])
ax.set_position(set_pos[4])
cb = colorbar.ColorbarBase(ax,
                           cmap=cmap_ovlp,
                           norm=norm_ovlp,
                           orientation='vertical',
                           label='Density',
                           ticks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
cb.ax.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=8)

# plt.show()

# fig1.savefig('Fig_4.jpg',
#              format='jpg',
#              dpi=1200)


print('All Finished.')
