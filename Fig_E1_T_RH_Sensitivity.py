# Project: Regional Variation in the Role of Humidity on City-level Heat-Related Mortality
# Author: Dr.GUO Qiang, The University of Tokyo
# Contact: qiang@rainbow.iis.u-tokyo.ac.jp
# This script is developed and tested in Python 3.7.12 and Linux platform, and can be run also in Windows with
# proper environment settings.
# Description:
# This script is the to plot the sensitivity of indicators to Tair and RH (Extended Data Fig.1)


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import gridspec
from matplotlib import colorbar
import matplotlib as mpl

# Self Defined Packages
import HSI_Calculation_Algorithms as HSI
import newcolors as ncl

mpl.use('Qt5Agg')
plt.close('all')
one_clm = 3.46
two_clm = 7.08


hsi_tick_adjust = [r'T$_{\mathrm{w}}$', r'T$_{\mathrm{S}}$', r'T$_{\mathrm{WBG}}$', r'T$_{\mathrm{sWBG}}$',
                   r'H$_{\mathrm{x}}$', 'APT', 'UTCI', 'HI']

# Some fixed constant and variables (air pressure, wind speed, and solar radiation)
# for HSIs that consider other variables
ws = 5.0  # m/s
ssrd = 1000.0  # W/m2
pres = 101325.0  # Pa


#######################################################################################################################
# Prepare the range of T and H, and the array of HSIs
T_20_50 = np.linspace(20, 50, 61)  # Celsius degree, ranges from 20 degree to 50 degree
H_0_100 = np.linspace(0, 100, 101)  # RH%, ranges from 0 to 100%
HSIs_contour = np.zeros((8, 61, 101)) * np.NaN  # Tw, Ts, Twbg, Tswbg, Xx, APT, UTCI, HI

# Loop for T and RH, and calculating corresponding values for HSIs
for tt in range(61):
    for hh in range(101):

        at = T_20_50[tt]
        rh = H_0_100[hh]

        HSIs_contour[0, tt, hh] = HSI.get_tw_stull(at, rh)   # calculate Tw
        HSIs_contour[1, tt, hh] = HSI.get_ts(HSIs_contour[0, tt, hh], rh)  # calculate Ts
        HSIs_contour[2, tt, hh] = HSI.get_wbgt(at, rh, ws, ssrd)  # calculate Twbg
        HSIs_contour[3, tt, hh] = HSI.get_swbgt(at, rh, pres)  # calculate Tswbg
        HSIs_contour[4, tt, hh] = HSI.get_hx(at, rh, pres)  # calculate Hx
        HSIs_contour[5, tt, hh] = HSI.get_apt(at, rh, ws, pres)  # calculate APT
        HSIs_contour[6, tt, hh] = HSI.get_utci(at, at, ws, rh, pres)  # calculate UTCI
        HSIs_contour[7, tt, hh] = HSI.get_hi(at, rh)  # calculate HI

print('HSI matrix generation finished.')

#######################################################################################################################
# set the colormap
cmap_hsi = ncl.newcmp('t10_r')
bounds_hsi = np.linspace(10, 50, 9)
norm_hsi = mpl.colors.BoundaryNorm(bounds_hsi, cmap_hsi.N)
# set the subplot number
sp_od = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']


# define the function to plot the contour
def plot_contour(hsi, hsi_array, ref_array, title, yl, xl):

    # plot the shaded color for the corresponding HSI
    plt.pcolormesh(T_20_50,
                   H_0_100,
                   hsi_array.T,
                   norm=norm_hsi,
                   cmap=cmap_hsi,
                   alpha=1,
                   shading='auto')

    # plot the contour line of Tw, and used as a reference
    ct = plt.contour(T_20_50,
                     H_0_100,
                     ref_array.T,
                     levels=np.linspace(10, 50, 9, dtype=int),
                     colors='k',
                     linewidths=1.1)

    # add value labels for the Tw contour
    fmt = {}
    nums = np.linspace(10, 50, 9, dtype=int)
    for l, s in zip(ct.levels, nums):
        fmt[l] = s
    plt.clabel(ct, ct.levels, inline=True, fmt=fmt, fontsize=8)

    # add title and No. for subplot
    plt.title(title, fontsize=11)
    plt.text(19, 103, sp_od[hsi], fontsize=9)

    # add ticklabels
    if yl:
        plt.ylabel('RH (%)', fontsize=9)
        plt.yticks(np.linspace(0, 100, 5, dtype=int),
                   np.linspace(0, 100, 5, dtype=int),
                   fontsize=8)
    else:
        plt.yticks([])

    if xl:
        plt.xlabel(r'T$_{\mathrm{air}}$ ($^{\mathrm{o}}$C)', fontsize=9)
        plt.xticks(np.linspace(20, 50, 4, dtype=int),
                   np.linspace(20, 50, 4, dtype=int),
                   fontsize=8)
    else:
        plt.xticks([])


# Start to plot the figure
fig1 = plt.figure(1, figsize=(two_clm, 5.4), constrained_layout=True)
gs = gridspec.GridSpec(figure=fig1,
                       nrows=3,
                       ncols=5,
                       height_ratios=[1,
                                      1,
                                      1],
                       width_ratios=[1, 1, 1, 1, 0.08])

# Contour
for hsi in range(8):
    if hsi % 4 == 0:
        yl = True
    else:
        yl = False

    if hsi > 3:
        xl = True
    else:
        xl = False

    ax = fig1.add_subplot(gs[hsi // 4, hsi % 4])
    plot_contour(hsi, HSIs_contour[hsi, :, :], HSIs_contour[0, :, :],
                 hsi_tick_adjust[hsi], yl, xl)

# colobar for contour
ax = fig1.add_subplot(gs[:2, 4])
cb = colorbar.ColorbarBase(ax,
                           cmap=cmap_hsi,
                           norm=norm_hsi,
                           orientation='vertical',
                           extend='both')
cb.set_label(label=r'HSI value ($^{\mathrm{o}}$C)', size=9)
cb.ax.set_yticks(np.linspace(10, 50, 9),
                 np.linspace(10, 50, 9, dtype='int'),
                 size=8)

# Bar
sns_rh = np.zeros((5, 8)) * np.NaN
sns_rh[0, :] = np.nanmin(HSIs_contour, axis=(1, 2))  # get the minimum values for each HSI in the whole range
sns_rh[1, :] = np.nanmax(HSIs_contour, axis=(1, 2))  # get the maximum values for each HSI in the whole range

# get the value difference when RH change from 0 to 100% when Tair=25 degree for each HSI
sns_rh[2, :] = np.nanmax(HSIs_contour[:, 10, :], axis=1) \
               - np.nanmin(HSIs_contour[:, 10, :], axis=1)

# get the value difference when RH change from 0 to 100% when Tair=30 degree for each HSI
sns_rh[3, :] = np.nanmax(HSIs_contour[:, 20, :], axis=1) \
               - np.nanmin(HSIs_contour[:, 20, :], axis=1)

# get the value difference when RH change from 0 to 100% when Tair=35 degree for each HSI
sns_rh[4, :] = np.nanmax(HSIs_contour[:, 30, :], axis=1) \
               - np.nanmin(HSIs_contour[:, 30, :], axis=1)

ax = fig1.add_subplot(gs[2, :4])
# get the relative value on HSIs change by dividing the whole range for each HSI
sns = {
    r'T$_{\mathrm{air}}$ = 25 $^{\mathrm{o}}$C': sns_rh[2, :] / (sns_rh[1, :] - sns_rh[0, :]) * 100,
    r'T$_{\mathrm{air}}$ = 30 $^{\mathrm{o}}$C': sns_rh[3, :] / (sns_rh[1, :] - sns_rh[0, :]) * 100,
    r'T$_{\mathrm{air}}$ = 35 $^{\mathrm{o}}$C': sns_rh[4, :] / (sns_rh[1, :] - sns_rh[0, :]) * 100,
}

x = np.arange(len(hsi_tick_adjust))
width = 0.27
multiplier = 0

# plot the results using barchart
for attribute, measurement in sns.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute, edgecolor='k', linewidth=0.2)
    ax.bar_label(rects, label_type='edge', fontsize=8, fmt='%.0f')
    multiplier += 1

# add the ticklabel, legend, etc.
ax.set_ylabel('Fraction (%)', fontsize=9)
plt.xticks(x + width, hsi_tick_adjust, fontsize=8)
plt.yticks(np.linspace(0, 75, 4, dtype=int),
              np.linspace(0, 75, 4, dtype=int), fontsize=8)
ax.legend(ncol=3, fontsize='small', loc='upper right')
ax.set_ylim(0, 75)
ax.text(-0.5, 76, 'i', fontsize=9)


plt.show()
# fig1.savefig('Fig_E1.jpg',
#              format='jpg',
#              dpi=1200)

print('All Finished.')


