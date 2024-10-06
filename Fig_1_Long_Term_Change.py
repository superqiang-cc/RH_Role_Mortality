# Project: Regional Variation in the Role of Humidity on City-level Heat-Related Mortality
# Author: Dr.GUO Qiang, The University of Tokyo
# Contact: qiang@rainbow.iis.u-tokyo.ac.jp
# This script is developed and tested in Python 3.7.12 and Linux platform, and can be run also in Windows with
# proper environment settings.
# Description:
# This script is the to plot the long-term trend of HSIs (Figure 1 in the manuscript)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
from matplotlib import colorbar
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from mpl_toolkits.basemap import maskoceans

# Self Defined Packages
import newcolors as ncl

mpl.use('Qt5Agg')
plt.close('all')
one_clm = 3.46
two_clm = 7.08

# Load some pre-calculated data for plotting the figure
with np.load('fig1_long_term.npz') as file:
    # size of dec_trend and lin_sig: (9, 2, 2, 721, 1440), each dimension represents
    # axis=0: variables ('at', 'tw', 'swbgt', 'hx',  'apt', 'hi', 'utci', 'sh', 'rh'),
    # axis=1: quantile (95th, 99th),
    # axis=2: max/mean,
    # axis=3: lat,
    # axis=4: lon,
    dec_trend = file['dec_trend']
    lin_sig = file['lin_sig']
    # lon and lat of the land surface grids
    lon_mesh = file['lon_mesh']
    lat_mesh = file['lat_mesh']

########################################################################################################################
# plot trend the HSI consistency:


def plot_trend(to_show, sig, cmp, nrm, tt, od, is_hsi, s=0.5):

    # mask the ocean
    to_show_mask = maskoceans(lon_mesh, lat_mesh, to_show, inlands=True, resolution='h')

    # set the global map
    ax.set_global()
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), edgecolor=[0.2, 0.2, 0.2], linewidth=0.2, zorder=2)
    ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor=[0.8, 0.8, 0.8], edgecolor=[0.2, 0.2, 0.2], linewidth=0.3,
                   zorder=1)
    ax.add_feature(cfeature.OCEAN.with_scale('50m'), facecolor=[0.7, 0.7, 0.7], zorder=1)

    # plot the shaded color
    plt.pcolormesh(lon_mesh[40: 600, :],
                   lat_mesh[40: 600, :],
                   to_show_mask[40: 600, :],
                   cmap=cmp,
                   norm=nrm,
                   shading='auto',
                   zorder=2,
                   transform=ccrs.PlateCarree())

    lon_mesh_sig = lon_mesh.copy()
    lat_mesh_sig = lat_mesh.copy()
    lon_mesh_sig[to_show_mask.mask] = np.NaN
    lat_mesh_sig[to_show_mask.mask] = np.NaN

    # plot the dot, which indicates the significance
    if is_hsi == 1:
        lon_mesh_sig[sig < 0.5] = np.NaN
        lat_mesh_sig[sig < 0.5] = np.NaN

    else:
        lon_mesh_sig[sig > 0.05] = np.NaN
        lat_mesh_sig[sig > 0.05] = np.NaN

    plt.scatter(lon_mesh_sig[40:600:10, ::10],
                lat_mesh_sig[40:600:10, ::10],
                s=s,
                c='k',
                zorder=3,
                facecolor='k',
                edgecolor='none',
                transform=ccrs.PlateCarree())

    # add title and subplot number
    plt.title(tt, fontsize=10)
    plt.text(-1.65e+07, 8.1e+06, od, fontsize=9)


# Start to plot
fig1 = plt.figure(1, figsize=(two_clm, 4.8))
gs = gridspec.GridSpec(figure=fig1,
                       nrows=4,
                       ncols=6,
                       height_ratios=[1,
                                      0.04,
                                      1.5,
                                      0.02],
                       width_ratios=[1, 1, 1, 1, 1, 1])

od_no = ['a', 'b', 'c', 'd', 'e', 'f', 'g']

set_pos = [[0.01, 0.71, 0.32, 0.25],     # a 0
           [0.34, 0.71, 0.32, 0.25],  # b 1
           [0.67, 0.71, 0.32, 0.25],  # c 2
           [0.01, 0.07, 0.48, 0.5],      # d 3
           [0.51, 0.07, 0.48, 0.5],      # e 4
           [0.015, 0.66, 0.31, 0.015],     # cb Tair 5
           [0.345, 0.66, 0.31, 0.015],     # cb Q 6
           [0.675, 0.66, 0.31, 0.015],      # cb RH 7
           [0.05, 0.09, 0.9, 0.02],  # cb HSI 8
           ]

pt = 1  # 0: 95th, 1: 99th
mm = 0  # 0: max,  1: mean

# Obtain the HSI vote
sign_dec_trend = dec_trend.copy()
sign_dec_trend[dec_trend > 0] = 1
sign_dec_trend[dec_trend < 0] = -1

# Check the significance of trends
sigt_dec_trend = lin_sig.copy()
sigt_dec_trend[lin_sig < 0.05] = 1
sigt_dec_trend[lin_sig >= 0.05] = 0

# daily mean temperature: 99th
cmap_t = ncl.newcmp('t10_r')
bounds_t = np.linspace(-1, 1, 11)
norm_t = mpl.colors.BoundaryNorm(bounds_t, cmap_t.N)

ax = fig1.add_subplot(gs[0, 0: 2], projection=ccrs.Robinson())
ax.set_position(set_pos[0])

to_show = dec_trend[0, pt, mm, :, :]
sig = lin_sig[0, pt, mm, :, :].copy()
plot_trend(to_show, sig, cmap_t, norm_t, r'T$_{\mathrm{air}}$ trend', 'a', s=0.5, is_hsi=0)

# colorbar for T
ax = fig1.add_subplot(gs[1, 0: 2])
ax.set_position(set_pos[5])
cb = colorbar.ColorbarBase(ax,
                           cmap=cmap_t,
                           norm=norm_t,
                           orientation='horizontal',
                           extend='both',
                           ticks=[-0.8, -0.4, 0, 0.4, 0.8])
cb.ax.set_xticklabels([-0.8, -0.4, 0, 0.4, 0.8], fontsize=8)
cb.set_label(label=r'$^{\mathrm{o}}$C dec$^{-1}$', size=9)


# specific humidity: 99th
cmap_sh = ncl.newcmp('p10')
bounds_sh = np.linspace(-1, 1, 11)
norm_sh = mpl.colors.BoundaryNorm(bounds_sh, cmap_sh.N)

ax = fig1.add_subplot(gs[0, 2: 4], projection=ccrs.Robinson())
ax.set_position(set_pos[1])

to_show = dec_trend[7, pt, mm, :, :] * 1000  # * 1000: change kg/kg to g/kg
sig = lin_sig[7, pt, mm, :, :].copy()
plot_trend(to_show, sig, cmap_sh, norm_sh, 'Q trend', 'b', s=0.5, is_hsi=0)

# colorbar for Q
ax = fig1.add_subplot(gs[1, 2: 4])
ax.set_position(set_pos[6])
cb = colorbar.ColorbarBase(ax,
                           cmap=cmap_sh,
                           norm=norm_sh,
                           orientation='horizontal',
                           extend='both',
                           ticks=[-0.8, -0.4, 0, 0.4, 0.8])
cb.ax.set_xticklabels([-0.8, -0.4, 0, 0.4, 0.8], fontsize=8)
cb.set_label(label='g/kg dec$^{-1}$', size=9)


# relative humidity: 99th
cmap_rh = ncl.newcmp('p10')
bounds_rh = np.linspace(-5, 5, 11)
norm_rh = mpl.colors.BoundaryNorm(bounds_rh, cmap_rh.N)

ax = fig1.add_subplot(gs[0, 4: 6], projection=ccrs.Robinson())
ax.set_position(set_pos[2])

to_show = dec_trend[8, pt, mm, :, :] * 100  # * 100: change 0-1 to 0-100%
sig = lin_sig[8, pt, mm, :, :].copy()
plot_trend(to_show, sig, cmap_rh, norm_rh, 'RH trend', 'c', s=0.5, is_hsi=0)

# colorbar for RH
ax = fig1.add_subplot(gs[1, 4: 6])
ax.set_position(set_pos[7])
cb = colorbar.ColorbarBase(ax,
                           cmap=cmap_rh,
                           norm=norm_rh,
                           orientation='horizontal',
                           extend='both',
                           ticks=[-4, -2, 0, 2, 4])
cb.ax.set_xticklabels([-4, -2, 0, 2, 4], fontsize=8)
cb.set_label(label='% dec$^{-1}$', size=9)

# the consistency among 6 HSIs (mean)
cmap_cst = plt.cm.get_cmap('RdYlBu_r', 7)
bounds_cst = np.linspace(-7, 7, 8)
norm_cst = mpl.colors.BoundaryNorm(bounds_cst, cmap_cst.N)

ax = fig1.add_subplot(gs[2, 0: 3], projection=ccrs.Robinson())
ax.set_position(set_pos[3])

to_show = np.nansum(sign_dec_trend[1:7, pt, 1, :, :], axis=0)  # sum the HSI vote
sig_s = np.nansum(sigt_dec_trend[1:7, pt, 1, :, :], axis=0)
sig_s[sig_s > 0.5] = 1  # at least trend for one HSI reaches significance
plot_trend(to_show, np.NaN, cmap_cst, norm_cst, r'HSI trend (day$_{\mathrm{mean}}$)', 'd', s=0.5, is_hsi=1)

# the consistency among 6 HSIs (max)
ax = fig1.add_subplot(gs[2, 3: 6], projection=ccrs.Robinson())
ax.set_position(set_pos[4])

to_show = np.nansum(sign_dec_trend[1:7, pt, 0, :, :], axis=0)
sig_s = np.nansum(sigt_dec_trend[1:7, pt, 0, :, :], axis=0)
sig_s[sig_s > 0.5] = 1
plot_trend(to_show, np.NaN, cmap_cst, norm_cst, r'HSI trend (day$_{\mathrm{max}}$)', 'e', s=0.5, is_hsi=1)

# colorbar for consistency
ax = fig1.add_subplot(gs[3, 0: 6])
ax.set_position(set_pos[8])
cb = colorbar.ColorbarBase(ax,
                           cmap=cmap_cst,
                           norm=norm_cst,
                           orientation='horizontal',
                           extend='neither',
                           ticks=[-6, -4, -2, 0, 2, 4, 6])
cb.ax.set_xticklabels([-6, -4, -2, 0, 2, 4, 6], fontsize=9)
cb.ax.set_xticklabels(['-6', '-4', '-2', '0', '2', '4', '6'])
cb.set_label(label='HSI Vote', size=9)

plt.show()
# fig1.savefig('Fig_1.jpg',
#              dpi=1200,
#              format='jpg')

print('All Finished.')

