# Project: Regional Variation in the Role of Humidity on City-level Heat-Related Mortality
# Author: Dr.GUO Qiang, The University of Tokyo
# Contact: qiang@rainbow.iis.u-tokyo.ac.jp
# This script is developed and tested in Python 3.7.12 and Linux platform, and can be run also in Windows with
# proper environment settings.
# Description:
# This script is the to plot the difference of peak time (Figure 2 in the manuscript).


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
from matplotlib import colorbar
from matplotlib import colors
from mpl_toolkits.basemap import Basemap, maskoceans
from scipy.stats import gaussian_kde
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Self Defined Packages
import newcolors as ncl


mpl.use('Qt5Agg')
plt.close('all')
one_clm = 3.46
two_clm = 7.08

hsi_mcc_tick = [r'T$_{\mathrm{air}}$',
                r'T$_{\mathrm{w}}$', r'T$_{\mathrm{S}}$', r'T$_{\mathrm{WBG}}$', r'T$_{\mathrm{sWBG}}$',
                r'H$_{\mathrm{x}}$', 'APT', 'UTCI', 'HI']
hsi_grid_tick = [r'T$_{\mathrm{air}}$',
                 r'T$_{\mathrm{w}}$', r'T$_{\mathrm{sWBG}}$',
                 r'H$_{\mathrm{x}}$', 'APT', 'UTCI', 'HI']

# Load some pre-calculated data for plotting the figure
with np.load('fig2_peak_time.npz', allow_pickle=True) as file:
    grid_ptime = file['grid_ptime']  # peak time for land surface grid
    lon_mesh = file['lon_mesh']  # lon for land surface grid
    lat_mesh = file['lat_mesh']  # lat for land surface grid
    htst_city = file['htst_city']  # hottest 10 days for 4 cities
    city_lon = file['city_lon']  # lon for 4 cities
    city_lat = file['city_lat']  # lat for 4 cities
    city_name = file['city_name']  # name of 4 cities


########################################################################################################################
# function to mask ocean
def maskoc(to_show):

    to_show_mask = maskoceans(lon_mesh, lat_mesh, to_show, inlands=True, resolution='l')
    return to_show_mask


# function to plot density (Gaussian Kernel Density)
def plot_density(city_num, htst_10, cmap_ovlp, norm_ovlp, od):
    for i in range(9):
        if i == 0:
            x = np.ones(400) * (i + 1)  # hottest 10 days for each year, totalling 400 days
            y = htst_10[i, :, :, city_num].reshape(400)
            z = gaussian_kde(y)(y)
            idx = z.argsort()
            y, z = y[idx], z[idx]
        else:
            x = np.append(x, np.ones(400) * (i + 1))
            y_temp = htst_10[i, :, :, city_num].reshape(400)
            z_temp = gaussian_kde(y_temp)(y_temp)
            idx = z_temp.argsort()
            y_temp, z_temp = y_temp[idx], z_temp[idx]
            y = np.append(y, y_temp)
            z = np.append(z, z_temp)

    plt.scatter(y, x, c=z, s=20, cmap=cmap_ovlp, norm=norm_ovlp)
    plt.axis([0, 366, 0.5, 9.5])
    plt.yticks(np.linspace(1, 9, 9), hsi_mcc_tick, fontsize=8)
    plt.title('(' + str(od) + ') ' + city_name[city_num], fontsize=9)
    plt.gca().invert_yaxis()


# Start to plot
# colormap for subplot a
cmap_pat = plt.cm.get_cmap('plasma', 73)
bounds_pat = np.linspace(0, 365, 74)
norm_pat = mpl.colors.BoundaryNorm(bounds_pat, cmap_pat.N)

# colormap for subplot b-g
cmap_dif = ncl.newcmp('t11_r')
bounds_dif = np.linspace(-77, 77, 12, dtype='int')
norm_dif = mpl.colors.BoundaryNorm(bounds_dif, cmap_dif.N)

# colormap for subplot h-k
cmap_ovlp = plt.cm.get_cmap('Spectral_r')
norm_ovlp = colors.Normalize(vmin=0.000,
                             vmax=0.025)

fig1 = plt.figure(1, figsize=(two_clm, 7))
gs = gridspec.GridSpec(figure=fig1,
                       nrows=5,
                       ncols=5,
                       height_ratios=[1,
                                      1,
                                      1,
                                      1,
                                      0.07],
                       width_ratios=[1, 0.05, 1, 0.05, 0.8])

od_no = ['a', 'b', 'c', 'd', 'e', 'f', 'g']

set_pos = [[0.01, 0.72, 0.5, 0.25],  # a 0
           [0.01, 0.48, 0.31, 0.23],  # b 1
           [0.34, 0.48, 0.31, 0.23],  # c 2
           [0.01, 0.275, 0.31, 0.23],  # d 3
           [0.34, 0.275, 0.31, 0.23],  # e 4
           [0.01, 0.07, 0.31, 0.23],  # f 5
           [0.34, 0.07, 0.31, 0.23],  # g 6
           [0.72, 0.795, 0.25, 0.175],  # h 7
           [0.72, 0.58, 0.25, 0.175],  # i 8
           [0.72, 0.365, 0.25, 0.175],  # j 9
           [0.72, 0.15, 0.25, 0.175],  # k 10
           [0.53, 0.72, 0.015, 0.25],  # cb1 11
           [0.01, 0.060, 0.6, 0.015],  # cb2 12
           [0.72, 0.060, 0.25, 0.015]]  # cb3 13

# main colormap
for i in range(7):
    # subplot a
    if i == 0:
        cmp = cmap_pat
        bds = bounds_pat
        nrm = norm_pat
        to_show = maskoc(np.nanmean(grid_ptime[0, 1, :, :, :], axis=0))  # hsi, max/mean, year, lat, lon
        ax = fig1.add_subplot(gs[0, 0: 3], projection=ccrs.Robinson())

    # subplot b-g
    else:
        cmp = cmap_dif
        bds = bounds_dif
        nrm = norm_dif
        to_show = maskoc(np.nanmean(grid_ptime[i, 1, :, :, :] - grid_ptime[0, 1, :, :, :], axis=0))
        time_1 = np.nanmean((365 - grid_ptime[0, 1, :, :, :]) + grid_ptime[i, 1, :, :, :], axis=0)
        time_2 = np.nanmean((360 - grid_ptime[i, 1, :, :, :]) + grid_ptime[0, 1, :, :, :], axis=0)

    if (i == 1) | (i == 2):
        ax = fig1.add_subplot(gs[1, ((i - 1) % 2) * 2: ((i - 1) % 2) * 2 + 2], projection=ccrs.Robinson())
    elif (i == 3) | (i == 4):
        ax = fig1.add_subplot(gs[2, ((i - 1) % 2) * 2: ((i - 1) % 2) * 2 + 2], projection=ccrs.Robinson())
    elif (i == 5) | (i == 6):
        ax = fig1.add_subplot(gs[3, ((i - 1) % 2) * 2: ((i - 1) % 2) * 2 + 2], projection=ccrs.Robinson())

    ax.set_position(set_pos[i])

    # set the global map
    ax.set_global()
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), edgecolor=[0.2, 0.2, 0.2], linewidth=0.2, zorder=2)
    ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor=[0.8, 0.8, 0.8],
                   edgecolor=[0.2, 0.2, 0.2], linewidth=0.3, zorder=1)
    ax.add_feature(cfeature.OCEAN, facecolor=[0.7, 0.7, 0.7], zorder=4)

    # plot the shaded color
    plt.pcolormesh(lon_mesh[40: 600, :],
                   lat_mesh[40: 600, :],
                   to_show[40: 600, :],
                   cmap=cmp,
                   norm=nrm,
                   shading='auto',
                   zorder=2,
                   transform=ccrs.PlateCarree())
    plt.title(hsi_grid_tick[i], fontsize=10)

    # add the location of 4 cities
    if i == 0:
        plt.scatter(city_lon, city_lat,
                    s=8,
                    c='k',
                    zorder=5,
                    transform=ccrs.PlateCarree())
        plt.text(city_lon[0] * 0.9, city_lat[0] * 1.05, '(1)', fontsize=8, transform=ccrs.PlateCarree(), zorder=5)
        plt.text(city_lon[1] * 0.9, city_lat[1] * 1.05, '(2)', fontsize=8, transform=ccrs.PlateCarree(), zorder=5)
        plt.text(city_lon[2] * 0.9, city_lat[2] * 1.05, '(3)', fontsize=8, transform=ccrs.PlateCarree(), zorder=5)
        plt.text(city_lon[3] * 0.9, city_lat[3] * 1.05, '(4)', fontsize=8, transform=ccrs.PlateCarree(), zorder=5)
    plt.text(-1.65e+07, 8.1e+06, od_no[i], fontsize=9)


# colorbar for subplot a
ax = fig1.add_subplot(gs[0, 3])
ax.set_position(set_pos[11])
cb = colorbar.ColorbarBase(ax,
                           cmap=cmap_pat,
                           norm=norm_pat,
                           orientation='vertical',
                           extend='neither',
                           ticks=[0, 60, 120, 180, 240, 300, 360])
cb.ax.set_yticklabels(list([0, 60, 120, 180, 240, 300, 360]), fontsize=8)
cb.set_label(label='Day of year', size=9)

# colorbar for subplot a-g
ax = fig1.add_subplot(gs[4, 0: 4])
ax.set_position(set_pos[12])
cb = colorbar.ColorbarBase(ax,
                           cmap=cmap_dif,
                           norm=norm_dif,
                           orientation='horizontal',
                           extend='both',
                           ticks=np.linspace(-77, 77, 12, dtype='int'))
cb.ax.set_xticklabels(np.linspace(-77, 77, 12, dtype='int'), fontsize=8)
cb.set_label(label='Number of days', size=9)

# hot density, Austin
ax = fig1.add_subplot(gs[0, 4])
ax.set_position(set_pos[7])
plot_density(0, htst_city, cmap_ovlp, norm_ovlp, 1)
plt.xticks([])
plt.text(5, 0, 'h', fontsize=9)

# hot density, Brasilia
ax = fig1.add_subplot(gs[1, 4])
ax.set_position(set_pos[8])
plot_density(1, htst_city, cmap_ovlp, norm_ovlp, 2)
plt.xticks([])
plt.text(5, 0, 'i', fontsize=9)

# hot density, London
ax = fig1.add_subplot(gs[2, 4])
ax.set_position(set_pos[9])
plot_density(2, htst_city, cmap_ovlp, norm_ovlp, 3)
plt.xticks([])
plt.text(5, 0, 'j', fontsize=9)

# hot density, Bangkok
ax = fig1.add_subplot(gs[3, 4])
ax.set_position(set_pos[10])
plot_density(3, htst_city, cmap_ovlp, norm_ovlp, 4)
plt.xticks([0, 100, 200, 300], fontsize=8)
plt.xlabel('Day of year', fontsize=9)
plt.text(5, 0, 'k', fontsize=9)

# colorbar for subplot h-k
ax = fig1.add_subplot(gs[4, 4])
ax.set_position(set_pos[13])
cb = colorbar.ColorbarBase(ax,
                           orientation='horizontal',
                           cmap=cmap_ovlp,
                           norm=norm_ovlp,
                           ticks=[0.00, 0.01, 0.02])
cb.ax.set_xticklabels([0.00, 0.01, 0.02], fontsize=8)
cb.set_label(label='Occurrence frequency', size=9)

plt.show()
# fig1.savefig('Fig_2.jpg',
#              dpi=1200,
#              format='jpg')

print('All finished.')
