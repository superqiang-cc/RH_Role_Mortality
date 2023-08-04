# Project: Regional Variation in the Role of Humidity on City-level Heat-Related Mortality
# Author: Dr.GUO Qiang, The University of Tokyo
# Contact: qiang@rainbow.iis.u-tokyo.ac.jp
# This script is developed and tested in Python 3.7.12 and Linux platform, and can be run also in Windows with
# proper environment settings.
# Description:
# This script is the to plot the heat-related deaths (Figure 6 in the manuscript)
# The four cities shown are: Miami (U.S.), Bristol (UK), Hochiminh (Vietnam), Taipei (Taiwan)


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl

mpl.use('Qt5Agg')
plt.close('all')
one_clm = 3.46
two_clm = 7.08

hsi_tick_humi_weight = [r'T$_{\mathrm{air}}$',
                        'HI', 'UTCI', 'APT', r'H$_{\mathrm{x}}$',
                        r'T$_{\mathrm{sWBG}}$', r'T$_{\mathrm{WBG}}$',
                        r'T$_{\mathrm{s}}$', r'T$_{\mathrm{w}}$']

bfi_4 = [6, 7, 6, 6]  # Twbg, Ts, Twbg, Twbg, the best fit indicator for 4 cities

# Load some pre-calculated data for plotting the figure
with np.load('fig6_heat_death.npz') as file:
    rr_365_4 = file['rr_365_4']  # The time series of RR for 365 days of 4 cities
    cen_8 = file['cen_8']  # The optimal value for Tair and BFI of 4 cities
    afabs_hot_4 = file['afabs_hot_4']  # The heat-related death estimated for 4 cities
    hd_365_4 = file['hd_365_4']  # The indicator values extracted for the warm season (warmest 6 consecutive months)
    date_365 = file['date_365']  # The date matrix (i.e., year, month, day) for a year
    warmsea_4 = file['warmsea_4']  # The warm season (warmest consecutive 6 months) for 4 cities
    mia_t = file['mia_t']  # The Tair-mortality risks association of Miami
    mia_h = file['mia_h']  # The BFI-mortality risks association of Miami
    bst_t = file['bst_t']  # The Tair-mortality risks association of Bristol
    bst_h = file['bst_h']  # The BFI-mortality risks association of Bristol
    hoc_t = file['hoc_t']  # The Tair-mortality risks association of Ho Chi Minh City
    hoc_h = file['hoc_h']  # The BFI-mortality risks association of Ho Chi Minh City
    tpi_t = file['tpi_t']  # The Tair-mortality risks association of Taipei
    tpi_h = file['tpi_h']  # The BFI-mortality risks association of Taipei


# #######################################################################################################################


# plot the RR variation
def plot_time_series(ax, rr_365, hd_365, af, bfi, date, wm, cen, name, od):

    x = date[:, 3]
    plt.title(name, fontsize=9)

    # remove the data out of the warm season
    rr_365[~np.in1d(date[:, 1], wm), :] = np.NaN
    hd_365[~np.in1d(date[:, 1], wm), :] = np.NaN

    # rr for tair and BFI
    l1 = ax.plot(x, rr_365[:, 0], color=tc, label=r'T$_{\mathrm{air}}$', lw=0.8)
    l2 = ax.plot(x, rr_365[:, 1], color=hc, label='BFI', lw=0.8)

    # indicate the hot days (indicator value > optimal value)
    ax.fill_between(x, 1.06, where=hd_365[:, 0] > cen[0], facecolor=tc, alpha=0.3)
    ax.fill_between(x, 1.06, where=hd_365[:, 1] > cen[1], facecolor=hc, alpha=0.3)

    # add the text of heat-related deaths
    ax.text(360, 1.05,
            r"T$_{\mathrm{air}}$: " + str(format(af[0, 0], '.2f')) + "%"
            + " (" + str(format(af[0, 2], '.2f')) + ' - ' + str(format(af[0, 1], '.2f')) + ")"
            , fontsize=8, c=tc, horizontalalignment='right')
    ax.text(360, 1.04,
            bfi + ": " + str(format(af[1, 0], '.2f')) + "%"
            + " (" + str(format(af[1, 2], '.2f')) + ' - ' + str(format(af[1, 1], '.2f')) + ")"
            , fontsize=8, c=hc, horizontalalignment='right')

    # add the ticks, etc.
    ax.axis([0, 365, 1, 1.06])
    plt.yticks(np.linspace(1.00, 1.06, 3), [1.00, 1.03, 1.06], fontsize=8)
    ax.set_ylabel('RR', fontsize=9)

    ax.text(0, 1.062, od, fontsize=9)

    if name == 'Taipei':
        plt.xticks([0, 50, 100, 150, 200, 250, 300, 350],
                   [0, 50, 100, 150, 200, 250, 300, 350],
                   fontsize=8)

        ax.set_xlabel('Day of year', fontsize=9)
    else:
        ax.set_xticks([])


def plot_rr_value(ax, t_rr, h_rr, name, od, cen, bfi):

    # plot the exposure-response curve
    ax.plot(t_rr[:, 0], t_rr[:, 1], c=tc, label=r'T$_{\mathrm{air}}')
    ax.plot(h_rr[:, 0], h_rr[:, 1], c=hc, label='BFI')

    # plot the uncertainty, the shaded area
    ax.fill_between(t_rr[:, 0], y1=t_rr[:, 2], y2=t_rr[:, 3], facecolor=tc, alpha=0.3)
    ax.fill_between(h_rr[:, 0], y1=h_rr[:, 2], y2=h_rr[:, 3], facecolor=hc, alpha=0.3)

    # plot the 95th quantile value of the indicator, the dashed line
    ax.axvline(np.quantile(t_rr[:, 0], 0.95), color=tc, linewidth=0.5, linestyle='--')
    ax.axvline(np.quantile(h_rr[:, 0], 0.95), color=hc, linewidth=0.5, linestyle='--')

    # adjust the axis settings, etc.
    id_value = np.hstack((t_rr[:, 0], h_rr[:, 0]))
    ax.axis([np.min(id_value), np.max(id_value), 1, 2.2])
    ax.set_ylabel('RR', fontsize=9)
    plt.yticks(np.linspace(1.00, 2.0, 3), [1.0, 1.5, 2.0], fontsize=8)
    plt.title(name, fontsize=9)
    ax.text(np.min(id_value), 2.24, od, fontsize=9)

    plt.xticks(np.round(np.linspace(np.min(id_value), np.max(id_value), 4)),
               np.round(np.linspace(np.min(id_value), np.max(id_value), 4)),
               fontsize=8)

    if name == 'Taipei':
        ax.set_xlabel('Indicator value', fontsize=9)

    # add the optimal value (minimum mortality risk) for both indicators
    ax.text(np.min(id_value) + np.ptp(id_value) * 0.03, 2.01,
            r"T$_{\mathrm{air}}$: " + str(format(cen[0], '.1f')) + r' $^{\mathrm{o}}$C', fontsize=8, c=tc)
    ax.text(np.min(id_value) + np.ptp(id_value) * 0.03, 1.81,
            bfi + ": " + str(format(cen[1], '.1f')) + r' $^{\mathrm{o}}$C', fontsize=8, c=hc)


# color for Tair and BFI
tc = 'k'
hc = 'r'

# Start to plot
fig1 = plt.figure(1, figsize=(two_clm, 6.5), constrained_layout=True)
gs = gridspec.GridSpec(figure=fig1,
                       nrows=4,
                       ncols=2,
                       height_ratios=[1,
                                      1,
                                      1,
                                      1],
                       width_ratios=[1, 1])

# Miami rr, subplot a
ax = fig1.add_subplot(gs[0, 0])
plot_rr_value(ax, mia_t, mia_h, 'Miami', 'a', cen_8[0: 2], hsi_tick_humi_weight[bfi_4[0]])

# Miami time series, subplot b
ax = fig1.add_subplot(gs[0, 1])
plot_time_series(ax, rr_365_4[0, :, :], hd_365_4[0, :, :], afabs_hot_4[:, 0, :],
                 hsi_tick_humi_weight[bfi_4[0]], date_365, warmsea_4[0, :],
                 cen_8[0: 2], 'Miami', 'b')

# Bristol rr, subplot c
ax = fig1.add_subplot(gs[1, 0])
plot_rr_value(ax, bst_t, bst_h, 'Bristol', 'c', cen_8[2: 4], hsi_tick_humi_weight[bfi_4[1]])

# Bristol time series, subplot d
ax = fig1.add_subplot(gs[1, 1])
plot_time_series(ax, rr_365_4[1, :, :], hd_365_4[1, :, :], afabs_hot_4[:, 1, :],
                 hsi_tick_humi_weight[bfi_4[1]], date_365, warmsea_4[1, :],
                 cen_8[2: 4], 'Bristol', 'd')

# Ho chi Minh rr, subplot e
ax = fig1.add_subplot(gs[2, 0])
plot_rr_value(ax, hoc_t, hoc_h, 'Ho Chi Minh City', 'e', cen_8[4: 6], hsi_tick_humi_weight[bfi_4[2]])

# Ho chi Minh time series, subplot f
ax = fig1.add_subplot(gs[2, 1])
plot_time_series(ax, rr_365_4[2, :, :], hd_365_4[2, :, :], afabs_hot_4[:, 2, :],
                 hsi_tick_humi_weight[bfi_4[2]], date_365, warmsea_4[2, :],
                 cen_8[4: 6], 'Ho Chi Minh City', 'f')

# Taipei rr, subplot g
ax = fig1.add_subplot(gs[3, 0])
plot_rr_value(ax, tpi_t, tpi_h, 'Taipei', 'g', cen_8[6: 8], hsi_tick_humi_weight[bfi_4[3]])

# Taipei time series, subplot h
ax = fig1.add_subplot(gs[3, 1])
plot_time_series(ax, rr_365_4[3, :, :], hd_365_4[3, :, :], afabs_hot_4[:, 3, :],
                 hsi_tick_humi_weight[bfi_4[3]], date_365, warmsea_4[3, :],
                 cen_8[6: 8], 'Taipei', 'h')

plt.show()

# fig2.savefig('Fig_6.jpg',
#              format='jpg',
#              dpi=1200)

print('All Finished.')


