# Project: Global Apparent Temperature Change
# Project start time: 2021/08/10
# Author: Dr.GUO Qiang, The University of Tokyo
# Contact: qiang@rainbow.iis.u-tokyo.ac.jp
# Description:
# This script is the function summary of metrics calculation

import numpy as np
import math
from metpy.calc import dewpoint_from_relative_humidity
from metpy.calc import mixing_ratio_from_relative_humidity
from metpy.calc import vapor_pressure
from metpy.calc import heat_index
from metpy.units import units
from pythermalcomfort.models import utci_optimized
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
heatstress = importr('HeatStress')


def show_time():
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print("now = ", dt_string)


# 1 Calculate the wet-bulb temperature, method of Stull (2011)
def get_tw_stull(at, rh):
    # unit: Celsius degree, rh%

    tw_stull = at * np.arctan(0.151977 * (rh + 8.313659) ** 0.5) \
               + np.arctan(at + rh) \
               - np.arctan(rh - 1.676331) \
               + 0.00391838 * (rh ** 1.5) * np.arctan(0.023101 * rh) \
               - 4.686035

    return tw_stull


# 2 Calculate the Ts
def get_ts(tw, rh):

    rh = rh / 100
    ts = tw + 4.5 * (1 - rh ** 2)

    return ts


# 3 Calculate the wbgt, Liljegren et al., 2018
def get_wbgt(at, rh, ws, rad):
    # unit: Celsius degree, m/s, w/m2
    at_c = units.Quantity(at, "degC")
    rh = rh / 100
    dewt = dewpoint_from_relative_humidity(at_c, rh)

    utc_time = ["2021-03-21"]
    lon = 0.0
    lat = 0.0

    wbgt_cc = heatstress.wbgt_Liljegren(tas=robjects.FloatVector(np.array(at)),
                                        dewp=robjects.FloatVector(np.array(dewt.m)),
                                        wind=robjects.FloatVector(np.array(ws)),
                                        radiation=robjects.FloatVector(np.array(rad)),
                                        dates=robjects.StrVector(utc_time),
                                        lon=robjects.FloatVector(np.array(lon)),
                                        lat=robjects.FloatVector(np.array(lat)),
                                        noNAs=False,
                                        hour=False)
    wbgt = np.array(wbgt_cc)[0]
    return wbgt


# 4 Calculate the swbgt, method from Australian Bureau of Meteorology
def get_swbgt(at, rh, pres):

    at = units.Quantity(at, "degC")
    pres = units.Quantity(pres, "Pa")
    rh = rh / 100

    mx = mixing_ratio_from_relative_humidity(pres, at, rh)
    vp = vapor_pressure(pres, mx).m  # Pa

    swbgt = 0.567 * at.m + 0.393 * vp / 100 + 3.94

    return swbgt


# 5 Calculate the Humidex
def get_hx(at, rh, pres):

    at = units.Quantity(at, "degC")
    pres = units.Quantity(pres, "Pa")
    rh = rh / 100

    mx = mixing_ratio_from_relative_humidity(pres, at, rh)
    vp = vapor_pressure(pres, mx).m  # Pa

    hx = at.m + 5 / 9 * (vp / 100 - 10)

    return hx


# 6 Apparent Temperature, ABM
def get_apt(at, rh, ws, pres):

    at = units.Quantity(at, "degC")
    pres = units.Quantity(pres, "Pa")
    rh = rh / 100

    mx = mixing_ratio_from_relative_humidity(pres, at, rh)
    vp = vapor_pressure(pres, mx).m  # Pa

    apt = at.m + 0.33 * vp / 100 - 0.7 * ws - 4

    return apt


# 7 Calculate the Heat Index, method of Steadman (1979), Rothfusz (1990), Anderson (2013)
def get_hi(at, rh):

    rh = rh / 100
    at = units.Quantity(at, "degC")

    hi = heat_index(at,
                    rh,
                    mask_undefined=False).to(units.degC)

    return hi.m


# 8 Calculate the UTCI
def get_utci(tdb, tr, v, rh, pres):

    # def exponential(t_db):
    #     g = [
    #         -2836.5744,
    #         -6028.076559,
    #         19.54263612,
    #         -0.02737830188,
    #         0.000016261698,
    #         (7.0229056 * (10 ** (-10))),
    #         (-1.8680009 * (10 ** (-13))),
    #     ]
    #     tk = t_db + 273.15  # air temp in K
    #     es = 2.7150305 * math.log1p(tk)
    #     for count, i in enumerate(g):
    #         es = es + (i * (tk ** (count - 2)))
    #     es = math.exp(es) * 0.01  # convert Pa to hPa
    #     return es

    # if (
    #         (tdb.any() < -50.0)
    #         or (tdb.any() > 50.0)
    #         or ((tr - tdb).any() < -30.0)
    #         or ((tr - tdb).any() > 70.0)
    #         or (v.any() < 0.5)
    #         or (v.any() > 17)
    # ):
    #     raise ValueError(
    #         "The value you entered are outside the equation applicability limits"
    #     )

    at = units.Quantity(tdb, "degC")
    pres = units.Quantity(pres, "Pa")
    rh = rh / 100

    mx = mixing_ratio_from_relative_humidity(pres, at, rh)
    eh_pa = vapor_pressure(pres, mx).m  # Pa

    # eh_pa = exponential(tdb) * rh

    delta_t_tr = tr - tdb
    pa = eh_pa / 1000  # convert vapour pressure to kPa

    utci_approx = utci_optimized(tdb, v, delta_t_tr, pa)

    if (tdb < -50) | (tdb > 50) | (v > 30.3) | (pa > 5):
        utci_approx = np.NaN

    return utci_approx





