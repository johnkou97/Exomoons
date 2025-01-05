import astropy.constants as const
import numpy as np

# Usefull constants in S.I.

C = const.c.value
G = const.G.value
DAY = 24*60*60                          # in seconds 

SOLAR_MASS = const.M_sun.value
SOLAR_RAD = const.R_sun.value
ASTR_UNIT = const.au.value
JUP_MASS = const.M_jup.value
JUP_RAD = const.R_jup.value
EAR_MASS = const.M_earth.value
EAR_RAD = const.R_earth.value

IO_MASS = 89.3e21
IO_RAD = 1.8216e6
IO_PERIOD_DAY = 1.769                   # in days

BETA_MASS = 11.9*JUP_MASS               # taken from exoplanets.eu
BETA_RAD = 1.65*JUP_RAD
BETA_DAY = 8.1 * 60 * 60                # in seconds
BETA_VEQ = 2*np.pi*(BETA_RAD)/BETA_DAY  # in m/s