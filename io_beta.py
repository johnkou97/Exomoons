if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    import exoplanet as xo
    from my_plot import TEX_FONTS, LATEX_PLOT, PlotStyle
    from constants import *
    
    t = np.linspace(0, 20, 5000)    # 20 days with 5000 points

    # Define the orbit for Io orbiting Jupiter
    io_orbit = xo.orbits.KeplerianOrbit(
        period=IO_PERIOD_DAY,  # All times are in days
        incl=(1/2) * np.pi,  # All angles are in radians
        m_planet=IO_MASS/SOLAR_MASS,  # All masses and distances are in Solar units
        m_star=JUP_MASS/SOLAR_MASS,
        r_star=JUP_RAD/SOLAR_RAD,
        t0=0
    )

    rv_io = io_orbit.get_radial_velocity(t)

    # Define the orbit for Earth-mass moon around beta Pic b with P=10days
    beta_orbit = xo.orbits.KeplerianOrbit(
        period=10,  
        incl=(1/2) * np.pi,  
        m_planet=EAR_MASS/SOLAR_MASS, 
        m_star=BETA_MASS/SOLAR_MASS,
        r_star=BETA_RAD/SOLAR_RAD,
        t0=1.5
    )

    rv_beta = beta_orbit.get_radial_velocity(t)

    with PlotStyle('science', TEX_FONTS, LATEX_PLOT):
        plt.plot(t, rv_beta.eval(), label=r'Mass earth moon around $\beta$ Pic b',c='seagreen')
        plt.plot(t, rv_io.eval(), label=r'Io around Jupiter',c='crimson')
        plt.xlabel("Time (days)")
        plt.ylabel("Radial Velocity (m/s)")
        plt.xlim(t[0], t[-1])
        plt.legend()
        plt.xticks(np.arange(min(t),max(t)+1,1))
        plt.yticks(np.arange(-int(max(abs(rv_beta.eval())))-1,int(max(abs(rv_beta.eval())))+2,1))
        plt.grid()
        plt.savefig('plots/Io_betaPicb.png',bbox_inches='tight',dpi=300)