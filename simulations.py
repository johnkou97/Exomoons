import os
import exoplanet as xo
import numpy as np
import pymc3 as pm
import pymc3_ext as pmx
import arviz as az
import aesara_theano_fallback.tensor as tt
from constants import *

if not os.path.exists('traces'):
    os.makedirs('traces')

def SimulateOrbit(mass: int = 50, period: int = 10, error: int = 500) -> tuple:
    '''
    Simulate a dataset for a moon orbiting a planet with an injected error.
    '''
    # Define the orbit for Earth-mass moon around beta Pic b with P=10days
    beta_orbit = xo.orbits.KeplerianOrbit(
        period=period,  
        incl=(1/2) * np.pi,  
        m_planet=mass*EAR_MASS/SOLAR_MASS, 
        m_star=BETA_MASS/SOLAR_MASS,
        r_star=BETA_RAD/SOLAR_RAD,
        t0=1.5)

    # Create dataset
    random = np.random.default_rng(1234)
    t_plot = np.linspace(0, 140, 500)
    t = np.linspace(1,130,25)
    rv_err = error + np.zeros_like(t)
    rv_obs = beta_orbit.get_radial_velocity(t).eval() + np.sqrt(0.05**2 + rv_err**2) * random.normal(size=len(t))

    return t, rv_obs, rv_err, t_plot

def SimulateDetection(t: np.ndarray, rv_obs: np.ndarray, rv_err: np.ndarray, t_plot: np.ndarray, init_period: int = 10, name: str = 'detection') -> az.InferenceData:
    '''
    Simulate a detection of a moon orbiting a planet with Bayesian inference.
    '''
    with pm.Model():
        # Period, semi-amplitude, and eccentricity
        log_period = pm.Normal("log_period", mu=np.log(init_period), sigma=1.0) # Normal prior
        period = pm.Deterministic("period", tt.exp(log_period))                 # Deterministic variable
        log_semiamp = pm.Normal("log_semiamp", mu=np.log(750), sd=50)           # Normal prior
        semiamp = pm.Deterministic("semiamp", tt.exp(log_semiamp))              # Deterministic variable
        ecc = pm.Uniform("ecc", lower=0, upper=1)                               # Uniform prior

        # At low eccentricity, omega and the phase of periastron (phi) are
        # correlated so it can be best to fit in (omega ± phi) / 2
        plus = pmx.Angle("plus")                            # Angle prior
        minus = pmx.Angle("minus")                          # Angle prior
        phi = pm.Deterministic("phi", plus + minus)         # Deterministic variable
        omega = pm.Deterministic("omega", plus - minus)     # Deterministic variable

        # Jitter & the system mean velocity offset
        log_jitter = pm.Normal("log_jitter", mu=np.log(0.05), sd=5.0)   # Normal prior
        zero_point = pm.Normal("zero_point", mu=0, sd=10.0)             # Normal prior

        # Then we define the orbit
        tperi = pm.Deterministic("tperi", period * phi / (2 * np.pi))
        orbit = xo.orbits.KeplerianOrbit(period=period, t_periastron=tperi, ecc=ecc, omega=omega)

        # And define the RV model
        rv_model = zero_point + orbit.get_radial_velocity(t, K=semiamp)

        # Finally add in the observation model
        err = tt.sqrt(rv_err**2 + tt.exp(2 * log_jitter))
        pm.Normal("obs", mu=rv_model, sigma=rv_err, observed=rv_obs)

        # We'll also track the model just for plotting purposes
        pm.Deterministic("rv_plot", zero_point + orbit.get_radial_velocity(t_plot, K=semiamp))

        soln = pmx.optimize(vars=[plus, minus, ecc])
        soln = pmx.optimize(soln)
        trace = pmx.sample(
            tune=1000,
            draws=1000,
            cores=2,
            chains=2,
            start=soln,
            return_inferencedata=True)

    az.to_netcdf(trace, f'traces/{name}_trace.nc')
    az.summary(trace, var_names=["^(?!rv_plot).*"], filter_vars="regex").to_csv(f'traces/{name}_summary.csv')

    return trace

def PlotResults(trace: az.InferenceData, t: np.ndarray, rv_obs: np.ndarray, rv_err: np.ndarray, t_plot: np.ndarray, name: str = 'detection') -> None:
    '''
    Plot the results of the Bayesian inference.
    '''
    # Radial velocity data and model
    with PlotStyle('science', TEX_FONTS, LATEX_PLOT):
        plt.errorbar(t, rv_obs, yerr=rv_err, fmt=".k", label="data")
        plt.plot(t_plot, q50)
        plt.fill_between(t_plot, q16, q84, alpha=0.3, label="posterior")
        plt.xlim(t_plot[0], t_plot[-1])
        # plt.xlim(0, 50)   # Uncomment to zoom in
        plt.legend()    # ,bbox_to_anchor=(0.30,0.86))
        plt.xlabel(r'Time (days)')
        plt.ylabel(r'Radial Velocity (m/s)')
        plt.savefig(f'plots_detection/{name}_rv-t.png',bbox_inches='tight')

    # Corner plot
    with PlotStyle('science', TEX_FONTS):
        fig=corner(trace, var_names=['period','ecc'],labels=[r'period (days)',r'eccentricity'],**CORNER_KWARGS)
        fig.set_size_inches(LATEX_PLOT[0],LATEX_PLOT[1]*1.5)
        plt.savefig(f'plots_detection/{name}_corner.png',bbox_inches='tight')


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from corner import corner
    from my_plot import TEX_FONTS, CORNER_KWARGS, LATEX_PLOT, PlotStyle

    if not os.path.exists('plots_detection'):
        os.makedirs('plots_detection')

    error = 250                     # in m/s
    period = 10                     # in days
    masses = [60, 40, 20]           # in Earth masses
    for mass in masses:
        name = f'{error}_{period}_{mass}'
        t, rv_obs, rv_err, t_plot = SimulateOrbit(mass, period, error)
        trace = SimulateDetection(t, rv_obs, rv_err, t_plot, period, name)
        rv_plot = trace.posterior["rv_plot"].values
        q16, q50, q84 = np.percentile(rv_plot, [16, 50, 84], axis=(0, 1))
        PlotResults(trace, t, rv_obs, rv_err, t_plot, name)

    error = 250                     # in m/s
    period = 20                     # in days
    masses = [90, 70, 50]           # in Earth masses
    for mass in masses:
        name = f'{error}_{period}_{mass}'
        t, rv_obs, rv_err, t_plot = SimulateOrbit(mass, period, error)
        trace = SimulateDetection(t, rv_obs, rv_err, t_plot, period, name)
        rv_plot = trace.posterior["rv_plot"].values
        q16, q50, q84 = np.percentile(rv_plot, [16, 50, 84], axis=(0, 1))
        PlotResults(trace, t, rv_obs, rv_err, t_plot, name)

    error = 500                     # in m/s
    period = 10                     # in days
    masses = [80, 60, 40]           # in Earth masses
    for mass in masses:
        name = f'{error}_{period}_{mass}'
        t, rv_obs, rv_err, t_plot = SimulateOrbit(mass, period, error)
        trace = SimulateDetection(t, rv_obs, rv_err, t_plot, period, name)
        rv_plot = trace.posterior["rv_plot"].values
        q16, q50, q84 = np.percentile(rv_plot, [16, 50, 84], axis=(0, 1))
        PlotResults(trace, t, rv_obs, rv_err, t_plot, name)

    error = 500                     # in m/s
    period = 20                     # in days
    masses = [170, 150, 130, 110]   # in Earth masses
    for mass in masses:
        name = f'{error}_{period}_{mass}'
        t, rv_obs, rv_err, t_plot = SimulateOrbit(mass, period, error)
        trace = SimulateDetection(t, rv_obs, rv_err, t_plot, period, name)
        rv_plot = trace.posterior["rv_plot"].values
        q16, q50, q84 = np.percentile(rv_plot, [16, 50, 84], axis=(0, 1))
        PlotResults(trace, t, rv_obs, rv_err, t_plot, name)