import numpy as np

def SineFit(t: np.ndarray, a: float, b: float, c: float) -> np.ndarray:
    '''
    Simple sine function to fit the data
    '''
    return   a* np.sin(b * t + c)

def InversePeriod(x: float, a: float) -> float:
    '''
    Simple relation between number of spots and period
    '''
    return a/(x)


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    import starry
    import astropy.units as u
    from scipy.optimize import curve_fit
    from scipy.interpolate import interp1d 
    from my_plot import TEX_FONTS, LATEX_PLOT, PlotStyle 
    from constants import *

    starry.config.lazy = False
    starry.config.quiet = True

    beta_mass = 11.9    # in u.Mjup
    beta_rad = 1.65     # in u.Rjup
    theta = np.linspace(-180, 180, 1000)                    # in degrees
    time_ros = np.linspace(0,BETA_DAY,len(theta))/(60*60)   # in hours

    # Generate the data

    n=7
    results=np.empty((6*n,4))
    for k in range(6):
        latit=k*15
        for i in range(n):
            A = starry.Primary(
                starry.Map(ydeg=10, udeg=2, rv=True, amp=1, veq=BETA_VEQ, alpha=0),
                r=beta_rad,
                m=beta_mass,
                length_unit=u.Rjup,
                mass_unit=u.Mjup)
            map_beta=A.map
            map_beta[1] = 0.5
            map_beta[2] = 0.25
            number_of_spots=4+i
            for j in range(number_of_spots):
                map_beta.spot(contrast=.95 , radius=15 ,lat=latit, lon=j*(360/number_of_spots))
            rv_ros = map_beta.rv(theta=theta)
            
            val=0
            for j in range(len(time_ros)):
                if abs(rv_ros[j+1])-abs(rv_ros[j])<0:
                    val=j
                    break

            fit_freq=2*np.pi/(4*time_ros[val])
            fit_ampl=np.amax(rv_ros)
            
            params, params_covariance = curve_fit(SineFit, time_ros, rv_ros , p0=[fit_ampl,fit_freq,0])
            results[k*n+i,0]=latit
            results[k*n+i,1]=number_of_spots
            results[k*n+i,2]=params[1]
            results[k*n+i,3]=abs(params[0])

    results=results.T
    latitudes=results[0]
    n_spots=results[1]
    frequencies=results[2]
    periods=2*np.pi/frequencies
    amplitudes=results[3]

    nsp=np.unique(n_spots)      # number of spots
    ltd=np.unique(latitudes)    # number of latitudes

    # Analysis of Periodicity

    # Create a mask for each number of spots
    mask4=n_spots==4
    mask5=n_spots==5
    mask6=n_spots==6
    mask7=n_spots==7
    mask8=n_spots==8
    mask9=n_spots==9
    mask10=n_spots==10

    with PlotStyle('science', TEX_FONTS, LATEX_PLOT):
        plt.imshow(periods.reshape(len(ltd),len(nsp)).T)#,extent=[0,75,4,10])
        x = np.arange(0,76,15) # the grid to which your data corresponds
        nx = x.shape[0]
        no_labels = 5 # how many labels to see on axis x
        step_x = int(nx / (no_labels - 1)) # step between consecutive labels
        x_positions = np.arange(0,nx,step_x) # pixel count at label position
        x_labels = x[::step_x] # labels you want to see
        plt.xticks(x_positions, x_labels)

        y = np.arange(4,11,1) # the grid to which your data corresponds
        ny = y.shape[0]
        no_labels = 7 # how many labels to see on axis x
        step_y = int(ny / (no_labels - 1)) # step between consecutive labels
        y_positions = np.arange(0,ny,step_y) # pixel count at label position
        y_labels = y[::step_x] # labels you want to see
        plt.yticks(y_positions, y_labels)
        plt.ylabel('Number of Spots')
        plt.xlabel(r'Latitude (degrees)')
        cb = plt.colorbar()
        cb.set_label('Period (days)')
        plt.savefig('plots_rossiter/periods.png',bbox_inches='tight',dpi=300)


    with PlotStyle('science', TEX_FONTS, LATEX_PLOT):
        plt.scatter(latitudes[mask4],periods[mask4],label='4 Spots')
        plt.scatter(latitudes[mask5],periods[mask5],label='5 Spots')
        plt.scatter(latitudes[mask6],periods[mask6],label='6 Spots')
        plt.scatter(latitudes[mask7],periods[mask7],label='7 Spots')
        plt.scatter(latitudes[mask8],periods[mask8],label='8 Spots')
        plt.scatter(latitudes[mask9],periods[mask9],label='9 Spots')
        plt.scatter(latitudes[mask10],periods[mask10],label='10 Spots')
        plt.ylabel('Period (hours)')
        plt.xlabel(r'Latitude (degrees)')
        legend=plt.legend(loc=1,bbox_to_anchor=(1.325,0.85),frameon=True)
        plt.savefig('plots_rossiter/period-lat.png',bbox_inches='tight',dpi=300)

    params, params_covariance = curve_fit(InversePeriod, n_spots, periods)

    with PlotStyle('science', TEX_FONTS, LATEX_PLOT):
        x = np.linspace(4, 20, 100)
        plt.scatter(n_spots,periods,s=20,label='Simulated Data',c='black')
        plt.plot(x, InversePeriod(x, params[0]), color='r',label='Fitted function')
        plt.ylabel('Period (hours)')
        plt.xlabel('Number of spots')
        plt.legend()
        plt.savefig('plots_rossiter/periods_fit.png',bbox_inches='tight',dpi=300)

    print(f'Fitted function: Period = {params[0]:.2f}/x')
    print(f'Error: {np.sqrt(params_covariance[0,0]):.2f}')

    # Analysis of Amplitudes

    plt.figure(figsize=(5,5))
    plt.imshow(amplitudes.reshape(len(ltd),len(nsp)).T,norm=SymLogNorm(linthresh=.0001))#,extent=[0,75,4,10])
    plt.gca().invert_yaxis()

    x = np.arange(0,76,15) # the grid to which your data corresponds
    nx = x.shape[0]
    no_labels = 5 # how many labels to see on axis x
    step_x = int(nx / (no_labels - 1)) # step between consecutive labels
    x_positions = np.arange(0,nx,step_x) # pixel count at label position
    x_labels = x[::step_x] # labels you want to see
    plt.xticks(x_positions, x_labels)

    y = np.arange(4,11,1) # the grid to which your data corresponds
    ny = y.shape[0]
    no_labels = 7 # how many labels to see on axis x
    step_y = int(ny / (no_labels - 1)) # step between consecutive labels
    y_positions = np.arange(0,ny,step_y) # pixel count at label position
    y_labels = y[::step_x] # labels you want to see
    plt.yticks(y_positions, y_labels)
    plt.ylabel('Number of Spots')
    plt.xlabel('Latitude')
    cb = plt.colorbar()
    cb.set_label('Amplitude (m/s)')
    plt.savefig('plots_rossiter/amplitudes.png',bbox_inches='tight',dpi=300)

    with PlotStyle('science', TEX_FONTS, LATEX_PLOT):
        plt.scatter(latitudes[mask4],amplitudes[mask4],label='4 spots')
        plt.scatter(latitudes[mask5],amplitudes[mask5],label='5 spots')
        plt.scatter(latitudes[mask6],amplitudes[mask6],label='6 spots')
        plt.scatter(latitudes[mask7],amplitudes[mask7],label='7 spots')
        plt.scatter(latitudes[mask8],amplitudes[mask8],label='8 spots')
        plt.scatter(latitudes[mask9],amplitudes[mask9],label='9 spots')
        plt.scatter(latitudes[mask10],amplitudes[mask10],label='10 spots')
        plt.ylabel('Amplitude')
        plt.xlabel('Latitude')
        plt.legend()
        plt.savefig('plots_rossiter/ampl_lat.png',bbox_inches='tight',dpi=300)

    size=(2*np.pi*BETA_RAD*np.cos(np.deg2rad(latitudes))*(15/360))/JUP_RAD #in u.Jup

    with PlotStyle('science', TEX_FONTS, LATEX_PLOT):
        plt.scatter(size[mask4],amplitudes[mask4],label='4 spots')
        plt.scatter(size[mask5],amplitudes[mask5],label='5 spots')
        plt.scatter(size[mask6],amplitudes[mask6],label='6 spots')
        plt.scatter(size[mask7],amplitudes[mask7],label='7 spots')
        plt.scatter(size[mask8],amplitudes[mask8],label='8 spots')
        plt.scatter(size[mask9],amplitudes[mask9],label='9 spots')
        plt.scatter(size[mask10],amplitudes[mask10],label='10 spots')
        plt.ylabel('Amplitude (m/s)')
        plt.xlabel('Size of Spots (Jupiter radii)')
        plt.legend()
        plt.savefig('plots_rossiter/ampl_size.png',bbox_inches='tight',dpi=300)

    # Interpolation of the data
    fun_1 = interp1d(size[mask10], amplitudes[mask10], fill_value='extrapolate',kind='linear')
    fun_2 = interp1d(size[mask10], amplitudes[mask10], fill_value='extrapolate',kind='quadratic')
    fun_3 = interp1d(size[mask10], amplitudes[mask10], fill_value='extrapolate',kind='cubic')

    plt.figure(figsize=(12, 5))
    plt.title('All 10 Spots cases')
    plt.scatter(size[mask10],amplitudes[mask10])
    x = np.linspace(0.1, .6, 100)
    plt.plot(x,fun_1(x))
    plt.plot(x,fun_2(x))
    plt.plot(x,fun_3(x))
    plt.ylabel('Amplitude [m/s]')
    plt.xlabel('Total Size of Spots [Jupiter radius]')
    plt.savefig('plots_rossiter/ampl_size_fit.png',bbox_inches='tight',dpi=300)

    lat_0=latitudes==0
    lat_15=latitudes==15
    lat_30=latitudes==30
    lat_45=latitudes==45
    lat_60=latitudes==60
    lat_75=latitudes==75

    with PlotStyle('science', TEX_FONTS, LATEX_PLOT):
        plt.scatter(n_spots[lat_0],amplitudes[lat_0],label='0 latitude')
        plt.scatter(n_spots[lat_15],amplitudes[lat_15],label='15 latitude')
        plt.scatter(n_spots[lat_30],amplitudes[lat_30],label='30 latitude')
        plt.scatter(n_spots[lat_45],amplitudes[lat_45],label='45 latitude')
        plt.scatter(n_spots[lat_60],amplitudes[lat_60],label='60 latitude')
        plt.scatter(n_spots[lat_75],amplitudes[lat_75],label='75 latitude')
        plt.ylabel('Amplitude')
        plt.xlabel('Number of Spots')
        plt.legend()
        plt.savefig('plots_rossiter/ampl_spots.png',bbox_inches='tight',dpi=300)

    # Interpolation of the data
    f_1 = interp1d(n_spots[lat_0],amplitudes[lat_0], fill_value='extrapolate',kind='linear')
    f_2 = interp1d(n_spots[lat_0],amplitudes[lat_0], fill_value='extrapolate',kind='quadratic')
    f_3 = interp1d(n_spots[lat_0],amplitudes[lat_0], fill_value='extrapolate',kind='cubic')

    plt.figure(figsize=(12, 5))
    plt.title('All 0 latitude cases')
    plt.scatter(n_spots[lat_0],amplitudes[lat_0],c='black',label='Data')
    x = np.linspace(4, 30, 100)
    plt.plot(x,f_1(x),label='Linear')
    plt.plot(x,f_2(x),label='Quadratic')
    plt.plot(x,f_3(x),label='Cubic')
    plt.legend()
    plt.ylim(-100,200)
    plt.ylabel('Amplitude [m/s]')
    plt.xlabel('Number of Spots')
    plt.savefig('plots_rossiter/ampl_spots_fit.png',bbox_inches='tight',dpi=300)