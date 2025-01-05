from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from constants import *
from my_plot import TEX_FONTS, LATEX_PLOT, PlotStyle

def VelocityMoon(mass: np.ndarray, period: np.ndarray) -> np.ndarray:
	'''
	Compute the radial velocity of a moon around a planet given the mass and period of the moon.
	'''
	vel = np.zeros((len(mass),len(period)))
	for i,m in enumerate(tqdm(mass)):
		for j,p in enumerate(period):
			vel[i,j] = ((2*np.pi*G)/(BETA_MASS**2 * p))**(1/3) * m
	return vel

def RochePeriod(Mmoon: np.ndarray) -> np.ndarray:
	'''
	Compute the Roche period of a moon around a planet given the mass of the moon
	Created by Elina Kleisioti
	'''
	Rmoon = 129362 *10**3* (Mmoon/JUP_MASS)**0.55 # if ms<0.39 Mjup
	Roche_radius = Rmoon * (2*BETA_MASS/Mmoon)**(1/3)
	T = (Roche_radius**3*4*np.pi**2/(G*(BETA_MASS)))**(1/2)/(60*60*24) 
	return T 


if __name__ == '__main__':
	#set resolution
	size_mass=1000
	size_period=10000

	moon_mass=np.linspace(.1*(EAR_MASS/SOLAR_MASS),160*(EAR_MASS/SOLAR_MASS),size_mass)
	moon_period=np.linspace(.1,1e2,size_period)

	rv_max = VelocityMoon(moon_mass*SOLAR_MASS,moon_period*(24*60*60))
	t_roch = RochePeriod(moon_mass*SOLAR_MASS)
	
	hill_sphere = 3909.08	# Provided by Elina Kleisioti

	x_pos = [4.6,10,15,20]
	y_pos = [40,50,80,100]

	Kepler1708 = [4.6, 37] 
	Kepler1625 = [22, 19]

	x = moon_period
	y = moon_mass*SOLAR_MASS/EAR_MASS

	with PlotStyle('science', TEX_FONTS, LATEX_PLOT):
		plt.contourf(x,y,rv_max)
		cb = plt.colorbar()
		cb.set_label(label='Peak Radial Velocity (m/s)')

		CS=plt.contour(x,y,rv_max,levels=[500],linestyles='--',colors='white')
		positions = [(1,35)]
		label = plt.clabel(CS,fmt='500m/s', inline=False,rightside_up=False,colors='C3',fontsize='15',manual=positions)
		for l in label:
			l.set_rotation(20)
		
		plt.fill_between(t_roch,moon_mass[-1]*SOLAR_MASS/EAR_MASS,moon_mass*SOLAR_MASS/EAR_MASS,color='grey')
		plt.text(.15,70,'Roche limit',fontsize='15',rotation=90)
		
		plt.scatter(Kepler1708[0],Kepler1708[1],label='Kepler-1708 b-i',marker='o',color='snow')
		plt.scatter(Kepler1625[0],Kepler1625[1],label='Kepler-1625 b-i',marker='o',color='gold')
		plt.scatter(x_pos,y_pos,label='RV detection limit',marker='+',color='cyan')
		
		plt.legend(frameon=True,loc='upper center')
		plt.xscale('log')
		plt.xlim(moon_period[0],moon_period[-1])
		plt.ylim(moon_mass[0]*SOLAR_MASS/EAR_MASS,moon_mass[-1]*SOLAR_MASS/EAR_MASS)
		plt.xlabel('Orbital Period (days)')
		plt.ylabel(r'Moon Mass ($M_\oplus$)')
		plt.tick_params(which='both',color='red')
		plt.savefig(f'plots/contour_no_hill.png',dpi=300,bbox_inches='tight')

	# increase moon_period limit to include hill sphere
	size_period_extension = 1000
	moon_period_extension = np.linspace(1e2, 7e3, size_period_extension)
	rv_max_extension = VelocityMoon(moon_mass*SOLAR_MASS,moon_period_extension*(24*60*60))
	rv_max = np.concatenate((rv_max, rv_max_extension), axis=1)
	moon_period = np.concatenate((moon_period, moon_period_extension))
	x = moon_period  # Should match the number of columns in rv_max
	y = moon_mass * SOLAR_MASS / EAR_MASS


	with PlotStyle('science', TEX_FONTS, LATEX_PLOT):
		plt.contourf(x,y,rv_max)
		cb = plt.colorbar()
		cb.set_label(label='Peak Radial Velocity (m/s)')
		
		CS=plt.contour(x,y,rv_max,levels=[500],linestyles='--',colors='white')
		positions = [(1,35)]
		label = plt.clabel(CS,fmt='500m/s', inline=False,rightside_up=False,colors='C3',fontsize='15',manual=positions)
		for l in label:
			l.set_rotation(20)
		
		plt.fill_between(t_roch,moon_mass[-1]*SOLAR_MASS/EAR_MASS,moon_mass*SOLAR_MASS/EAR_MASS,color='grey')
		plt.text(.15,70,'Roche limit',fontsize='13',rotation=90)

		plt.fill_between([hill_sphere,moon_period[-1]],moon_mass[0]*SOLAR_MASS/EAR_MASS,moon_mass[-1]*SOLAR_MASS/EAR_MASS,color='silver')
		plt.text(4e3,60,'Hill sphere limit',fontsize='13',rotation=90)

		plt.scatter(Kepler1708[0],Kepler1708[1],label='Kepler-1708 b-i',marker='o',color='snow')
		plt.scatter(Kepler1625[0],Kepler1625[1],label='Kepler-1625 b-i',marker='o',color='gold')
		plt.scatter(x_pos,y_pos,label='RV detection limit',marker='+',color='cyan')

		plt.legend(frameon=True,loc='upper center')
		plt.xscale('log')
		plt.xlim(moon_period[0],moon_period[-1])
		plt.ylim(moon_mass[0]*SOLAR_MASS/EAR_MASS,moon_mass[-1]*SOLAR_MASS/EAR_MASS)
		plt.xlabel('Orbital Period (days)')
		plt.ylabel(r'Moon Mass ($M_\oplus$)')
		plt.tick_params(which='both',color='red')
		plt.savefig(f'plots/contour_hill.png',dpi=300,bbox_inches='tight')

