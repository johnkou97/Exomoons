if __name__ == "__main__":
    '''
    This script is used to test starry's capabilities to model a zone by creating multiple spots.
    Starry seems to n ot be able to model more than 11 spots, making it impossible to model a zone
    as it requires more spots than that.
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    import starry
    import astropy.units as u
    from constants import *

    starry.config.lazy = False
    starry.config.quiet = True

    beta_mass = 11.9    # in u.Mjup
    beta_rad = 1.65     # in u.Rjup
    theta = np.linspace(-180, 180, 1000)
    time_ros = np.linspace(0,BETA_DAY,len(theta))/(60*60)

    # Testing simple addition of spots

    A = starry.Primary(
        starry.Map(ydeg=10, udeg=2, rv=True, amp=1, veq=BETA_VEQ, alpha=0),
        r=beta_rad,
        m=beta_mass,
        length_unit=u.Rjup,
        mass_unit=u.Mjup)
    map_beta=A.map
    map_beta[1] = 0.5
    map_beta[2] = 0.25
    map_beta.spot(contrast=.95 ,radius=16 ,lat=30, lon=0)
    rv_1_ros = map_beta.rv(theta=theta)

    A = starry.Primary(
        starry.Map(ydeg=10, udeg=2, rv=True, amp=1, veq=BETA_VEQ, alpha=0),
        r=beta_rad,
        m=beta_mass,
        length_unit=u.Rjup,
        mass_unit=u.Mjup)
    map_beta=A.map
    map_beta[1] = 0.5
    map_beta[2] = 0.25
    map_beta.spot(contrast=.95 ,radius=16 ,lat=-30, lon=-45)
    rv_2_ros = map_beta.rv(theta=theta)

    A = starry.Primary(
        starry.Map(ydeg=10, udeg=2, rv=True, amp=1, veq=BETA_VEQ, alpha=0),
        r=beta_rad,
        m=beta_mass,
        length_unit=u.Rjup,
        mass_unit=u.Mjup)
    map_beta=A.map
    map_beta[1] = 0.5
    map_beta[2] = 0.25
    map_beta.spot(contrast=.95 ,radius=16 ,lat=30, lon=0)
    map_beta.spot(contrast=.95 ,radius=16 ,lat=-30, lon=-45)
    rv_3_ros = map_beta.rv(theta=theta)

    # Plot the radial velocity
    plt.figure(figsize=(12, 5))
    #plt.title('Comparing the two methods')
    plt.plot(time_ros, rv_1_ros,label='spot 1')
    plt.plot(time_ros, rv_2_ros,label='spot 2')
    plt.plot(time_ros, rv_1_ros+rv_2_ros,label='spot 1 + spot 2')
    plt.plot(time_ros, rv_3_ros,label='spot 3')
    plt.legend()
    plt.xlabel("Time [hours]", fontsize=24)
    plt.ylabel("Radial velocity [m/s]", fontsize=24)
    plt.grid()
    plt.savefig('plots_rossiter/addition_1_spot.png')

    # Testing how many spots can starry handle

    A = starry.Primary(
        starry.Map(ydeg=10, udeg=2, rv=True, amp=1, veq=BETA_VEQ, alpha=0),
        r=beta_rad,
        m=beta_mass,
        length_unit=u.Rjup,
        mass_unit=u.Mjup)
    map_beta=A.map
    map_beta[1] = 0.5
    map_beta[2] = 0.25
    number_of_spots=4
    for i in range(number_of_spots):
        map_beta.spot(contrast=.95 , radius=15 ,lat=30, lon=i*(360/number_of_spots))
    rv_4_ros = map_beta.rv(theta=theta)


    A = starry.Primary(
        starry.Map(ydeg=10, udeg=2, rv=True, amp=1, veq=BETA_VEQ, alpha=0),
        r=beta_rad,
        m=beta_mass,
        length_unit=u.Rjup,
        mass_unit=u.Mjup)
    map_beta=A.map
    map_beta[1] = 0.5
    map_beta[2] = 0.25
    number_of_spots=6
    for i in range(number_of_spots):
        map_beta.spot(contrast=.95 ,radius=15 ,lat=30, lon=i*(360/number_of_spots))
    rv_5_ros = map_beta.rv(theta=theta)

    A = starry.Primary(
        starry.Map(ydeg=10, udeg=2, rv=True, amp=1, veq=BETA_VEQ, alpha=0),
        r=beta_rad,
        m=beta_mass,
        length_unit=u.Rjup,
        mass_unit=u.Mjup)
    map_beta=A.map
    map_beta[1] = 0.5
    map_beta[2] = 0.25
    number_of_spots=8
    for i in range(number_of_spots):
        map_beta.spot(contrast=.95 ,radius=15 ,lat=30, lon=i*(360/number_of_spots))
    rv_6_ros = map_beta.rv(theta=theta)

    A = starry.Primary(
        starry.Map(ydeg=10, udeg=2, rv=True, amp=1, veq=BETA_VEQ, alpha=0),
        r=beta_rad,
        m=beta_mass,
        length_unit=u.Rjup,
        mass_unit=u.Mjup)
    map_beta=A.map
    map_beta[1] = 0.5
    map_beta[2] = 0.25
    number_of_spots=10
    for i in range(number_of_spots):
        map_beta.spot(contrast=.95 ,radius=15 ,lat=30, lon=i*(360/number_of_spots))
    rv_7_ros = map_beta.rv(theta=theta)

    A = starry.Primary(
        starry.Map(ydeg=10, udeg=2, rv=True, amp=1, veq=BETA_VEQ, alpha=0),
        r=beta_rad,
        m=beta_mass,
        length_unit=u.Rjup,
        mass_unit=u.Mjup)
    map_beta=A.map
    map_beta[1] = 0.5
    map_beta[2] = 0.25
    number_of_spots=12
    for i in range(number_of_spots):
        map_beta.spot(contrast=.95 ,radius=15 ,lat=30, lon=i*(360/number_of_spots))
    rv_8_ros = map_beta.rv(theta=theta)

    A = starry.Primary(
        starry.Map(ydeg=10, udeg=2, rv=True, amp=1, veq=BETA_VEQ, alpha=0),
        r=beta_rad,
        m=beta_mass,
        length_unit=u.Rjup,
        mass_unit=u.Mjup)
    map_beta=A.map
    map_beta[1] = 0.5
    map_beta[2] = 0.25
    number_of_spots=14
    for i in range(number_of_spots):
        map_beta.spot(contrast=.95 ,radius=15 ,lat=30, lon=i*(360/number_of_spots))
    rv_9_ros = map_beta.rv(theta=theta)

    # Plot the radial velocity
    plt.figure(figsize=(12, 5))
    plt.plot(time_ros, rv_4_ros,label='4 spots' ,color='b')
    plt.plot(time_ros, rv_5_ros,label='6 spots' ,color='orange')
    plt.plot(time_ros, rv_6_ros,label='8 spots' ,color='purple')
    plt.plot(time_ros, rv_7_ros,label='10 spots' ,color='r')
    plt.plot(time_ros, rv_8_ros,label='12 spots',color='g')
    plt.plot(time_ros, rv_9_ros,label='16 spots',color='y')
    plt.legend()
    plt.xlabel("Time [hours]", fontsize=24)
    plt.ylabel("Radial velocity [m/s]", fontsize=24)
    plt.grid()
    plt.savefig('plots_rossiter/number_of_spots.png')

    # Plot the radial velocity
    plt.figure(figsize=(12, 5))
    plt.plot(time_ros, rv_6_ros,label='8 spots' ,color='purple')
    plt.plot(time_ros, rv_7_ros,label='10 spots' ,color='r')
    plt.plot(time_ros, rv_8_ros,label='12 spots',color='g')
    plt.plot(time_ros, rv_9_ros,label='16 spots',color='y')
    plt.xlabel("Time [hours]", fontsize=24)
    plt.ylabel("Radial velocity [m/s]", fontsize=24)
    plt.legend()
    plt.grid()
    plt.savefig('plots_rossiter/number_of_spots_2.png')

    # Plot the radial velocity
    plt.figure(figsize=(12, 5))
    plt.plot(time_ros, rv_8_ros,label='12 spots',color='g')
    plt.plot(time_ros, rv_9_ros,label='16 spots',color='y')
    plt.legend()
    plt.xlabel("Time [hours]", fontsize=24)
    plt.ylabel("Radial velocity [m/s]", fontsize=24)
    plt.grid()
    plt.savefig('plots_rossiter/number_of_spots_3.png')

    # Testing adding 4 spots at 30 and 4 spots at -30

    A = starry.Primary(
        starry.Map(ydeg=10, udeg=2, rv=True, amp=1, veq=BETA_VEQ, alpha=0),
        r=beta_rad,
        m=beta_mass,
        length_unit=u.Rjup,
        mass_unit=u.Mjup)
    map_beta=A.map
    map_beta[1] = 0.5
    map_beta[2] = 0.25
    number_of_spots=4
    for i in range(number_of_spots):
        map_beta.spot(contrast=.95 ,radius=15 ,lat=30, lon=i*(360/number_of_spots))
    rv_10_ros = map_beta.rv(theta=theta)

    A = starry.Primary(
        starry.Map(ydeg=10, udeg=2, rv=True, amp=1, veq=BETA_VEQ, alpha=0),
        r=beta_rad,
        m=beta_mass,
        length_unit=u.Rjup,
        mass_unit=u.Mjup)
    map_beta=A.map
    map_beta[1] = 0.5
    map_beta[2] = 0.25
    number_of_spots=4
    for i in range(number_of_spots):
        map_beta.spot(contrast=.95 ,radius=15 ,lat=-30, lon=(i+.5)*(360/number_of_spots))
    rv_11_ros = map_beta.rv(theta=theta)

    A = starry.Primary(
        starry.Map(ydeg=10, udeg=2, rv=True, amp=1, veq=BETA_VEQ, alpha=0),
        r=beta_rad,
        m=beta_mass,
        length_unit=u.Rjup,
        mass_unit=u.Mjup)
    map_beta=A.map
    map_beta[1] = 0.5
    map_beta[2] = 0.25
    number_of_spots=4
    for i in range(number_of_spots):
        map_beta.spot(contrast=.95 ,radius=15 ,lat=30, lon=i*(360/number_of_spots))
    for i in range(number_of_spots):
        map_beta.spot(contrast=.95 ,radius=15 ,lat=-30, lon=(i+.5)*(360/number_of_spots))
    rv_12_ros = map_beta.rv(theta=theta)

    # Plot the radial velocity
    plt.figure(figsize=(12, 5))
    plt.plot(time_ros, rv_10_ros,label='4 spots at 30',color='b')
    plt.plot(time_ros, rv_11_ros,label='4 spots at -30',color='r')
    plt.plot(time_ros, rv_12_ros,label='4 spots at 30 and 4 spots at -30',color='y')
    plt.plot(time_ros, rv_10_ros+rv_11_ros,label='4 spots at 30 plus 4 spots at -30',color='orange')
    plt.plot(time_ros, rv_6_ros,label='8 spots at 30' ,color='purple')
    plt.legend()
    plt.xlabel("Time [hours]", fontsize=24)
    plt.ylabel("Radial velocity [m/s]", fontsize=24)
    plt.grid()
    plt.savefig('plots_rossiter/addition_of_spots.png')

    # Plot the radial velocity
    plt.figure(figsize=(12, 5))
    plt.plot(time_ros, rv_12_ros,label='4 spots at -30 and 4 spots at -30',color='y')
    plt.plot(time_ros, rv_10_ros+rv_11_ros,label='4 spots at 30 plus 4 spots at -30',color='orange')
    plt.plot(time_ros, rv_6_ros,label='8 spots' ,color='purple',linestyle='--')
    plt.legend()
    plt.xlabel("Time [hours]", fontsize=24)
    plt.ylabel("Radial velocity [m/s]", fontsize=24)
    plt.grid()
    plt.savefig('plots_rossiter/addition_of_spots_2.png')