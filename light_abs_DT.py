from numpy import exp, pi
from pandas import read_csv

'''
I law of thermodynamics:
understanding how much of the incident power is converted to heat from photons absorption in the Si3N4 layer is very hard
assuming that the only energy exchange the layer has with the universe is light absorption, the I TDs law can help

Q = work + dU/dt

work = 0
Q=power absorbed

hence: dU/dt = P

where P is the power absorbed in the Si3N4

Q can be expressend like: Q = DT*m*cp

the absorption gets also integrated over the entire layer thickness.

There is no wavelenght integration at this stage.
'''

def T_absorption(P0=1, surf=140E-6*140E-6, wvl=10, z=300E-6, T0=299, t=10, path=path):
    
    # P0: power density [W/m^2]
    # surf: effective absorption area [m^2]
    # wvl: wavelenght of the incident light [um]
    # z: coordinate along the material thickness [m]
    # k: absorption coefficient for the Si3N4, it is tabulated against the wavelenght
    # T0: initial body temperature [K]
    # t: time [s]
    
    # mass *specific heat capacity of Si3N4 [J]
    mass_cp=1.65E-8
    # total incident power modulus
    Pp0=P0*surf

    '''
    absorption is extracted from a csv table in the local folder
    '''

    # try handle to pandas upload of the n&k data
    try:
        opt_param = read_csv(path)
        print("n&k data successfully read .. ")
    except:
        print("check n&k data .. ")

    # check for the absorption value related to the prompted wavelenght
    for i,j in opt_param['wl'].items():
        try:
            if j == wvl: k = i
        except: print("wavelenght not found .. ")

    # absorption formula (alpha)
    alpha = 4 * pi * k / wvl
    
    # lambert's law for light absorption in the material, linbo3 is considered 100% transparent
    P = Pp0 * exp(-alpha*z)

    # energy balance
    y = T0 + P*t/mass_cp

    return y

'''checked'''