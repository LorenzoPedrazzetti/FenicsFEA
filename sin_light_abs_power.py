from numpy import exp, pi
from pandas import read_csv

'''
I law of thermodynamics:
understanding how much of the incident power is converted to heat from photons absorption in the Si3N4 layer is very hard
assuming that the only energy exchange the layer has with the universe is light absorption, the I TDs law can help

'''

def T_absorption(P0=5, surf=0.02*0.02, thickness=0.01 , wvl=5, z=0.002):
    
    ''' it is still with surface power incident, shall we change it???'''

    
    # P0: power density [W/m^3]
    # surf: effective absorption area [m^2]
    # wvl: wavelenght of the incident light [um]
    # z: coordinate along the material thickness [m]
    # k: absorption coefficient for the material, it is tabulated against the wavelenght
    
    # total incident power modulus
    Pp0=P0*surf
    vol = surf*thickness
    '''
    absorption is extracted from a csv table in the local folder
    '''

    # try handle to pandas upload of the n&k data
    try:
        opt_param = read_csv('/home/lpedrazzetti/heat_equation_fenics/linbo_sin_3d/2d_problem/neumann_bc/SiN_nk.csv')
        print("n&k data successfully read .. ")
    except:
        print("check n&k data .. ")

    # check for the absorption value related to the prompted wavelenght
    for i,j in opt_param['wl'].items():
        try:
            if j == wvl: k = i
        except: print("wavelenght not found .. ")

    # absorption formula (alpha)
    alpha = 4 * pi * k / (wvl)

    # lambert's law for light absorption in the material, linbo3 is considered 100% transparent
    P = Pp0 * exp(-alpha*z) / vol # [W/m3]
    
    return P

print(T_absorption() , "W/m3")