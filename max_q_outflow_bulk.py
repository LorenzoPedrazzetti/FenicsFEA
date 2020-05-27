'''
this function describes how much heat flow the membrane supports depending on the extension of the "arms" 
that connect it to the bulk.

it is a simple implementetion of the heat equation:

q = -k * grad(T)[x,y]

it's a 1st degree differential equation, so I need to set IC and BC.

IC: this is the temperature of the membrane, after absorbing the light; it's obtained from the "sin_light_abs_DT.py"

BC: the BCs are mixed, since we have only outflow through a portion of the boundary (Neumann), while the rest of the system is kept adiabatic

since the linbo square is 2D, the solution of the equation is:

q = -2*k*dT/dx - 2*k*dT/dy ~ -4*k*dT/dx

when q is known as [J/W^2], we can multiply it by the surface of the arms.

'''

def max_q(T=301.5, T0=300, L):

    # L is the lenght over which the T changes, the extension of the linbo3 arms, used to calculate the gradient [cm]
    # n is the number of arms attached to the bulk
    # N is the integer that is used to divide the membrane side
    k = 0.056 # linbo2 thermal conductivity [Wcm^-1K^-1]
    # membrane T [K]
    # surroundings T0 [K]
    # arms lenght and width [cm]
    l = 140E-4
    w = 600E-7

    # heat equation [Wcm^-2]
    q = k * (T-T0)/(L) 
 
    # total heat flow [W]
    surf = 4*w*l
    Q = q * surf

    return Q

'''checked'''