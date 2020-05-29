

def max_q(T=350, T0=295, L=0.01):

    # L is the lenght over which the T changes, the extension of the linbo3 arms, used to calculate the gradient [cm]
    # n is the number of arms attached to the bulk
    # N is the integer that is used to divide the membrane side
    # linbo2 thermal conductivity [Wm^-1K^-1]
    k = 5.6
    # membrane T [K]
    # surroundings T0 [K]
    # arms lenght and width [m]
    l = 140E-6
    w = 600E-9

    # heat equation [W/m2]
    q = k * (T-T0)/(L)

    Q = q / l*l*w

    return Q

'''checked'''

print("heat outflow", max_q(), "W/m2")