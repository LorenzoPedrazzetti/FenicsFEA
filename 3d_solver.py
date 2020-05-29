#fenics and dolfin imports
from __future__ import print_function
from fenics import *
from dolfin import *
from dolfin.cpp.mesh import *
from mshr import *
#math libraries and plotting imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#custom methods imports
from max_q_outflow import max_q #detached membrane experiment
from sin_light_abs_power import T_absorption #row radiation power as source term, no energy balance

'''
Variable describing the specific prolem at hand.
All quantities in ## IS ## units!!!
'''

tol = 1E-14 #tolerance used to compute mesh points 
t = 0 #time begins [s]
freq = 10 #chopping frequency [Hz] 
dt = 1/freq #time step size [m]
iterations = 600
lenght = 0.02 #membrane sides [m]
thickness = 0.01 #membrane thickness [m]
centre = (lenght-tol, lenght-tol, thickness) #domain centre 
N = 2 #integer by which the membrane side is diveded
power = 80 #peak power DENSITY impinging on the membrane [W/m2]
temp0 = 295 #temperature of the membrane @ t=0 [K]

#create the mesh, 3D
mesh = BoxMesh(Point(0, 0, 0), Point(lenght, lenght, thickness), 10, 10, 10)
V = FunctionSpace(mesh, 'P', 1)

'''BCs are only Neumann: on one side the flux is proportional to the source term (power absorbed gets dissipated),
on all the other facets the flux is 0 (suspended membrane separated from bulk)'''

#subdomain on which the flux is 0
class Omega_0(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] > lenght/N+tol or x[1] > lenght/N+tol

#subdomain on which the flux is g
class Omega_1(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] <= lenght/N+tol or x[1] <=lenght/N+tol

boundary_markers = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
subdomain_0 = Omega_0()
subdomain_0.mark(boundary_markers, 0)

subdomain_1 = Omega_1()
subdomain_1.mark(boundary_markers, 1)

#redefine boundary integration measure
ds = Measure('ds', domain = mesh, subdomain_data = boundary_markers) 

#interpolate the initial temperature over the mesh to define the ICs
u_D = Constant(temp0)
u_n = interpolate(u_D, V) 

#material properties
#LN thermal diffusivity [m^2/s]
alpha = Constant(0.97E-6)
#LN thermal conductivity [W/mK] 
kappa = Constant(1.3)
#maxiumum heat outflow
Q = max_q(T=temp0) 
g = Expression('Q', degree=0, Q=Q)
#source term
source = T_absorption(P0=power)
gamma = freq
f = Expression('source * abs(sin(gamma * pi/2 * t))', degree=0, source=source, gamma=gamma, t=t)

#define variational problem functions
u = TrialFunction(V)
v = TestFunction(V)

#list of boundary conditions, a dictionary including all the markers we set above
boundary_conditions = {0 : {'Neumann' : 0}, 
                       1 : {'Neumann' : g}}

bcs=[]

#since the weak formulation features an integral over 'gammaN' implementation of a function thaat takes care of the integral splitting over the set is needed
integrals_N = []

for i in boundary_conditions:
    if 'Neumann' in boundary_conditions[i]:
        if boundary_conditions[i]['Neumann'] != 0:
            g = boundary_conditions[i]['Neumann']
            integrals_N.append(alpha/kappa*dt*g*v*ds(i))

F = u*v*dx + alpha*dt*dot(grad(u), grad(v))*dx - u_n*v*dx - alpha/kappa*f*dt*v*dx + sum(integrals_N)
a, L = lhs(F), rhs(F) 

'''
loop 0 calculation to init the first u value
'''

#call the solver
u = Function(V)
solve(a == L, u, bcs)

#loop counter
loop_counter = 0 
#init the arrays fot temp and time
temperature = [u(centre)]
time_passed = [0]

while loop_counter < iterations :

    #update loop counter
    loop_counter += 1

    #update time
    t = loop_counter*dt

    #update source term with new time step
    f.t=t

    #check number of loop and escape if it's > 5000
    if loop_counter == 5000:
        print("Maximum nuber of iterations reached .. try to increase dt.")
        break

    #call the solver
    solve(a == L, u, bcs)

    #compute the new Q outflow
    Q = max_q(T=u(centre))
    g.Q = Q

    #append temperature list
    temperature.append(u(centre))
    
    #append time list
    time_passed.append(loop_counter*dt)

    #update u_n to the current u value
    u_n.assign(u)

#export the last iteration as vtk file
vtkfile = File('/mnt/c/Users/LorenzoPedrazzetti/WSL/fenics_output_2d/arms/2d_solution.pvd')
vtkfile << u

#save data to csv
df_output_no_arms = pd.DataFrame(time_passed, columns=['time'])
df_output_no_arms['temperature'] = temperature
df_output_no_arms.to_csv('/mnt/c/Users/LorenzoPedrazzetti/WSL/fenics_output_2d/arms/output_arms.csv')

#print some infos
print("with N: ", N)
print("initial T is: ", temp0)
print("final T is: ", u(centre), "K")
print("chopping freq is: ", 1/dt, "Hz")
print("numer of iterations: ", loop_counter)

#plot the heating curve
plt.plot(time_passed, temperature)
plt.xlabel("time (s)")
plt.ylabel("Temperature (�K)")
plt.show()
