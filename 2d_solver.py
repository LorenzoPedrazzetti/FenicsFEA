from __future__ import print_function
import matplotlib.pyplot as plt
from fenics import *
from dolfin import *
from dolfin.cpp.mesh import *
from mshr import *
import numpy as np
import pandas as pd 
from max_q_outflow_bulk import max_q
from light_abs_DT import T_absorption
#import sympy as sym

'''
all quantities in IS units
'''

path = input("insert path for n&k file ")
wvl = input("insert absorption wavelenght ")
dt = 100 #time step size
L = input("insert square dimension ")
grad = input("insert distance for temperature gradient ")
centre = (L/2, L/2) #domain centre
tol = 1E-14 #tolerance

temp0 = T_absorption(surf=L*L, path=path) #temperature increase due to light absorption at 10 um, ~301
Q = max_q(T=temp0, L=grad) #maxiumum heat outflow

#create the mesh, 2D
domain = Rectangle(Point(0.0, 0.0), Point(L, L))
mesh = generate_mesh(domain, 15)
V = FunctionSpace(mesh, 'P', 1)

#Dirichlet BC, is the T of the membrane, extracted from "pyro_Si3N4.py"
u_D = Constant(temp0)

'''BCs are only Neumann: on one side the flux is proportional to the source term (power absorbed gets dissipated),
on all the other facets the flux is 0 (suspended membrane separated from bulk)'''

#subdomain on which the flux is g
class Omega_0(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0, tol) or near(x[0], 1, tol)

#subdomain on which the flux is g
class Omega_1(SubDomain):
    def inside(self,x,on_boundary):
        return near(x[1], 0, tol) or near(x[1], 1, tol) 

boundary_markers = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
subdomain_0 = Omega_0()
subdomain_0.mark(boundary_markers, 0)

subdomain_1 = Omega_1()
subdomain_1.mark(boundary_markers, 1)

#redefine boundary integration measure
ds = Measure('ds', domain = mesh, subdomain_data = boundary_markers) 

'''and here the Neumann expression, the outgoing flux, calculated by menial application of the heat conduction law'''
g = Constant(Q)
 
#define initial conditions interpolate the initial T on the mesh
u_n = interpolate(u_D, V) 

#define variational problem
u = TrialFunction(V)
v = TestFunction(V)

'''source term'''
#f = Expression('1E-6*sin(pi/2*t)', degree=0,t=0) 
f = Constant(0)

#LN thermal diffusivity [m2/s]
alpha = input("insert material thermal diffusivity")
#LN thermal conductivity [W/m°K]
kappa = input("insert material thermal condictivity")

#list of boundary conditions, a dictionary including all the markers we set above
boundary_conditions = {0 : {'Neumann' : g}, 1 : {'Neumann' : g}}

bcs=[] #empty BC list, no dirichlet

#since the weak formulation features an integral over 'gammaN' implementation of a function thaat takes care of the integral splitting over the set is needed
integrals_N = []

for i in boundary_conditions:
    if 'Neumann' in boundary_conditions[i]:
        if boundary_conditions[i]['Neumann'] != 0:
            g = boundary_conditions[i]['Neumann']
            integrals_N.append(1/kappa*g*v*ds(i))


#define the C++ expression
a = u*v*dx + alpha*dt*dot(grad(u), grad(v))*dx
L = u_n*v*dx - alpha*dt*sum(integrals_N) #all non-expression terms must be taken out of the integral over Neumann boundary 

#time stepping
u = Function(V)
#t = 1E-6

#loo counter
loop_counter = 0 

'''
loop 0 calculation to init the first u value
'''
#call the solver
solve(a == L, u, bcs)

#temperatures vector
temperature = [u(centre)]
#time vector
time_passed = [0]

#the condition is phrased in therms of % of T with respect to initial temperature, in the centre of the membrane
#20% of the initial DT (reduced by 80%)
percent = 0.9

while u(centre) > 300+(temp0-300)*percent: 

    #update loop counter
    loop_counter+=1

    #check number of loop and escape if it's > 100
    if loop_counter == 5000:
        print("Maximum nuber of iterations reached .. try to increase dt.")
        break

    #update the current time and BC
    #t+=dt
    #f.t=t

    #call the solver
    solve(a == L, u, bcs)

    #compute the new Q outflow
    Q = max_q(T=u(centre))
    g.assign(Q)

    #update vectors
    temperature.append(u(centre))
    time_passed.append(time_passed[0]+loop_counter*dt)

    #update u_n to the current u value
    u_n.assign(u)

#export the last iteration as vtk file
vtkfile = File('/mnt/c/Users/LorenzoPedrazzetti/WSL/fenics_output_2d/no-arms-membrane/2d_solution.pvd')
vtkfile << u

#save data to csv
df_output_no_arms = pd.DataFrame(time_passed, columns=['time'])
df_output_no_arms['temperature'] = temperature
df_output_no_arms.to_csv('/mnt/c/Users/LorenzoPedrazzetti/WSL/fenics_output_2d/no-arms-membrane/output_no_arms.csv')

print("initial T is: ", temp0)
print("final T is: ", u(centre), " K")
print("time to cool off "+str(percent*100)+"% of the temperature increase: ", time_passed[loop_counter], " s")
print("numer of iterations: ", loop_counter)

'''error evaluation'''

#compare error at vertices for last iteration
u_e = interpolate(u_n, V)
error = errornorm(u_e, u, norm_type='l2')
print('t = %.2f: error = %.3g' % (loop_counter,error))

plt.plot(time_passed, temperature)
plt.show()
