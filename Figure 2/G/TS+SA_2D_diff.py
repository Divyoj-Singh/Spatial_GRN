# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 21:08:56 2021
@author: N.V
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 18:47:06 2021
@author: N.V
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import os
# note that this must be executed before 'import numba'
os.environ['NUMBA_DISABLE_INTEL_SVML'] = '1'
from numba import njit
import time as process_time

 #Degradtion rate
gamma_A=0.1 ; gamma_B=0.1;
#Transcription rate
g_A=5 ; g_B=5;


# @njit
def hill(X,X0,l,n):
    #H+
    H1=(X**n)/(X0**n+X**n)
    #H-
    H2=(X0**n)/(X**n+X0**n)
    #Hill function =H-(+) lambda*H+
    H=H2+l*H1
    return H







def interactions_sa(X):
    dXdt=np.zeros(np.shape(X))
    X=np.asarray(X)
   
    #Degradtion rate
    gamma_A=0.1; gamma_B=0.1;
    #Transcription rate
    g_A=5 ; g_B=5;
    #Hills function threshold 
    A0B=120 ; B0A=120;
    A0A=80;B0B=80;
    
    #Co-operativity 
    nAtoB=1;
    nBtoA=1;
    nAtoA=4;
    nBtoB=4;
    
    
    #fold change
    lambda_AtoB =0.1
    lambda_BtoA=0.1
    lambda_AtoA=5
    lambda_BtoB=5
    
    #Equation=
    dXdt[0,:,:]=g_A*hill(X[0,:,:],A0A,lambda_AtoA,nAtoA)*hill(X[1,:,:],B0A,lambda_BtoA,nBtoA) -gamma_A*X[0,:,:]
    dXdt[1,:,:]=g_B*hill(X[1,:,:],B0B,lambda_BtoB,nBtoB)*hill(X[0,:,:],A0B,lambda_AtoB,nAtoB) -gamma_B*X[1,:,:]
    
    return dXdt





# dimsension size, mm
w = h = 10.
# intervals in x-, y- directions, mm
dx = dy = 0.1
#Diffusivity
D = 0.0001
alpha=10; gamma=0.5

u_initial, v_initial = 50, 50
T_max=250;
nx, ny = int(w/dx), int(h/dy)

dx2, dy2 = dx*dx, dy*dy
dt = 0.01; #dx2 * dy2 / (2 * D * (dx2 + dy2))
resolution=100/dt;
#Setting the initial conditions

u0 = u_initial * np.ones((nx, ny))+ 10*np.random.random((nx, ny))
v0 = v_initial * np.ones((nx, ny))+ 10*np.random.random((nx, ny))


u = u0.copy();v = v0.copy()






# @njit
def do_timestep(u0, u, v0, v):
    # Propagate with forward-difference in time, central-difference in space
     X0= [ u0[1:-1,1:-1] , v0[1:-1,1:-1] ]
    
     F=interactions_sa(X0)
    
     u[1:-1, 1:-1] = u0[1:-1, 1:-1] + D * dt * ((u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] + u0[:-2, 1:-1])/dx2+ (u0[1:-1, 2:] - 2*u0[1:-1, 1:-1] + u0[1:-1, :-2])/dy2 )+ dt*F[0,:,:] 
    
     v[1:-1, 1:-1] = v0[1:-1, 1:-1] + D * dt * ((v0[2:, 1:-1] - 2*v0[1:-1, 1:-1] + v0[:-2, 1:-1])/dx2+ (v0[1:-1, 2:] - 2*v0[1:-1, 1:-1] + v0[1:-1, :-2])/dy2 ) +dt*F[1,:,:]
    
    ## Boundary Condition along x axis:
     u[0,:],u[-1,:]=u[1,:],u[-2,:]
     v[0,:],v[-1,:]=v[1,:],v[-2,:]
    
    ## Boundary Condition along y axis:
     u[:,0],u[:,-1]=u[:,1],u[:,-2]
     v[:,0],v[:,-1]=v[:,1],v[:,-2]
    
    
     u0 = u.copy();v0 = v.copy()
     return u0, u, v0, v


def plotconcmap(u_k, k):
    # Clear the current plot figure
    plt.clf()
    plt.rcParams.update({'font.size':14})
    #plt.title(f"A at t = {k:.3f} min",weight='bold')
    plt.xlabel("x",weight='bold')
    plt.ylabel("y",weight='bold')

    # This is to plot u_k (u at time-step k)
    plt.pcolormesh(u_k, cmap=plt.cm.jet, vmin=0, vmax=4 ,shading='gouraud')
    plt.colorbar()

    return plt

def animate(k):
    plotconcmap(uk[k]/(g_A/gamma_A), k)


# Number of timesteps
nsteps = int(T_max/dt)

t_start=process_time.time()
uk = np.empty((nx, ny));vk=u.copy()
for m in range(nsteps):
    u0, u, v0, v = do_timestep(u0, u, v0, v)
    
    if (m%resolution)==0:
        print("At time = ",m*dt)
        vk=np.dstack((vk,v));
        uk=np.dstack((uk,u));
    
uk=np.moveaxis(uk, -1, 0)
vk=np.moveaxis(vk, -1, 0)
print("Simulation time (in min)= ",(process_time.time()-t_start)/60)
#%%

anim = animation.FuncAnimation(plt.figure(), animate,interval=100,frames=np.size(uk,0),repeat=False)
anim.save("TS_SA_bistable_03_09.mp4")

print("Total time (in min)= ",(process_time.time()-t_start)/60)
