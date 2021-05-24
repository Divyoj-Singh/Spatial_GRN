# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 18:17:38 2021

@author: csb
"""
#%% Importing libraries:
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

#%% Parameters:
alpha=1;
gamma=0.05;
t_max=150;
kth=2;

#%% functions:
def response_time(sol):
    x=np.squeeze(sol.y.T);time=np.array(sol.t)
    
    steady_state=x[-1]
    
    for i in range(np.size(x)):
        if x[i]> steady_state/2:
            return [time[i],steady_state]
            break

def simple_regulation(t, y): return alpha - gamma*y
def positive_regulation(t, y,nAR): return alpha*((y**nAR)/(kth**nAR+y**nAR)) - gamma*y
def negative_regulation(t, y,nAR): return alpha*((kth**nAR)/(kth**nAR+y**nAR)) - gamma*y

#%% Part A :
print ("part a")
t = np.linspace(0, t_max, 1000)

sol1 = solve_ivp(simple_regulation,   [0, t_max], [1],t_eval=t)
sol2 = solve_ivp(positive_regulation, [0, t_max], [1],args=(2,),t_eval=t)
sol3 = solve_ivp(negative_regulation, [0, t_max], [1],args=(2,),t_eval=t)

## Plotting:
f, ax = plt.subplots(figsize=(4,3))        
ax.plot(t, sol1.y.T/(alpha/gamma))
ax.plot(t, sol2.y.T/(alpha/gamma))
ax.plot(t, sol3.y.T/(alpha/gamma))
ax.set_ylabel("normalised levels",fontsize=10)
ax.set_xlabel("Time",fontsize=10)
ax.legend(['simple', 'positive','negative'], shadow=True)
f.suptitle('Different types of regulations')
plt.subplots_adjust(top=0.85, bottom=0.20, left=0.15, right=0.95, hspace=0.25,wspace=0.5)
plt.savefig("1a.png",dpi=200)
plt.close()
       
print("response time and non-normalised steady state of simple regulation:",response_time(sol1))        
print("response time and non-normalised steady state of positive regulation:",response_time(sol2))        
print("response time and non-normalised steady state of negative regulation:",response_time(sol3)) 


#%% Part B :

print ("part b")
t = np.linspace(0, t_max, 1000)
## positive regulation:
sol1 = solve_ivp(positive_regulation, [0, t_max], [3],args=(2,),t_eval=t)
sol2 = solve_ivp(positive_regulation, [0, t_max], [3],args=(5,),t_eval=t)


## Plotting:
f, ax = plt.subplots(figsize=(4,3))        
ax.plot(t, sol1.y.T/(alpha/gamma))
ax.plot(t, sol2.y.T/(alpha/gamma))
ax.set_ylabel("normalised levels",fontsize=10)
ax.set_xlabel("Time",fontsize=10)
ax.legend(['n_AR=2', 'n_AR=5'], shadow=True)
f.suptitle('positive regulation')
plt.subplots_adjust(top=0.85, bottom=0.20, left=0.15, right=0.95, hspace=0.25,wspace=0.5)
plt.savefig("1b_positive_regulation.png",dpi=200)
plt.close()
       
print("response time and non-normalised steady state of positive regulation:",response_time(sol1))        
print("response time and non-normalised steady state of strong positive regulation:",response_time(sol2))               

## negative regulation:
sol1 = solve_ivp(negative_regulation, [0, t_max], [1],args=(2,),t_eval=t)
sol2 = solve_ivp(negative_regulation, [0, t_max], [1],args=(5,),t_eval=t)

## Plotting:
f, ax = plt.subplots(figsize=(4,3))        
ax.plot(t, sol1.y.T/(alpha/gamma))
ax.plot(t, sol2.y.T/(alpha/gamma))
ax.set_ylabel("normalised levels",fontsize=10)
ax.set_xlabel("Time",fontsize=10)
ax.legend(['n_AR=2', 'n_AR=5'], shadow=True)
f.suptitle('negative_regulation')
plt.subplots_adjust(top=0.85, bottom=0.20, left=0.15, right=0.95, hspace=0.25,wspace=0.5)
plt.savefig("1b_negative_regulation.png",dpi=200)
plt.close()
       
print("response time and non-normalised steady state of negative_regulation regulation:",response_time(sol1))        
print("response time and non-normalised steady state of strong negative_regulation regulation:",response_time(sol2))               

#%% Part C :

print ("part c")
def mixed_regulation(t, y,nact,ninh,kth_act,kth_inh): return alpha*((kth_inh**ninh)/(kth_inh**ninh+y**ninh))+(alpha*(y**nact)/(kth_act**nact+y**nact)) - gamma*y

t = np.linspace(0, t_max, 1000)
## positive regulation:
sol1 = solve_ivp(mixed_regulation, [0, t_max], [1,2,4,5,10,20,25],args=(10,2,5,15),t_eval=t)

## Plotting:
f, ax = plt.subplots(figsize=(4,3))        
ax.plot(t, sol1.y.T/(alpha/gamma))
ax.set_ylabel("normalised levels",fontsize=10)
ax.set_xlabel("Time",fontsize=10)
ax.legend(['1','2','4','5','10','20','25'], shadow=True)
f.suptitle('both regulations together')
plt.subplots_adjust(top=0.85, bottom=0.20, left=0.15, right=0.95, hspace=0.25,wspace=0.5)
plt.savefig("1c.png",dpi=200)
plt.close()