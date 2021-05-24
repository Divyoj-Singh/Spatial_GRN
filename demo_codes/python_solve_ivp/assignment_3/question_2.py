# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 12:45:49 2021

@author: csb
"""

#%% Importing libraries:
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()
#%% Parameters:
k=4;
mu=1;
t_max=100;
n=2;

#%% functions:
def toggle_switch(t, z,k1,mu1,k2,mu2,n):
    x, y = z
    dxdt= k1/(1+y**n) - mu1*x;
    dydt= k2/(1+x**n) - mu1*y;
    return [dxdt,dydt] 
def state(sol):
    x=sol.y[0][-1];
    y=sol.y[1][-1];
    if x/y <1:
        return 0
    if x/y >1:
        return 1
    
#%% Part A :
print ("part a")
t = np.linspace(0, t_max, 100)

f, ax = plt.subplots(figsize=(4,3))        
for i in range(1000):
    
    sol = solve_ivp(toggle_switch, [0, t_max], 5*np.random.random(2),args=(k,mu,k,mu,n),t_eval=t)
    x=sol.y[0].T;
    y=sol.y[1].T;

    ## Plotting:
    ax.plot(t[5:], np.log2(x[5:]/y[5:]))

ax.set_ylabel("log2(x/y)",fontsize=10)
ax.set_xlabel("Time",fontsize=10)
f.suptitle("Bistable toggle switch")
plt.subplots_adjust(top=0.85, bottom=0.20, left=0.15, right=0.95, hspace=0.25,wspace=0.5)
plt.savefig("2a"+".png",dpi=200)    
plt.close()

#%% Part B :
print ("part b")
t = np.linspace(0, t_max, 100)
fraction_array=[];
f, ax = plt.subplots(figsize=(4,3))        
        
        
for experiment in range(10):
    fraction=0
    for cell in range(1000):
        sol = solve_ivp(toggle_switch, [0, t_max], 5*np.random.random(2),args=(k,mu,k,mu,n),t_eval=t)
        if state(sol)==0:
            fraction+=1/1000;
    fraction_array.append(fraction)        
ax.hist(fraction_array,color = "skyblue", ec="skyblue")
ax.set_xlim([0,1])
ax.set_ylabel("frequency distribution",fontsize=10)
ax.set_xlabel("fraction of cases where X/Y goes to 0",fontsize=10)

f.suptitle("frequency distribution for balanced circuit")
plt.subplots_adjust(top=0.85, bottom=0.20, left=0.15, right=0.95, hspace=0.25,wspace=0.5)
plt.savefig("2b"+".png",dpi=200)    
plt.close()



#%% Part C:
print ("part c")
t = np.linspace(0, t_max, 100)
fraction_array=[];
f, ax = plt.subplots(figsize=(4,3))        
        
        
for experiment in range(10):
    fraction=0
    for cell in range(1000):
        sol = solve_ivp(toggle_switch, [0, t_max], 5*np.random.random(2),args=(k+0.2*k,mu,k,mu,n),t_eval=t)
        if state(sol)==0:
            fraction+=1/1000;
    fraction_array.append(fraction)        
ax.hist(fraction_array,color = "skyblue", ec="skyblue")
ax.set_xlim([0,1])
ax.set_ylabel("frequency distribution",fontsize=10)
ax.set_xlabel("fraction of cases where X/Y goes to 0",fontsize=10)

f.suptitle("frequency distribution for unbalanced circuit")
plt.subplots_adjust(top=0.85, bottom=0.20, left=0.15, right=0.95, hspace=0.25,wspace=0.5)
plt.savefig("2c"+".png",dpi=200)    
plt.close()

#%% Part D :
print ("part d")
t = np.linspace(0, t_max, 100)

f, ax = plt.subplots(figsize=(4,3))        
for i in range(1000):
    ## changed n to 1.5
    sol = solve_ivp(toggle_switch, [0, t_max], 5*np.random.random(2),args=(k,mu,k,mu,1),t_eval=t)
    x=sol.y[0].T;
    y=sol.y[1].T;

    ## Plotting:
    ax.plot(t[5:], x[5:]/y[5:])

ax.set_ylabel("x/y",fontsize=10)
ax.set_xlabel("Time",fontsize=10)
f.suptitle("monostable toggle switch")
plt.subplots_adjust(top=0.85, bottom=0.20, left=0.15, right=0.95, hspace=0.25,wspace=0.5)
plt.savefig("2d"+".png",dpi=200)    
plt.close()
