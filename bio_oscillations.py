#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 15:06:55 2022

@author: Cianna Garibaldi, Arna Ramirez, Atom Ramirez, Xinyu Liu, Ivan Ruiz
"""

"""
Note:removed project comments
"""

#%%

import numpy as np
import scipy
from scipy import integrate
import matplotlib.pyplot as plt 

def rate (V, t):
    H, P, G = V
    n = 9; k1 = 0.2; k2 = 0.2; k3 = 0.2 
    
    dHdt = 1/(1+G**n)-k1*H
    dPdt = H - k2*P
    dGdt = P - k3*G
    
    return np.array([dHdt, dPdt, dGdt])

import math 
res = math.sqrt(2)

tmin = 0
tmax = 100
Nt = 2000

t = np.linspace(tmin, tmax, Nt+1)

V0 = np.array([0.1, 0.9, 2.0])
Ans = integrate.odeint(rate, V0, t)
V = Ans.T
H, P, G = V
plt.figure(1)
plt.plot(t,H,'r',t, P,'b',t,G,'g')




import numpy as np
import scipy
from scipy import integrate
import matplotlib.pyplot as plt 

def rate (V, t):
    H, P, G = V
    n = 3; k1 = 0.2; k2 = 0.2; k3 = 0.2; n1= 5
    
    dHdt = 1/(1+G**n)-k1*H
    dPdt = H - k2*P
    dGdt = P - k3*G
    
    return np.array([dHdt, dPdt, dGdt])

import math 
res = math.sqrt(2)

tmin = 0
tmax = 100
Nt = 2000

t = np.linspace(tmin, tmax, Nt+1)

V0 = np.array([0.1, 0.9, 2.0])
Ans = integrate.odeint(rate, V0, t)
V = Ans.T
H, P, G = V
plt.figure(2)
plt.plot(t,H,'r',t, P,'b',t,G,'g')



import numpy as np
import scipy
from scipy import integrate
import matplotlib.pyplot as plt 

def rate (V, t):
    H, P, G = V
    n = 5; k1 = 0.2; k2 = 0.2; k3 = 0.2 
    
    dHdt = 1/(1+G**n)-k1*H
    dPdt = H - k2*P
    dGdt = P - k3*G
    
    return np.array([dHdt, dPdt, dGdt])

import math 
res = math.sqrt(2)

tmin = 0
tmax = 100
Nt = 2000

t = np.linspace(tmin, tmax, Nt+1)

V0 = np.array([0.1, 0.9, 2.0])
Ans = integrate.odeint(rate, V0, t)
V = Ans.T
H, P, G = V
plt.figure(3)
plt.plot(t,H,'r',t, P,'b',t,G,'g')
