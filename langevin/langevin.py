#!/usr/bin/env python3.9
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import gaussian_kde
from sklearn.preprocessing import normalize
import sys

orig_stdout = sys.stdout

###defining stuff
T=300 #K
kT=0.5915 #kcal/mol
D=3 #A^2/fs 
Nsteps=1000000 #fs => 1 ns
dt=1 #timestep is 1 fs
y=kT/D
m=10 #amu
K=1*kT #A^-2
#U=
#F=-Kx


###reset 
if os.path.exists("traj.dat"):
  os.remove("traj.dat")
else:
  print("traj.dat does not exist")
x=0
v=0

###langevin loop
file=open("traj.dat","w")
sys.stdout = file
for b in range(1,1000000):
        R1 = np.random.rand(1)
        R2 = np.random.rand(1)
        g = np.sqrt((2*kT*y*0.001)/m**2) * np.sqrt(-2*np.log(R1)) * np.cos(2*np.pi*R2)
        F = -1*K*x
        x = x+v*0.001
        v = v+F/m*0.001-(y/m)*v*0.001+g
        print(str(x).lstrip('[').rstrip(']'),"\t",str(v).lstrip('[').rstrip(']'))
        #x_array=np.array([x])
        #v_array=np.array([v])
        #xv_array=np.array([x_array, v_array])
        #xv_array = xv_array.T
        #np.savetxt(file, xv_array, fmt=['%f','%f'])
sys.stdout = orig_stdout
file.close()

####plot x over time
xv_list=np.loadtxt("traj.dat")
plt.plot(xv_list[:,0],label='x')
plt.plot(xv_list[:,1],label='v')
plt.xlabel("Time (fs)")
#plt.ylabel("X")
plt.xlim([0, 25000])
plt.legend()
plt.show()

