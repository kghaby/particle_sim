#!/usr/bin/env python3.9
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import gaussian_kde
from sklearn.preprocessing import normalize

###resetting objects

###defining stuff
T=300 #K
kT=0.5915 #kcal/mol
D=0.1 #A^2/ps
Nsteps=1000 #ps => 1 ns
dt=1 #ps
#U=0.0025*x**4-0.5*x**2+0.2*x
#K=1 #for checking
#C=kT/K #for checking

###reset x 
if os.path.exists("x.txt"):
  os.remove("x.txt")
else:
  print("x.txt does not exist")
x=0

###brownian loop
file=open("x.txt","a")

for t in range(1,Nsteps*dt+1):
        R1 = np.random.rand(1)
        R2 = np.random.rand(1)
        g = np.sqrt(2*D*dt) * np.sqrt(-2*np.log(R1)) * np.cos(2*np.pi*R2)
        #F = K*x #for checking
        F = 0.01*x**3-x+0.2
        x = x-((D*F)/kT)*dt+g
        x_array=np.array([x])
        np.savetxt(file, x_array, fmt='%f')
        
file.close()

#### plot x over time
x_list=np.loadtxt("x.txt")
plt.plot(x_list)
plt.xlabel("Time (ps)")
plt.ylabel("X")
plt.show()

###for checking functionality of code with harmonic potential. avg_x^2=kT/K is expected
#avg_x=np.average(x_list)
#avg_x2=np.average(x_list**2)
#print(
#    "avg_x = ",avg_x,"\n"
#    "avg_x^2 = ",avg_x2,"\n"
#    "kT/K = ", C
#     )	

####

##replicates of brownian dynamics from above
if os.path.exists("x.txt"):
  os.remove("x.txt")
else:
  print("x.txt does not exist")

x=0
file=open("x.txt","a")
for i in range(1,1000):
    R1 = np.random.rand(1)
    R2 = np.random.rand(1)
    for t in range(1,Nsteps*dt):
            g = np.sqrt(2*D*dt) * np.sqrt(-2*np.log(R1)) * np.cos(2*np.pi*R2)
            F = 0.01*x**3-x+0.2
            x = x-((D*F)/kT)*dt+g
            x_array=np.array([x])
            np.savetxt(file, x_array, fmt='%f')
        
file.close()
x_list=np.loadtxt("x.txt")

#boltzmann dist
if os.path.exists("boltz.txt"):
  os.remove("boltz.txt")
else:
  print("boltz.txt does not exist")
file=open("boltz.txt","a")
b_array=np.linspace(-14,14,200)
for b in b_array:
        d=np.exp(-(0.0025*b**4-0.5*b**2+0.2*b)/kT)
        d_array=np.array([d])
        b_array2=np.array([b])
        bd_array=np.array([b_array2, d_array])
        bd_array = bd_array.T
        np.savetxt(file, bd_array, fmt=['%f','%f'])
file.close()
boltz_list=np.loadtxt("boltz.txt")
norm_boltz_list=boltz_list / np.linalg.norm(boltz_list)
plt.plot(boltz_list[:,0], norm_boltz_list[:,1],label='Boltzmann Distribution')
plt.yscale("log")
density_x = gaussian_kde(x_list)
xs = np.linspace(-17,17,200)
density_x.covariance_factor = lambda : 0.05
density_x._compute_covariance()
plt.plot(xs,density_x(xs),label='Brownian Dyn. Prob. Density')
#plt.hist(x_list, bins=100, density=True, label='Probability Density')
#plt.legend()
#plt.show()

##monte carlo
Nsims=1
Nsteps=1000
mc=np.random.rand(Nsims,Nsteps)
prob_mat=np.zeros((Nsims,Nsteps))
for i in range(1,Nsteps*dt):
    R=np.random.rand(1)
    if ((0.0025*mc[:,i]**4-0.5*mc[:,i]**2+0.2*mc[:,i]) <
        (0.0025*mc[:,i-1]**4-0.5*mc[:,i-1]**2+0.2*mc[:,i-1])):    
            prob_mat[:,i] = -np.exp(-((0.0025*mc[:,i]**4-0.5*mc[:,i]**2+0.2*mc[:,i])-(0.0025*mc[:,i-1]**4-0.5*mc[:,i-1]**2+0.2*mc[:,i-1])/kT))
    else:
            if np.exp(-((0.0025*mc[:,i]**4-0.5*mc[:,i]**2+0.2*mc[:,i])-(0.0025*mc[:,i-1]**4-0.5*mc[:,i-1]**2+0.2*mc[:,i-1])/kT)) > R:
                prob_mat[:,i] = np.exp(-((0.0025*mc[:,i]**4-0.5*mc[:,i]**2+0.2*mc[:,i])-(0.0025*mc[:,i-1]**4-0.5*mc[:,i-1]**2+0.2*mc[:,i-1])/kT))
            else:
                prob_mat[:,i] = prob_mat[:,i-1]
prob_mat[:,]*=8
density_prob = gaussian_kde(prob_mat)
xs = np.linspace(-17,17,200)
density_prob.covariance_factor = lambda : 0.15
density_prob._compute_covariance()
plt.plot(xs,density_prob(xs),label='Monte Carlo Prob. Density (scaled x)')
#plt.hist(prob_mat.T, bins=100, density=True, label='Probability Density')
plt.xlabel("x")
plt.ylabel("Count")
plt.legend()
plt.show()





