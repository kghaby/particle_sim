#!/usr/bin/env python3 
import numpy as np

## electron
px_e=0
py_e=0
pz_e=0
#vx_e=0
#vy_e=0
#vz_e=0
#ax_e=0
#ay_e=0
#az_e=0
Fx=0
Fy=0
Fz=0

## nucleus
#px_n=0
#py_n=0
#pz_n=0
#vx_n=0
#vy_n=0
#vz_n=0
#ax_n=0
#ay_n=0
#az_n=0

#T=300 #K
kT=0.5915 #@300 K; kcal/mol
D=0.1 #A^2/ps
NSteps=1000000
dt=1 #ps
#m=1.007 #proton amu
#C=.0001
### potential/force function 
	#F = (dU/dx + dU/dy + dU/dz)
	#  = -1*(0.004*px_e**3-0.1*px_e+0.03) 
		#worked well. represents basic form i think
	#  = C*-2*px_e 
		#harmonic 
	#  = -(2*m*px_e-T*np.sin(T*px_e)) 
	#  = -C*( (2*px_e/m-2*px_e*np.exp(-px_e**2/T)/m**2*T) )
		#diverged 
	#  = -1*(0.0004*px_e**3 - 0.01*px_e + (2*(px_e+5))*np.exp(-1*(px_e+5)**2) + (2*(px_e-5))*np.exp(-1*(px_e-5)**2)
	#  = -1*(0.0004*px_e**3 - 0.01*px_e + (4 *(px_e+5))*np.exp(-1*(px_e+5)**2) + (6*(px_e-5))*np.exp(-1*(px_e-5)**2))
		#uneven holes 
## values

### simulate
for i in range(0,NSteps):
    R1x = np.random.rand(1)
    R2x = np.random.rand(1)
    R1y = np.random.rand(1)
    R2y = np.random.rand(1)
    R1z = np.random.rand(1)
    R2z = np.random.rand(1)
    #print(i, "\t",str(px_e).lstrip('[').rstrip(']'),"\t",str(py_e).lstrip('[').rstrip(']'),"\t",str(pz_e).lstrip('[').rstrip(']'))
    #print(str(px_e).lstrip('[').rstrip(']'),"\t",str(py_e).lstrip('[').rstrip(']'),"\t",str(pz_e).lstrip('[').rstrip(']')) #without time
    print(str(px_e).lstrip('[').rstrip(']'),str(py_e).lstrip('[').rstrip(']'),str(pz_e).lstrip('[').rstrip(']')) #without time
    gx = np.sqrt(2*D*dt) * np.sqrt(-2*np.log(R1x)) * np.cos(2*np.pi*R2x)
    gy = np.sqrt(2*D*dt) * np.sqrt(-2*np.log(R1y)) * np.cos(2*np.pi*R2y)
    gz = np.sqrt(2*D*dt) * np.sqrt(-2*np.log(R1z)) * np.cos(2*np.pi*R2z)
    px_e = px_e + (D*Fx)/kT * dt + gx
    py_e = py_e + (D*Fy)/kT * dt + gy
    pz_e = pz_e + (D*Fz)/kT * dt + gz
    Fx = -1*(0.004*px_e**3-0.06*px_e-28*px_e*np.exp(-1.4*px_e**2))
    Fy = -1*(0.004*py_e**3+1.2*py_e)
    Fz = -1*(0.004*pz_e**3+1.2*pz_e)

