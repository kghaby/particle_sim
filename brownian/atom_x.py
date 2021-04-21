#!/usr/bin/env python3 
import numpy as np


## electron
px_e=.000001
py_e=0
pz_e=0
vx_e=0
vy_e=0
vz_e=0
ax_e=0
ay_e=0
az_e=0

## nucleus
px_n=0
py_n=0
pz_n=0
vx_n=0
vy_n=0
vz_n=0
ax_n=0
ay_n=0
az_n=0

### values
NSteps=10000
dt=.1
m=1 #
T=2 #speed/energy
C=0.01
### potential 
#U=px_e**2+py_e**2+pz_e**2
#U=(px_e**2/m+np.exp(-px_e**2/T)/m**2)
	#Fx = -C*( (2*px_e/m-2*px_e*np.exp(-px_e**2/T)/m**2*T) )
#U=m*px_e**2+np.cos(T*px_e)
	#Fx=-(2*m*px_e-T*np.sin(T*px_e))
#F=(dU/dx + dU/dy + dU/dz)
### simulate
for i in range(0,NSteps):
	R1x = np.random.rand(1)
	R2x = np.random.rand(1)
	Fx = -C*( 2*m*px_e - T*np.sin(T*px_e) )
	print(i, "\t",str(px_e).lstrip('[').rstrip(']'),"\t",str(vx_e).lstrip('[').rstrip(']'),str(Fx).lstrip('[').rstrip(']'))
	gx = np.sqrt(-2*np.log(R1x)) * np.cos(2*np.pi*R2x)
	px_e = px_e + vx_e * dt - 0.5*Fx/m * dt**2
	vx_e = vx_e + Fx/m * dt #+ gx

		
	
	
	
