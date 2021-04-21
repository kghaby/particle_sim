#!/usr/bin/env python3 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
import pandas as pd
#from sys import exit

###plot
def update_lines(num, data, line):
    # NOTE: there is no .set_data() for 3 dim data...
    line.set_data(data[0:2, :num])    
    line.set_3d_properties(data[2, :num])    
    return line

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)

# Reading the data from a CSV file using pandas
repo = pd.read_csv('traj.dat',sep=' ',header=None)
data = np.array((repo[0].values, repo[1].values, repo[2].values))
#print (data.shape[1])
#exit()

# Creating fifty line objects.
# NOTE: Can't pass empty arrays into 3d version of plot()
limit = 16.
line = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1])[0]

# Setting the axes properties
ax.set_xlim3d([-limit, limit])
ax.set_xlabel('X')
ax.set_ylim3d([-limit, limit])
ax.set_ylabel('Y')
ax.set_zlim3d([-limit, limit])
ax.set_zlabel('Z')
ax.set_title('Brownian motion from a 3D potential')

# Creating the Animation object
line_ani = animation.FuncAnimation(fig, update_lines, data.shape[1], fargs=(data, line), interval=1, blit=False)

plt.show()


