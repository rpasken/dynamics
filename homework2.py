#!/bin/python3
#
# Get the necessary libraries
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Tkagg")


#
# Ask for information about the wave structure
#

base_geopot = float(input("Base geoptential height? "))
mean_u = float(input("Mean u component? "))
mean_v = float(input("Mean v component? "))
cor = float(input("Beta? "))
Lx = int(input("X wavelength? "))
Ly = int(input("Y wavelength? "))

#
# Scaling factor to ake the plots look nice
#

A = 1e3

#
# Build a range of x/y values for the function
#

xx = range(-3000,3000,300)
yy = range(-1000,1000,100)

#
# Make the grid
#

x,y = np.meshgrid(xx,yy)

#
# To make life easier later get the length of the side of the grid
#

x_len = len(xx)
y_len = len(yy)

#
# Build empty lists (arrays)
#

geopot = np.empty((x_len,y_len))
vort = np.empty((x_len,y_len))
ugeo = np.empty((x_len,y_len))
vgeo = np.empty((x_len,y_len))
vort_advect_x = np.empty((x_len,y_len))
vort_advect_y = np.empty((x_len,y_len))
#
# Compute the  wave numbers
#

K = (2*np.pi) / Lx
L = (2*np.pi) / Ly

#
# Grid through all the grid points. Yes there are faster more
# obtuse python ways
#

for i in range(0,len(xx)):
    for j in range(len(yy)):
        geopot[i,j] = base_geopot - mean_u * cor * y[i,j] * A + mean_v * cor /K * A * np.sin(K*x[i,j]) * np.cos(L*y[i,j])
        vort[i,j] =  mean_v*((K**2 + L**2)/K) * np.sin(K*x[i,j]) * np.cos(L*y[i,j]) + cor
        ugeo[i,j] = A*(mean_u + mean_v * L/K * np.sin(K*x[i,j])*np.sin(L*y[i,j]))
        vgeo[i,j] = A*(mean_v * np.cos(K*x[i,j])*np.cos(L*y[i,j]))
        vort_advect_x[i,j] = ugeo[i,j]*mean_v*(K**2+L**2)*np.cos(K*x[i,j])*np.cos(L*y[i,j])
        vort_advect_y[i,j] = vgeo[i,j]*mean_v*(K**2+L**2)*np.cos(K*x[i,j])*np.cos(L*y[i,j])
        
        
fig = plt.figure(figsize=(8, 6))
plt.contour(x,y,geopot/9.8)
#b = plt.contourf(x,y,vort)
#plt.colorbar(b)
#plt.quiver(x,y,ugeo,vgeo)
plt.quiver(x,y,vort_advect_x,vort_advect_y)
plt.show()
