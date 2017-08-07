# -*- coding: utf-8 -*-
"""
Created on Sat Aug 05 15:57:46 2017

@author: Alexey
"""

#!/usr/bin/python
# Copyright (C) 2013 FlowKit Ltd, Lausanne, Switzerland
# E-mail contact: contact@flowkit.com
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License, either
# version 3 of the License, or (at your option) any later version.

#
# 2D flow around a cylinder
#

from numpy import *; from numpy.linalg import *
import matplotlib.pyplot as plt; from matplotlib import cm
###### Flow definition #########################################################
maxIter = 200000 # Total number of time iterations.
Re      = 300.0  # Reynolds number.
nx = 200; ny = 50; ly=ny-1.0; q = 9 # Lattice dimensions and populations.
cx = nx/10; cy=ny/2; r= ny / 10;           # Coordinates of the cylinder.
uLB     = 0.04                       # Velocity in lattice units.
nulb    = uLB*r/Re; omega = 1.0 / (3.*nulb+0.5); # Relaxation parameter.

###### Lattice Constants #######################################################
c = array([(x,y) for x in [0,-1,1] for y in [0,-1,1]]) # Lattice velocities.
t = 1./36. * ones(q)                                   # Lattice weights.
t[asarray([norm(ci)<1.1 for ci in c])] = 1./9.; t[0] = 4./9.
noslip = [c.tolist().index((-c[i]).tolist()) for i in range(q)] 
i1 = arange(q)[asarray([ci[0]<0  for ci in c])] # Unknown on right wall.
i2 = arange(q)[asarray([ci[0]==0 for ci in c])] # Vertical middle.
i3 = arange(q)[asarray([ci[0]>0  for ci in c])] # Unknown on left wall.

###### Function Definitions ####################################################
sumpop = lambda fin: sum(fin,axis=0) # Helper function for density computation.
def equilibrium(rho,u):              # Equilibrium distribution function.
    cu   = 3.0 * dot(c,u.transpose(1,0,2))
    usqr = 3./2.*(u[0]**2+u[1]**2)
    feq = zeros((q,nx,ny))
    for i in range(q): feq[i,:,:] = rho*t[i]*(1.+cu[i]+0.5*cu[i]**2-usqr)
    return feq

###### Setup: cylindrical obstacle and velocity inlet with perturbation ########
obstacle = fromfunction(lambda x,y: (x-cx)**2+(y-cy)**2<r**2, (nx,ny))
vel = fromfunction(lambda d,x,y: (1-d)*uLB*(1.0+1e-4*sin(y/ly*2*pi)),(2,nx,ny))
feq = equilibrium(1.0,vel); fin = feq.copy()


def streamline_plot(time, v, x_min, x_max, y_min, y_max):
    print(v.shape)
    fig = plt.figure(figsize = (20, 5))
    ax = plt.subplot(111)
    ax.set_xlabel(u'X-координата.')
    ax.set_ylabel(u'Y-координата.')
#    plt.title(u'Линии тока, Re=' + str(Re) + '.')
    
    x = np.array([i for i in range(x_min,x_max)])
    y = np.array([i for i in range(y_min, y_max)])
    
    vx = v[0].transpose()
    vy = v[1].transpose()
    
    
#    print(v.shape)
    vx = vx[y_min:y_max,x_min:x_max]
    vy = vy[y_min:y_max,x_min:x_max]
    strm = plt.streamplot(x, y, vx, vy, density = [0.6, 1], linewidth=2, color = '#1f77b4')  # color = 'k', linewidth=2, minlength = 0.1) #cmap=plt.cm.autumn
    #plt.colorbar(strm.lines)
    
    obstacle = plt.Circle((cx, cy), r, color='k')
    ax.add_artist(obstacle)
    
    out_name = '_'.join(['strml', 're' + str(Re), 't' + str(time)])
    out_name += '.png'
    plt.savefig(out_name)
    plt.clf()

def contour_plot(time, val, val_name):
    fig = plt.figure(figsize = (20, 5))
    ax = plt.subplot(111)
    ax.set_xlabel(u'X-координата.')
    ax.set_ylabel(u'Y-координата.')
#    plt.title(u'Поле средних скоростей, Re=' + str(Re) + '.')
    
    x = np.array([i for i in range(0, nx)])
    y = np.array([i for i in range(0, ny)])
    print('hhoo ---')
    print(val.transpose().shape)
    cont = plt.contourf(x, y, val.transpose(), cmap=plt.cm.jet)
    b = plt.colorbar(cont, orientation='vertical')
    
    obstacle = plt.Circle((cx, cy), r, color='k')
    ax.add_artist(obstacle)
    
    out_name = '_'.join([val_name, 'cont', 're' + str(Re), 't' + str(time)])
    out_name += '.png'
    plt.savefig(out_name)

###### Main time loop ##########################################################
for time in range(maxIter):
    print(str(time) + ' rho = '+ str(np.sum(rho)))    
    fin[i1,-1,:] = fin[i1,-2,:] # Right wall: outflow condition.
    rho = sumpop(fin)           # Calculate macroscopic density and velocity.
    u = dot(c.transpose(), fin.transpose((1,0,2)))/rho

    u[:,0,:] =vel[:,0,:] # Left wall: compute density from known populations.
    rho[0,:] = 1./(1.-u[0,0,:]) * (sumpop(fin[i2,0,:])+2.*sumpop(fin[i1,0,:]))

    feq = equilibrium(rho,u) # Left wall: Zou/He boundary condition.
    fin[i3,0,:] = fin[i1,0,:] + feq[i3,0,:] - fin[i1,0,:]
    fout = fin - omega * (fin - feq)  # Collision step.
    for i in range(q): fout[i,obstacle] = fin[noslip[i],obstacle]
    for i in range(q): # Streaming step.
        fin[i,:,:] = roll(roll(fout[i,:,:],c[i,0],axis=0),c[i,1],axis=1)
 
    if (time%1000==0 and time > 0): # Visualization
        v = u[0]**2 + u[1]**2
        print(v.shape)
        contour_plot(time, v,'vel')
#            contour_plot(time, rho, 'rho')
        streamline_plot(time, u, 0, nx, 0, ny)
#        
#        plt.clf(); plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(),cmap=cm.Reds)
#        plt.savefig("vel."+str(time/100).zfill(4)+".png")
#
