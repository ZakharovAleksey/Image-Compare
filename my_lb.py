# -*- coding: utf-8 -*-
"""
Created on Fri Aug 04 10:52:06 2017

@author: zWX460130
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from numpy.linalg import *
###### Flow definition #########################################################
maxIter = 200001 # Total number of time iterations.
Re      = 220.0  # Reynolds number.
nx = 200; 
ny = 70; 
Q = 9 # Lattice dimensions and populations.
uLB     = 0.01 # 4                       # Velocity in lattice units.
nulb    = uLB * nx / Re; 
omega = 1.0 / (3. * nulb + 0.5); # Relaxation parameter.

# circle center
cx = nx/4; cy=ny/2; r=ny/9; 


# Directions in D2Q9 model
c = np.array([[0,0] , [1,0] , [0,-1] , [-1, 0] , [0, 1], [1, -1], [-1, -1], [-1, 1], [1, 1]])
norm_dir = [i for i in range(0,Q)]
opposite_dir = [c.tolist().index((-c[i]).tolist()) for i in range(Q)]

rw = np.arange(Q)[np.asarray([ci[0]>0  for ci in c])] # Unknown on right wall.
mw = np.arange(Q)[np.asarray([ci[0]==0 for ci in c])] # Vertical middle.
lw = np.arange(Q)[np.asarray([ci[0]<0  for ci in c])] # Unknown on left wall.

# Weights vector
w = 1.0 / 36.0 * np.ones(Q)
w[0] = 4.0 / 9.0
w[1:5] = 1.0 / 9.0


boundary = np.fromfunction(lambda x, y : x == -1, (ny, nx), dtype=int)

def bb_top():
    boundary[0].fill(True)
    
def bb_bot():
    boundary[ny-1].fill(True)
    
def bb_left():
    boundary[1:ny-1, 0].fill(True)
    
def bb_right():
    boundary[1:ny-1, nx-1].fill(True)


def obst():
    for x in range(nx):
        for y in range(ny):
            if (x-cx)**2+(y-cy)**2<r**2:
                boundary[y,x] = True
    
 
    
########## func
def feq_calc(rho, v):
    vSq = v[0]**2 + v[1]**2
    feq = np.zeros((Q, ny, nx))
    
    for q in range(0, Q):
        ve = c[q][0] * v[0] + c[q][1] * v[1]
        feq[q] = w[q] * rho * (1.0 + 3.0*ve + 4.5*ve**2 - 1.5*vSq)
    return feq

def calc_rho(fin):
    R = np.zeros((ny, nx))
    for q in range(Q):
        R += fin[q]
    return R

def calc_v(fin):
    V = np.zeros((2, ny, nx))
    for j in range(2):
        for q in range(Q):
            V[j] += c[q, j] * fin[q]
    return V

sumpop = lambda fin: np.sum(fin,axis=0)

# ---------------------- MAIN CAVITY CODE -----------------

#if __name__ == "__main__":
#    fig = plt.figure(figsize = (10, 5))
#    plt.xlabel(u'X-координата')
#    plt.ylabel(u'Y-координата')
#    plt.title(u'Линии уровня Re=' + str(Re) + '.')
#    
#    rho = np.ones((ny, nx))
#    v = np.zeros((2, ny, nx))
#    v[1, 1:ny-1, 0] = uLB
#    
#    feq = np.zeros((Q, ny, nx))
#    
#    
#    feq = feq_calc(rho, v)
#    fin = feq.copy()
#    
#    
#    bb_top()
#    bb_bot()
#    #bb_left()
#    bb_right()
#    
#    for time in range(0, maxIter):
#        rho = calc_rho(fin)
#        v = calc_v(fin) / rho
#        print(str(time) + ' rho = '+ str(np.sum(rho)))    
#        
#        v[1, 1:ny-1, 0] = uLB
#        v[0, 1:ny-1, 0] = 0.0
#        v[0, boundary] = 0.0
#        v[1, boundary] = 0.0
#        rho[1:ny-1, 0] = 1./ (1.-v[0, 1:ny-1, 0]) * (sumpop(fin[mw,1:ny-1, 0]) +2.*sumpop(fin[lw,1:ny-1, 0]))
#        
#        feq = feq_calc(rho, v)
#        fin[rw,1:ny-1, 0] = fin[lw,1:ny-1, 0] + feq[rw,1:ny-1, 0] - fin[lw,1:ny-1, 0]            
#        
#        fout = fin - omega * (fin - feq)
#        
#        for q in range(Q): 
#            fout[norm_dir[q], boundary] = fin[opposite_dir[q],boundary]
#        
#        for q in range(Q):
#            fin[q,:,:] = np.roll(np.roll(fout[q,:,:],c[q,0],axis=1),c[q,1],axis=0)
#        
#        if (time%500==0 and time > 0):
#            plt.clf()
#            x =  np.array([j for j in range(0, nx)])
#            y =  np.array([j for j in range(0, ny)])
#
#            strm = plt.streamplot(x, y, v[0],v[1], color= np.sqrt(v[0]**2+v[1]**2), linewidth=2, cmap=plt.cm.autumn, minlength = 0.1)
#            plt.colorbar(strm.lines)
#            plt.savefig("vel."+str(time).zfill(4)+".png")
#            plt.clf()
#            
#            cont = plt.contourf(x,y, np.sqrt(v[0]**2 + v[1]**2))
#            b = plt.colorbar(cont, orientation='vertical')
#            plt.savefig("cont."+str(time).zfill(4)+".png")
#            plt.clf()
            
#            x =  np.array([j for j in range(nx - 50, nx)])
#            y =  np.array([j for j in range(ny - 40, ny)])
#            strm = plt.streamplot(x, y, v[0, ny-40:ny, nx-50:nx],v[1, ny-40:ny, nx-50:nx],color='k', linewidth=2, minlength = 0.1)
#            plt.savefig("top_vel."+str(time).zfill(4)+".png")
            
#           --------- MAIN FLOW CODE-----

if __name__ == "__main__":
    fig = plt.figure(figsize = (20, 10))
    plt.xlabel(u'X-координата')
    plt.ylabel(u'Y-координата')
    plt.title(u'Линии уровня Re=' + str(Re) + '.')
    
    rho = np.ones((ny, nx))
    v = np.zeros((2, ny, nx))
    v_r = [ uLB * np.sin(np.pi * o / (ny - 1) ) for o in range(0, ny)]
    
    v[0, 1:ny-1, 0] = v_r[1:ny-1]# uLB
    
    
    feq = np.zeros((Q, ny, nx))
    
    
    feq = feq_calc(rho, v)
    fin = feq.copy()
    
    
    #bb_top()
    #bb_bot()
    #bb_left()
    #bb_right()
    
    obst()
    
    for time in range(0, maxIter):
        #fin[lw,:,-1] = fin[lw,:,-2] # Right wall: outflow condition.
        rho = calc_rho(fin)
        v = calc_v(fin) / rho
        print(str(time) + ' rho = '+ str(np.sum(rho)))    
        
        v[0, 1:ny-1, 0] = v_r[1:ny-1]# uLB
        v[1, 1:ny-1, 0] = 0.0
        v[0, boundary] = 0.0
        v[1, boundary] = 0.0
        rho[1:ny-1, 0] = 1./ (1.-v[0, 1:ny-1, 0]) * (sumpop(fin[mw,1:ny-1, 0]) +2.*sumpop(fin[lw,1:ny-1, 0]))
        
        feq = feq_calc(rho, v)
        fin[rw,1:ny-1, 0] = fin[lw,1:ny-1, 0] + feq[rw,1:ny-1, 0] - fin[lw,1:ny-1, 0]            
        
        fout = fin - omega * (fin - feq)
        
        for q in range(Q): 
            fout[norm_dir[q], boundary] = fin[opposite_dir[q],boundary]
        
        for q in range(Q):
            fin[q,:,:] = np.roll(np.roll(fout[q,:,:],c[q,0],axis=1),c[q,1],axis=0)
        
        if (time%1000==0 and time > 0):
            fig1 = plt.figure(figsize = (20, 10))
            plt.clf()
            x =  np.array([j for j in range(0, nx)])
            y =  np.array([j for j in range(0, ny)])

            strm = plt.streamplot(x, y, v[0],v[1], color= np.sqrt(v[0]**2+v[1]**2), linewidth=2, cmap=plt.cm.autumn, minlength = 0.1)
            plt.colorbar(strm.lines)
            plt.savefig("vel."+str(time).zfill(4)+".png")
            plt.clf()
            
            cont = plt.contourf(x,y, np.sqrt(v[0]**2 + v[1]**2))
            b = plt.colorbar(cont, orientation='vertical')
            plt.savefig("cont."+str(time).zfill(4)+".png")
            plt.clf()
            
#            fig1 = plt.figure(figsize = (10, 10))
#            ax = plt.subplot(111)
#            y =  np.array([j for j in range(0, ny)])
#            v_r = [0.002 / (4. * nulb * nx) * ((ny/2.)**2 - np.abs(i - ny/2)**2) for i in range(0, ny)]
#            ax.plot(y, v[0, 0:ny, nx/2], 'bo')
#            ax.plot(y, v_r, 'r--')
#            plt.savefig("profile_"+str(time).zfill(4)+".png")
#            plt.clf()
            