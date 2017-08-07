# -*- coding: utf-8 -*-
"""
Created on Fri Aug 04 10:52:06 2017
@author: zWX460130
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from numpy.linalg import *
import os

###### Flow definition #########################################################
maxIter = 10001 # Total number of time iterations.
Re      = 100.0  # Reynolds number.
nx = 200; 
ny = 50; 
Q = 9 # Lattice dimensions and populations.
uLB     = 0.1 # 4                       # Velocity in lattice units.
nulb    = uLB * nx / Re; 
omega = 1.0 / (3. * nulb + 0.5); # Relaxation parameter.

# circle center
cx = nx/10; cy=ny/2; r= ny / 10; 


# Directions in D2Q9 model
c = np.array([[0,0] , [1,0] , [0,-1] , [-1, 0] , [0, 1], [1, -1], [-1, -1], [-1, 1], [1, 1]])
norm_dir = [i for i in range(0,Q)]
opposite_dir = [c.tolist().index((-c[i]).tolist()) for i in range(Q)]

rw = np.arange(Q)[np.asarray([ci[0]>0  for ci in c])] # Unknown on right wall.
mw = np.arange(Q)[np.asarray([ci[0]==0 for ci in c])] # Vertical middle.
lw = np.arange(Q)[np.asarray([ci[0]<0  for ci in c])] # Unknown on left wall.
print(lw)
print(rw)
print(mw)
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

# ----------------- VISUALISATION

def streamline_plot(time, v, x_min, x_max, y_min, y_max, name):
    fig = plt.figure(figsize = (20, 8))
    ax = plt.subplot(111)
    ax.set_xlabel(u'X-координата.')
    ax.set_ylabel(u'Y-координата.')
#    plt.title(u'Линии тока, Re=' + str(Re) + '.')
    
    x = np.array([i for i in range(x_min, x_max)])
    y = np.array([i for i in range(y_min, y_max)])
    
    vx = v[0,y_min:y_max,x_min:x_max]
    vy = v[1,y_min:y_max,x_min:x_max]
    
    strm = plt.streamplot(x, y, vx, vy, density = [0.6, 1], linewidth=2, color = '#1f77b4')  # color = 'k', linewidth=2, minlength = 0.1) #cmap=plt.cm.autumn
    #plt.colorbar(strm.lines)
    
#    obstacle = plt.Circle((cx, cy), r, color='k')
#    ax.add_artist(obstacle)
    
    out_name = '_'.join([name, 'strml', 're' + str(Re), 't' + str(time)])
    out_name += '.png'
    plt.savefig(out_name)
    plt.clf()

def contour_plot(time, val, val_name):
    fig = plt.figure(figsize = (20, 8))
    ax = plt.subplot(111)
    ax.set_xlabel(u'X-координата.')
    ax.set_ylabel(u'Y-координата.')
#    plt.title(u'Поле средних скоростей, Re=' + str(Re) + '.')
    
    x = np.array([i for i in range(0, nx)])
    y = np.array([i for i in range(0, ny)])
    
    data = []
    if len(val.shape) == 3:
        data = np.sqrt(val[0]**2 + val[1]**2)
    elif len(val.shape) == 2:
        data = val
    else:
        print('Dimension error in contour plot!')
#    lin = np.linspace(-0.02, 0.02, 25, endpoint=True)
    cont = plt.contourf(x, y, data, cmap=plt.cm.jet) # np.rotate() ''', norm=plt.Normalize(vmax=0.02, vmin=-0.02)'''
    b = plt.colorbar(cont, orientation='vertical')
    
#    obstacle = plt.Circle((cx, cy), r, color='k')
#    ax.add_artist(obstacle)
    
    out_name = '_'.join([val_name, 'cont', 're' + str(Re), 't' + str(time)])
    out_name += '.png'
    plt.savefig(out_name)

def pois_prof_comp(time, val):
      fig1 = plt.figure(figsize = (8, 8))
      ax = plt.subplot(111)
      ax.set_xlabel(u'Y-координата.')
      ax.set_ylabel(u'Величина скорости вдоль оси Ox.')
      plt.grid(True)
      
      y = [j + 1 for j in range(0, ny)]
      print(len(y))
      
      rad = (ny - 1.) / 2.
      print(rad)
      v_teoretical = [0.037 / (4. * nulb * nx) * (rad**2 - np.abs(i - rad)**2) for i in range(0, ny)]
      print(len(v_teoretical))
      v_analitical = val[0, 0:ny, nx/2]
      
      ax.plot(y, v_analitical, label = u'Экспериментальная зависимость',color = '#ff7f0e', linewidth=2., marker='<' )
      ax.plot(y, v_teoretical, label = u'Теоритическая зависимость.', color = '#1f77b4', linewidth=2., marker='o' )
      
      plt.legend()
      plt.savefig("profile_"+str(time).zfill(4)+".png")
      plt.clf()


def cavity_flow(time, val):
    fig1 = plt.figure(figsize = (8,8))
    ax = plt.subplot(111)
    # Veocity X
    ax.set_xlabel(u'Y-координата.')
    ax.set_ylabel(u'Величина скорости вдоль оси X.')
    plt.grid(True)
    
    y = [i+1 for i in range(0, ny)]
    vx = val[0, 0:ny, nx/2]
    vx = vx * -10
    ax.plot(y, vx, label=u'Результаты моделирования.', color='#1f77b4', linewidth=2., marker='o')
    
#    U. Ghia, K.N. Ghia, C.T.Shin. DATA RE 100
    y_ch = [1,9,10,11,13,21,30,31,65,104,111,117,122,123,124,125,129]
#    vx_ch = [0., 0.09233, 0.10091, 0.10890, 0.12317, 0.16077, 0.17507, 0.17527, 0.05454, -0.24533, -0.22445, -0.16914, -0.10313, -0.08864, -0.07391, -0.05906, 0.0]
#    RE 400
#    vx_ch = [0., 0.18360, 0.19713, 0.20920, 0.22965, 0.28124, 0.30203, 0.30174, 0.05186, -0.38598, -0.44993, -0.34827, -0.22847, -0.19254, -0.15663, -0.12146, 0.]
    vx_ch = [0., 0.27485, 0.29012, 0.30353, 0.32627, 0.37095, 0.33075, 0.32235, 0.02526, -0.31966, -0.42665, -0.51550, -0.39188, -0.33714, -0.27669, -0.21388, 0.]
    ax.plot(y_ch, vx_ch, 'kx', label='U. Ghia, K.N. Ghia, C.T.Shin.', markersize=10. )
    
    plt.legend()
    name = '_'.join(['cavity', 'vx', 're' + str(Re), 't' + str(time)])
    name +='.png'
    plt.savefig(name)
    plt.clf()
    
    # Velocity Y
    fig1 = plt.figure(figsize = (8,8))
    ax = plt.subplot(111)
    ax.set_xlabel(u'Величина скорости вдоль оси Y.')
    ax.set_ylabel(u'X-координата.')
    plt.grid(True)
    x = [nx - i for i in range(0, nx)]
    vy = val[1, ny/2]
    for i in range(0, nx/2):
        vy[i], vy[nx-1-i] = vy[nx-1-i], vy[i]
    vy = vy *10
    ax.plot(vy, y, label=u'Результаты моделирования.', color='#ff7f0e', linewidth=2., marker='<')
    
#    U. Ghia, K.N. Ghia, C.T.Shin. DATA RE 100
    x_ch = [1,8,9,10,14,23,37,59,65,80,95,110,123,124,125,126,129]
#    vy_ch = [0.0, -0.03717,-0.04192,-0.04775,-0.06434, -0.10150,-0.15662, -0.21090,-0.20581, -0.13641, 0.00332, 0.23151, 0.68717, 0.73722, 0.78871, 0.84123, 1.]
#    RE 400
#    vy_ch = [0., -0.08186, -0.09266, -0.10338, -0.14612, -0.24299, -0.32726, -0.17119, -0.11477, 0.002135, 0.16256, 0.29093, 0.55892, 0.61756, 0.68439, 0.75837, 1.]
    vy_ch = [0., -0.18109, -0.20196, -0.22220, -0.29730, -0.38289, -0.27805, -0.10648, -0.06080, 0.05702, 0.18719, 0.33304, 0.46604, 0.51117, 0.57492, 0.65928, 1.]
    ax.plot(vy_ch, x_ch, 'kx', label='U. Ghia, K.N. Ghia, C.T.Shin.', markersize=10.)
    
    plt.legend()
    name = '_'.join(['cavity', 'vy', 're' + str(Re), 't' + str(time)])
    name +='.png'
    plt.savefig(name)
    plt.clf()
    

# ---------------------- MAIN CAVITY CODE -----------------

#if __name__ == "__main__":
#    if not os.path.exists('cav_flow'):
#        os.makedirs('cav_flow')
#        
#    rho = np.ones((ny, nx))
#    v = np.zeros((2, ny, nx))
#    v[1, 1:ny-1, 0] = uLB
#    
#    feq = np.zeros((Q, ny, nx))   
#    feq = feq_calc(rho, v)
#    fin = feq.copy()
#       
#    bb_top()
#    bb_bot()
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
#        if (time%2000==0 and time > 0):
#            plt.clf()
#            x =  np.array([j for j in range(0, nx)])
#            y =  np.array([j for j in range(0, ny)])
#            
#            contour_plot(time, v, 'vel')
#            streamline_plot(time, v, 0, nx, 0, ny, 'vel')
#            streamline_plot(time, v, nx/2, nx, ny/2, ny, 'vel_top')
#            streamline_plot(time, v, nx/2, nx, 0, ny/2, 'vel_bot')
##            cavity_flow(time, v)
#

            
#           --------- FLOW AROUND CYLINDER

if __name__ == "__main__":
    rho = np.ones((ny, nx))
    v = np.zeros((2, ny, nx))
    # Velocity on Left wall
    v_r = np.array([uLB for i in range(0, ny-1)]) #[ uLB * np.sin(np.pi * o / (ny - 1) ) for o in range(0, ny)]
    v[0, 1:ny-1, 0] = v_r[1:ny-1]
    # Dencity Left wall
#    rho[1:ny-1, 0] = 1.001
    # Dencity Right wall
#    rho[1:ny-1, -1] = 1.0
    
    feq = np.zeros((Q, ny, nx))
    feq = feq_calc(rho, v)
    fin = feq.copy()
    
    
     # Set Bounce-Back bounaries    
    bb_top()
    bb_bot()
#    bb_left()
#    bb_right()
#    obst()
    
    for time in range(0, maxIter):
        rho = calc_rho(fin)
        v = calc_v(fin) / rho
        print(str(time) + ' rho = '+ str(np.sum(rho)))    
        
        # Left boundary velocities
        v[0, 1:ny-1, 0] = v_r[1:ny-1]
        v[1, 1:ny-1, 0] = 0.0
        rho[1:ny-1, 0] = 1./ (1.-v[0, 1:ny-1, 0]) * (sumpop(fin[mw,1:ny-1, 0]) +2.*sumpop(fin[lw,1:ny-1, 0]))
        # Left boundary density
#        rho[1:ny-1, 0] = 1.001
#        v[0,1:ny-1,0] = 1. - (sumpop(fin[mw,1:ny-1, 0]) +2.*sumpop(fin[lw,1:ny-1, 0])) / rho[1:ny-1, 0]
        #Right boundary density
#        rho[1:ny-1,-1] = 1.0
#        v[0,1:ny-1,-1] = (sumpop(fin[mw,1:ny-1, -1]) +2.*sumpop(fin[rw,1:ny-1, -1])) / rho[1:ny-1, -1] - 1.
        
        # Bounce-back updates
#        v[0, boundary] = 0.0
#        v[1, boundary] = 0.0
        
        # Calculation of equilibrium distribution function for Zo-He BCs
        feq = feq_calc(rho, v)
        # Left wall
        fin[rw,1:ny-1, 0] = fin[lw,1:ny-1, 0] + feq[rw,1:ny-1, 0] - fin[lw,1:ny-1, 0]
        #Right wall
#        fin[3,1:ny-1,-1] = fin[1,1:ny-1,-1] - 2./3. * rho[1:ny-1,-1] * v[0,1:ny-1,-1]
#        temp = fin[4,1:ny-1,-1]-fin[2,1:ny-1,-1] / 2.
#        fin[6,1:ny-1,-1] = fin[8,1:ny-1,-1] + temp - rho[1:ny-1,-1] * v[0,1:ny-1,-1] / 6.
#        fin[7,1:ny-1,-1] = fin[5,1:ny-1,-1] - temp - rho[1:ny-1,-1] * v[0,1:ny-1,-1] / 6.
#        fin[lw,1:ny-1, -1] = fin[rw,1:ny-1,-1] + feq[lw,1:ny-1,-1] - fin[rw,1:ny-1,-1]
        
        # Collision
        fout = fin - omega * (fin - feq)
        
        # Bounce-back BCs
        for q in range(Q): 
            fout[norm_dir[q], boundary] = fin[opposite_dir[q],boundary]
        # Streaming
        for q in range(Q):
            fin[q,:,:] = np.roll(np.roll(fout[q,:,:],c[q,0],axis=1),c[q,1],axis=0)
        # Write information
        if (time%1000==0 and time > 0):

            contour_plot(time, v ,'vel')
            streamline_plot(time, v, 0, nx / 4, 0, ny, 'vel')
            pois_prof_comp(time, v)