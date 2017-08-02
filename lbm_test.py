# -*- coding: utf-8 -*-
"""
Created on Wed Aug 02 12:50:50 2017

@author: zWX460130
"""

import numpy as np
import matplotlib.pyplot as plt

Q = 9

x = 32
y = 32
time = 5000
tau = 1.0

rho = np.zeros((y, x), dtype = float)
vx = np.zeros((y, x), dtype = float)
vy = np.zeros((y, x), dtype = float)

f = np.zeros((Q, y, x), dtype = float)
feq = np.zeros((Q, y, x), dtype = float)

g = np.zeros((Q, y, x), dtype = float)

w = np.array([ 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 ])
ex = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1])
ey = np.array([0, 0, -1, 0, 1, -1, -1, 1, 1])


def set_up_rho(rho):
    for y in range (1, rho.shape[0] - 1):
        for x in range (1, rho.shape[1] - 1):
            rho[y][x] += 1.0
    

def f_eq_calc(feq, rho, vx, vy):
    
    for q in range(0, Q):
        v = ex[q] * vx + ey[q] * vy
        feq[q] = w[q] * rho * (1.0 + 3.0 * v + 4.5 * v * v - 1.5 * (vx * vx + vy* vy))


def g_calc():
    gx = 0.00001
    gy = 0.0
    g[0] = (1.0 - 0.5 / tau) * w[0] * (3.0 * (vx * gx + vy * gy) )
    
    g[1] = (1.0 - 0.5 / tau) * w[1] * (3.0 * ((1.0 - vx) * gx - vy * gy) )  + 9.0 * (vx * gx);
    g[2] = (1.0 - 0.5 / tau) * w[2] * (3.0 * (-vx * gx + ( 1.0 - vy) * gy) ) + 9.0 * (vy * gy);
    g[3] = (1.0 - 0.5 / tau) * w[3] * (3.0 * (( - 1.0 - vx) * gx - vy * gy) )  + 9.0 * (vx * gx);
    g[4] = (1.0 - 0.5 / tau) * w[4] * (3.0 * (-vx * gx + (- 1.0 - vy) * gy) ) + 9.0 * (vy * gy);
    
    g[5] = (1.0 - 0.5 / tau) * w[5] * (3.0 * ( (1.0 - vx) * gx + (1.0 - vy) * gy) ) + 9.0 * (vx + vy) * (gx + gy)
    g[6] = (1.0 - 0.5 / tau) * w[6] * (3.0 * ((-1.0 - vx) * gx + (1.0 - vy) * gy) )  + 9.0 * (vx - vy) * (gx - gy)
    g[7] = (1.0 - 0.5 / tau) * w[7] * (3.0 * (( -1.0 - vx) * gx + (-1.0 - vy) * gy) )  + 9.0 * (vx + vy) * (gx + gy)
    g[8] = (1.0 - 0.5 / tau) * w[8] * (3.0 * (( 1.0 - vx) * gx + (-1.0 - vy) * gy) )  + 9.0 * (vx - vy) * (gx - gy)
    

def collision():
    for q in range(0, Q):
        f[q] += (feq[q] - f[q]) / tau  #+ g[q]

def streaming():  
    for q in range(0, Q):
        temp = f[q].copy()
        f[q].fill(0.0)
        for y in range (1, rho.shape[0] - 1):
            for x in range (1, rho.shape[1] -1):
                f[q][y + ey[q]][x + ex[q]] = temp[y][x]

def bb_top():
    top_y = 0
    
    f2 = f[2][top_y].copy()
    f[4][top_y + ey[4]] = f2.copy()
    f[2][top_y].fill(0.0)
    
    f5 = [f[5][top_y][i] for i in range(1, x)]
    f5.append(0.0)
    f[7][top_y + ey[7]] = f5
    f[5][top_y].fill(0.0)
    
    
    f6 = [f[6][top_y][i] for i in range(0, x - 1)]
    f6 = [0.0] + f6
    f[8][top_y + ey[8]] = f6
    f[6][top_y].fill(0.0)

def bb_bot():
    bot_y = y - 1
    
    f4 = f[4][bot_y].copy()
    f[2][bot_y + ey[2]] = f4.copy()
    f[4][bot_y].fill(0.0)
    
    f7 = [f[7][bot_y][i] for i in range(0, x - 1)]
    f7 = [0.0] + f7
    f[5][bot_y + ey[5]] = f7
    f[7][bot_y].fill(0.0)
    
    
    f8 = [f[8][bot_y][i] for i in range(1, x)]
    f8.append(0.0)
    f[6][bot_y + ey[6]] = f8
    f[8][bot_y].fill(0.0)

def bb_left():
    left_x = 0
    
    f3 = f[3][:, left_x].copy()
    f[1][:, left_x + ex[1]] = f3
    f[3][:, left_x].fill(0.0)
    
    f6 = f[6][1:y-2, left_x]
    f[8][2:y-1, left_x + ex[8]] = f6
    f[6][:, left_x].fill(0.0)
    
    f7 = f[7][2:y-1, left_x]
    f[5][1:y-2, left_x + ex[5]] = f7
    f[7][:, left_x].fill(0.0)
    

def bb_right():
    right_x = x - 1
    
    f1 = f[1][:, right_x].copy()
    f[3][:, right_x + ex[3]] = f1
    f[1][:, right_x].fill(0.0)
    
    f8 = f[8][2:y-1, right_x]
    f[6][1:y-2, right_x + ex[6]] = f8
    f[8][:, right_x].fill(0.0)
    
    f5 = f[5][1:y-2, right_x]
    f[7][2:y-1, right_x + ex[7]] = f5
    f[5][:, right_x].fill(0.0)


def per_left():
    left_x = 0
    right_x = x - 2
    
    f3 = f[3][:, left_x].copy()
    
    f[3][:, right_x] = f3
    f[3][:, left_x].fill(0.0)
    
    f6 = f[6][1:y-2, left_x]
    
    f[6][1:y-2, right_x] = f6
    f[6][:, left_x].fill(0.0)
    
    f7 = f[7][2:y-1, left_x]
    
    f[7][2:y-1, right_x] = f7
    f[7][:, left_x].fill(0.0)

def per_right():
    left_x = 1
    right_x = x - 1
    
    f1 = f[1][:, right_x].copy()
    f[1][:, left_x] = f1
    f[1][:, right_x].fill(0.0)
    
    f8 = f[8][2:y-1, right_x]
    f[8][2:y-1, left_x] = f8
    f[8][:, right_x].fill(0.0)
    
    f5 = f[5][1:y-2, right_x]
    f[5][1:y-2, left_x] = f5
    f[5][:, right_x].fill(0.0)


def recalculate(rho, vx, vy):
    rho.fill(0.0)
    vx.fill(0.0)
    vy.fill(0.0)
    for q in range(0, Q):
        rho += f[q].copy()
        vx += ex[q] * f[q]
        vy += ey[q] * f[q]


def von_neumann_left(vx0, vy0):
    rho0 = 1.0
    for q in range(0, Q):
        v = ex[q] * vx0 + ey[q] * vy0
        f[q][1:y-1, 1] = w[q] * rho0 * (1.0 + 3.0 * v + 4.5 * v * v - 1.5 * (vx0 * vx0 + vy0* vy0))
    
    
#    left_x = 1
#    
#    f3 = f[3][1:y-1, left_x].copy()    # 2:y-2
#    f6 = f[6][0:y-2, left_x].copy()    # 1:y-3
#    f7 = f[7][2:y, left_x].copy()     # 3:y-1
#    
#    left_x = 1
#    
#    f0 = f[0][1:y-1, left_x].copy()    
#    f2 = f[2][0:y-2, left_x].copy()    
#    f4 = f[4][2:y, left_x].copy()
#    
#    
#    rho0 = (f0 + f2 + f4 + 2.0 * (f3 + f6 + f7)) / (1.0 - vx0)
#    print(rho0)
#    
#    f[1][1:y-1, left_x] = f3 + 2.0 / 3.0 * rho0 * vx0
#    temp = (f2 - f4) / 2.0
#    print(temp)
#    f[5][1:y-1, left_x] = f7 - temp - (vx0/6.0 - vy0/2.0) * rho0
#    f[8][1:y-1, left_x] = f6 + temp - (vx0/6.0 - vy0/2.0) * rho0
    

if __name__ == "__main__":
   
    g_calc()
    vy[1:y-1, 1].fill(0.01)
    set_up_rho(rho)
    f_eq_calc(feq,rho, vx, vy)
    f = feq
    
    print(rho)
    print(vy)
    
    for i in range(0, time):
        collision()
        streaming()
        
        von_neumann_left(0.0, 0.01)
        bb_top()
        bb_bot()
        #per_left()
        #per_right()
        #bb_left()
        bb_right()
        
        
        
        for q in range(0, Q):
            f[q][:, 0].fill(0.0)
       
#        for q in range(0, Q):
#            print('-----' + str(q) + '------')
#            print(f[q])
    
        
        recalculate(rho, vx, vy)
        f_eq_calc(feq,rho, vx, vy)
        
#        print(rho)
#        print(vy)
        
        print(str(i) + ' rho = ' + str( np.sum(rho))        )
        if((i % 100 == 0) and i > 0):
            print('Plot data')
            
            fig = plt.figure(figsize = (8,6))
            
#            x_ = np.array([j for j in range(x-5, x)])
#            y_ = np.array([j for j in range(y-5, y)])
            
            x_ = np.array([j for j in range(0, x)])
            y_ = np.array([j for j in range(0, y)])
            
            #c = plt.contourf(x_,y_,vx * vx + vy * vy)
            #b = plt.colorbar(c, orientation='vertical')
            name = 'heatmap_' + str(i) + '.png'
            
#            vx_ = np.array(vx[y-5:y, x-5:x])
#            vy_ = np.array(vy[y-5:y, x-5:x])
            
            vx_ = np.array(vx)
            vy_ = np.array(vy)
            strm = plt.streamplot(x_, y_, vx_, vy_, color= vx_ * vx_ + vy_ * vy_, linewidth=2, cmap=plt.cm.autumn)
            plt.colorbar(strm.lines)
            
            
            print(name)
            plt.savefig(name)
            plt.clf()
            
    