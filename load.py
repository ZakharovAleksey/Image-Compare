#! usr/bin/env python

from PIL import Image
import numpy
from numpy import linalg as LA
import math
import random 
import sys


def make_another_image(filename):
    im = Image.open(filename)
    (xDim, yDim) = im.size
    for j in range(yDim):
        for i in range(xDim):
            cur_pix_color = im.getpixel((i,j)) + random.randint(90, 100)
            im.putpixel((i,j), cur_pix_color)
    im.show()
    im.save('another.gif', 'GIF')

def show_picture(filename, _pix_vec):
    pix_vec = _pix_vec

    im = Image.open(filename)
    (xDim, yDim) = im.size

    '''max = numpy.amax(pix_vec)
    for j in range(xDim * yDim):
        pix_vec[j] = pix_vec[j] / max *255.0'''

    pix_map = numpy.zeros([yDim, xDim])

    for i in range(xDim):
        for j in range(yDim):
            pix_map[j][i] = pix_vec[j + i * yDim]
            im.putpixel((i,j), pix_map[j][i] * 10000) # 100000
    
    im.show()

class Solver:
    def __init__ (self, filename, pic_count):
        self.filename = filename
        self.pic_count = pic_count
        
        im = Image.open(self.filename)
        (self.xDim, self.yDim) = im.size
        
        self.size = self.xDim * self.yDim

        self.image = numpy.zeros([self.yDim, self.xDim])
        
        for j in range(self.yDim):
            for i in range(self.xDim):
                self.image[j][i] = im.getpixel((i,j))
        
        vec_im = []
        for i in range(self.xDim):
            for j in range(self.yDim):
                vec_im.append(self.image[j][i])
        
        self.vec_im = numpy.zeros([self.size, 1])    
        for j in range(self.size):
            self.vec_im[j] = vec_im[j]
        
        norma = LA.norm(self.vec_im)
        for j in range(self.size):
            self.vec_im[j] = self.vec_im[j] / norma

        print(filename, 'loaded successfully!')

            
    def display(self):
        print('-'*8, self.filename, '-'*8)
        print("image xDim = ", self.xDim)
        print("image yDim = ", self.yDim)

        print(self.image)

    def getSize(self):
        return self.size

    def getPicCount():
        return self.pic_count

    def getVecIm(self, j):
        return self.vec_im[j]

    def getFileName(self):
        return self.filename

    def printInFile(self):
        file_ = open("data1.txt", 'w')
        for j in range(self.yDim):
            for i in range(self.xDim):
                print(self.image[j][i], end = ' ', file = file_)
            print('\n', file = file_)
        file_.close()

class Solv:
    def __init__ (self, s1,s2,s3,s4,s5,s6, test, pic_numb, dim_form): # input s6

        self.dim_form = dim_form    # number of eigen vectors
        self.pic_numb = pic_numb    # number of load pictures
        self.size = s1.getSize()    # size of picture image
    
        self.F = numpy.zeros([self.size, self.pic_numb])

        for j in range(self.size):
            self.F[j][0] = s1.getVecIm(j)
            self.F[j][1] = s2.getVecIm(j)
            self.F[j][2] = s3.getVecIm(j)
            self.F[j][3] = s4.getVecIm(j)
            self.F[j][4] = s5.getVecIm(j)
            self.F[j][5] = s6.getVecIm(j)
        
        self.Tau = numpy.zeros([self.pic_numb, self.pic_numb])
        
        self.eigen_glos = {}
        self.tau_glos = {}
        
    def Solve_algorithm(self, test):
        self.Tau = numpy.dot(self.F.T, self.F)
        
        eigen_val, eigen_vec = numpy.linalg.eig(self.Tau)#numpy.LA.eig(self.Tau)
        print(eigen_val)
        print(eigen_vec)
        
        for i in range(self.dim_form):
            print(eigen_val[i]* eigen_vec[i])
            print('--'*8)
            print(numpy.dot(self.Tau, eigen_vec[i]))
            print('--'*8)
            
        print('\n'*3, 'in save values:')
        for i in range(self.dim_form):
            temp = numpy.zeros([self.pic_numb, 1])
            for j in range(self.pic_numb):
                temp[j] = eigen_vec[numpy.argmax(eigen_val) + i][j]
            
            self.eigen_glos[numpy.amax(eigen_val)] = temp
            eigen_val = numpy.delete(eigen_val, numpy.argmax(eigen_val))
        
        for key in self.eigen_glos:
            self.tau_glos[key] = numpy.dot(self.F, self.eigen_glos[key])
            self.tau_glos[key]/= math.sqrt(key)
            self.tau_glos[key]/= LA.norm(self.tau_glos[key])
            
            

        '''for key in self.tau_glos:
            temp = numpy.zeros([self.size, 1])
            for j in range(self.size):
                temp[j] = numpy.absolute(self.tau_glos[key][j])
            show_picture("for_draw.gif", temp)'''

        g = numpy.zeros([self.size, 1])
        for j in range(self.size):
            g[j] = test.getVecIm(j)
        
        Pg = numpy.zeros([self.size, 1])
        
        for key in self.tau_glos:
            a = numpy.dot(g.T, self.tau_glos[key])
            Pg += self.tau_glos[key]*a;

  
        Pg/= LA.norm(Pg)
        
        norma = LA.norm(g - Pg)
        print(norma)
        if norma < 0.21:
            print("Image", test.getFileName()," has the same form!")
        else:
            print("Image", test.getFileName()," does NOT the match the form!")

#----------------------------------------------------------- 25 звонить профилак
if __name__ == '__main__':
    print('-' * 8, "TEST", '-' * 8)
    
    s1 = Solver("s1.gif",1)
    s2 = Solver("s2.gif",2)
    s3 = Solver("s3.gif",3)
    s4 = Solver("s4.gif",4)
    s5 = Solver("s5.gif",5)
    s6 = Solver("s6.gif",6)
    test = Solver("test.gif", 4)
    wrong = Solver("wrong.gif", 4)
    

    solv = Solv(s1,s2,s3,s4,s5,s6,test, 6, 2) # s6
    solv.Solve_algorithm(wrong)
    #solv.Solve_algorithm(test)
