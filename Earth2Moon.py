'''
Created on Apr 28, 2017

@author: Yunshi Zhao
'''
import random
import math
#import numpy as np
#from numpy import cross
#from scipy.integrate import ode
#import matplotlib.pyplot as plt

"""
No Uncertainties
    m_p = 1000 #kg - Mass of the projectile 
    rad_p = 1 #m - Radius of the projectile
    A_p = rad_p*2*pi() #m^2 - sectional area of the projectile
    Thrust = 4*math.pow(10,6) #N - Constant thrust for the first 30 seconds
    T_thrust = 30 #s - initial thrust time 
    rad_e = 6371*1000 #m - Radius of the earth 
    rad_m = 1737*1000 #m - Radius of the moon
    a_m = 384405*1000 #m - Semimajor axis of the moon orbit


Have Uncertainties
    m_m =  7.3476*math.pow(10,22) #kg - Mass of the moon
    T_m = 27.32*24*60*60 #s - Orbital period of the moon    
    #moon location at launch
        theta_m  #moon location at launch
        x_m  #moon location at launch
        y_m  #moon location at launch
    m_e = 5.972*math.pow(10,24) #kg - Mass of the earth   
    alpha #uniform dist from 2 to 2.5
    C_d = 0.5 #Drag coefficient, problem statement didn't say if it has uncertainty
    densities = [1.22500,0.736116,0.412707,0.193674,0.0880349,0.0394658,0.0180119,0.00821392,0.00385101,0.00188129,0.000977525]
"""

def gravitation(r1,r2,m1,m2):
    """Calculates 2 body gravitation vector
    Arguments: r1, r2 position vectors
    m1, m2 masses 
    G gravitational constant"""
    G = 6.67408 * math.pow(10,-11)
    direction = [(r2[0] - r1[0]),(r2[1] - r1[1])]
    d2 = math.pow(direction[0],2)+math.pow(direction[1],2)
    Fmag = G * m1 * m2 / d2
    Fg = [direction[0] * Fmag / math.sqrt(d2), direction[1] * Fmag / math.sqrt(d2)]
    return Fg

class moon:
    """exact"""
    a_m = 384405*1000 #m - Semimajor axis of the moon orbit
    rad_m = 1737*1000 #m - Radius of the moon
    """uncertainties"""
    T_m = 27.32*24*60*60 #s - Orbital period of the moon
    m_m =  7.3476*math.pow(10,22) #kg - Mass of the moon
    def __init__(self):
        """initial position"""
        #can replaced by different assumption
        self.theta_0 = random.random()* (0.5*math.pi) #moon location at launch
        self.x0 = self.a_m*math.cos(self.theta_0) #moon location at launch
        self.x0 = self.a_m*math.sin(self.theta_0) #moon location at launch
    def Set_m_m(self,sigma):
        #place holder, can be updated to be any distribution
        self.m_m = random.gauss(self.m_m, sigma)
    def Set_T_m(self,sigma):
        self.T_m = random.gauss(self.T_m, sigma)
    def Get_moon_position(self,t):
        """calculate moon's position"""
        self.t = t
        self.theta = self.theta_0 + 2 * math.pi * (self.t % self.T_m) * 1.0 / self.T_m
        self.x = self.a_m*math.cos(self.theta)
        self.y = self.a_m*math.sin(self.theta)

class earth():
    """exact"""
    rad_e = 6371*1000 #m - Radius of the earth
    """#uncertainties"""    
    m_e = 5.972*math.pow(10,24) #kg - Mass of the earth
    def Set_m_e(self,sigma):
        self.m_e = random.gauss(self.m_e, sigma)        
        
class projectile():
    #exact
    m_p = 1000 #kg - Mass of the projectile 
    rad_p = 1 #m - Radius of the projectile
    A_p = rad_p*2*math.pi #m^2 - sectional area of the projectile
    Thrust = 4*math.pow(10,6) #N - Constant thrust for the first 30 seconds
    T_thrust = 30 #s - initial thrust time
    #uncertainties
    alpha = 2 #uniform dist from 2 to 2.5
    C_d = 0.05 #Drag coefficient, problem statement didn't say if it has uncertainty
    densities = [1.22500,0.736116,0.412707,0.193674,0.0880349,0.0394658,0.0180119,0.00821392,0.00385101,0.00188129,0.000977525]
    rho_index = [[0]*len(densities),[0]*len(densities)]
    for i in range(len(densities)):
        rho_index[0][i] = i*5000 
        rho_index[1][i] = densities[i] 
    def __init__(self):
        self.x,self.y = 0, 6371000
        self.xd,self.yd = 0, 3074.6
        self.t = 0
    def Set_alpha(self,sigma):
        self.alpha = random.uniform(self.alpha-2.5,self.alpha+2.5)
    def Set_C_d(self,sigma):
        self.C_d = random.gauss(self.D_d, sigma)
    def Set_rho(self,sigma):
        self.rho = random.gauss(self.rho, sigma)
    def Get_rho(self):
        self.h = math.sqrt(math.pow(self.x,2)+math.pow(self.y,2)) - e.rad_e
        if self.h > max(self.rho_index[0]):
            self.rho = 0
        else:
            try:
                h_up = int(math.ceil(self.h/5)) * 5
                h_down = int(math.floor(self.h/5)) * 5
                self.rho = self.rho_index[1][self.rho_index[0].index(h_down)] + (self.h - h_down ) / 5000 * (self.rho_index[1][self.rho_index[0].index(h_up)]-self.rho_index[1][self.rho_index[0].index(h_down)])
            except:
                self.rho = self.rho_index[1][0]
    def Get_resistance(self):
            self.Get_rho()
            velsquared = math.pow(self.xd, 2) + math.pow(self.yd, 2)
            vel = math.sqrt(velsquared)
            velunitvec = [self.xd / vel, self.yd / vel]
            Fd = 0.5 * self.rho * math.pow(vel,self.alpha) * self.C_d * self.A_p
            self.Fd = [-Fd*velunitvec[0], -Fd*velunitvec[1]]
    def Get_acceleration(self,t):
        self.Get_resistance()
        m.Get_moon_position(t)
        
        self.Fge = gravitation([self.x,self.y],[0,0],self.m_p,e.m_e)
        self.Fgm = gravitation([self.x,self.y],[m.x,m.y],self.m_p,m.m_m)
        #print(self.Fge,self.Fgm,self.Fd)
        if t <= 30:
            self.xdd = (self.Fge[0] + self.Fgm[0] + self.Fd[0]) / self.m_p
            self.ydd = (self.Fge[1] + self.Fgm[1] + self.Fd[1] + p.Thrust) / self.m_p    
        else:
            self.xdd = (self.Fge[0] + self.Fgm[0] + self.Fd[0]) / self.m_p
            self.ydd = (self.Fge[1] + self.Fgm[1] + self.Fd[1]) / self.m_p
    def Get_velocity(self,td):
        self.xd = self.xd + self.xdd * td
        self.yd = self.yd + self.ydd * td
    def Get_location(self,td):
        self.x = self.x + self.xd * td
        self.y = self.y + self.yd * td   
m=moon()
e=earth()
p=projectile()
file = open('moon.csv', 'w')
for i in range(5000):
    p.Get_acceleration(i)
    p.Get_velocity(1)
    p.Get_location(1)
    string=str(p.h)+','+str(p.Fge[0])+','+str(p.Fge[1])+',' + str(p.Fgm[0])+',' + str(p.Fgm[1]) + ',' + str(p.Fd[0])+','  + str(p.Fd[1])+',' + str(p.x) + ',' + str(p.y) + ',' + str(p.xd) + ',' + str(p.yd) + ',' + str(p.xdd) + ',' + str(p.ydd)
    file.write(string+'\n')

