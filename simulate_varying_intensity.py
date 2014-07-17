#!/usr/bin/pythonfrom numpy 
import *from scipy 
import integrate
import matplotlib.pyplot as plt
import os

class laserProfile:	
"""This is the class with the details of the laser."""	

def __init__(self, param, omega, area, t, dt, gaussian):		
self.param = param		
self.omega = omega		
self.totalTime = t		
self.dt = dt		
self.isGaussian = gaussian		
self.time = linspace(0,2*self.totalTime,2*self.totalTime/self.dt)		
if gaussian:			
self.amplitude = exp(-param*(self.time-self.totalTime)**2)		
else:			
self.amplitude = 1/pi * param/(param**2+(self.time-self.totalTime)**2)		
self.amplitude = area/trapz(self.amplitude,dx=self.dt) * self.amplitude * cos(omega*(self.time-self.totalTime))	

def plot(self):		
plt.plot(self.time,self.amplitude)		
plt.ylabel("Amplitude")		
plt.xlabel("Time")		
plt.title("Laser Profile")		
if self.isGaussian:			
file_name = "Gaussian"		
else:			
file_name = "Lorentzian"		
file_name+="_"+str(self.param)+"_"+str(self.omega)+".png"		
plt.savefig(file_name)		
plt.show()def compare(laser1,laser2):	
diff = 0;	
for i in range(len(laser1.time)):
diff+=(laser1.amplitude[i]-laser2.amplitude[i])**2	
return diff


class propagate:	
def __init__(self,d,hbar,laser,omega0):		
self.hbar = hbar		
self.d = d		
self.omega0 = omega0		
self.laser = laser		
self.cg=[]		
self.ce=[]	
def timePropagate(self,cg,ce):
cg0 = cg		
ce0 = ce		
for time in range(len(self.laser.time)):			
self.cg.append(abs(cg0)**2/(abs(cg0)**2+abs(ce0)**2))			
self.ce.append(abs(ce0)**2/(abs(cg0)**2+abs(ce0)**2))			
oldcg0 = cg0			
cg0 = cg0 + 1j * self.laser.dt/self.hbar * self.d * self.laser.amplitude[time] * ce0 * exp(-1j * self.omega0 * (self.laser.time[time]-self.laser.totalTime))			
ce0 = ce0 + 1j * self.laser.dt/self.hbar * self.d * self.laser.amplitude[time] * oldcg0 * exp(1j * self.omega0 * (self.laser.time[time]-self.laser.totalTime))	

def plot(self):		
plt.plot(self.laser.time,self.cg,label="Ground State")		
plt.plot(self.laser.time,self.ce,label="Excited State")		
plt.legend(loc='upper left')		
plt.xlabel("Time")		
plt.ylabel("Population")		
plt.title("Population Progression")		
if self.laser.isGaussian:			
file_name = "Gaussian"		
else:			
file_name = "Lorentzian"		
file_name+="_"+str(self.laser.param)+"_"+str(self.laser.omega)+";"+str(self.omega0)+".png"		
plt.savefig(file_name)		
plt.show()if __name__=="__main__":	
omega = 0.5	
Time = 50	
dt = 0.01	
omega0 = omega	
d = 1	
hbar = 1	
param = 0.009	
cg = 1	
ce = 0	
cg_list = []	
ce_list = []	
area = []	
for i in arange(0,500,10):		
gaussian=laserProfile(param, omega, 2*pi*i/500, Time, dt, True)		
area.append(2*pi*i/500)		
if i==250:			
gaussian.plot()	#gaussian = laserProfile(param,omega,1,Time,dt,True)		
p = propagate(d,hbar,gaussian,omega0)		
p.timePropagate(cg,ce)		
cg_list.append(p.cg[-1])		
#print p.cg[-1]		
ce_list.append(p.ce[-1])	
plt.plot(arange(1,500,10),cg_list,label="Ground State")	
plt.plot(arange(1,500,10),ce_list,label="Excited State")	
plt.ylabel("Population")	
plt.xlabel("Area")	
fil = str(param)+"_"+str(omega)+".png"	
plt.savefig(fil)	
plt.show()#	
gaussian.plot()#	
p.plot()#	
p1 = propagate(d,hbar,lorentzian,omega0)#	
p1.timePropagate(cg,ce)#	
p1.plot()