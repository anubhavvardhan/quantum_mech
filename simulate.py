#!/usr/bin/python
from numpy import *
from scipy import integrate
import matplotlib.pyplot as plt

class laserProfile:
	"""This is the class with the details of the laser."""

	def __init__(self, param, omega, t, dt, gaussian):
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
		self.amplitude = pi/trapz(self.amplitude,dx=self.dt) * self.amplitude * cos(omega*(self.time-self.totalTime))

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
		plt.show()

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
		plt.show()

if __name__=="__main__":
	laser=laserProfile(0.009,10,50,0.01,True)
	laser.plot()
	p = propagate(1,1,laser,10)
	p.timePropagate(1,0)
	p.plot()
