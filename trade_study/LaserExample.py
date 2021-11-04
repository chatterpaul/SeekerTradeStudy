# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 11:18:18 2021

@author: Chatterpaul
# This show basic Laser Link Power Budget calculation
"""

import numpy as np
import matplotlib.pyplot as plt


ang = np.linspace(-np.pi, np.pi,256)
focalSpot = np.abs(np.sinc(ang))

plt.figure(1,figsize = (9,6))
plt.plot(focalSpot,color='k')
plt.axis('off')

# laser performance example 
P = 1e3 # power (watt)
tx = .8 # transmission fraction tx through the air/hae/clouds/etc
F = 20e3 # range to target (m)
D = .1 # beam aperture size (m)
L = 1e-6 # wavelength (m)
BQ = 2 # beam quality
DL = L/D # diffraction limit divergence angle (radians - small angle approximation)
d = F*DL # ideal spot size (m)
rd = BQ*d  # realistic beam quality focal spot
convCm2 = 100*100 # convert the area from meter-squared to centimeter-squared 
I = (P*tx)/rd**2 *1/convCm2 # intensity (W/cm^2)

print(f'This {P*1e-3:.1f} killowatt laser has a divergence angle'
      f' {DL*1e6:.1f} microradians,\n focal spot size {rd*1e2:.1f} centimeters'
      f' and intensity {I:.1f} W/cm^2.')

Frng = np.linspace(100,20e3,100) # range array (m)
d = Frng*DL # ideal spot size (m)
rd = BQ*d  # realistic beam quality focal spot
I = (P*tx)/rd**2 *1/convCm2 # intensity (W/cm^2)

plt.figure(2,figsize = (9,6))
plt.plot(Frng,I)
plt.xlabel('Range F (km)')
plt.ylabel('Intensity (W/cm$^{2}$)')
plt.title(f'Laser Example:\nP={P*1e-3}(kW), tx={tx}, D={D*100}(cm),'
          f' $\lambda$={L*1e6:.1f}(um), BQ={BQ}')
plt.yscale('log')
plt.grid()

#lidar/ladar range equation

PW = 10e-9 # 10e-12 # pulse-width (sec)
# Note Energy(Joules) = Power(Watts)*Time(Seconds)
# *** Assumption: Using previous values EP would only be 10 uW, is that reasonable?
EP = .01 #P*PW # pulse-energy (J),
K = 1/np.exp(1) # illumination constant
BQ = 2 # beam quality
L = 1e-6 # wavelength (m)
D = 1#.1 # beam aperature size (set for both tx and rx)
R = np.linspace(100,20e3,100)
Nsys = .7 # transmittance efficiency (set for both tx and rx)
rho = .1 # reflectivity of target
# *** Assumption: is using  a 40 cm focal spot?
ocs = 4 * rho * .4**2 # optical cross section (m^2)
alp = 0.1 # atmospheric extinction coefficient (loss/km)
atmL =  np.exp(-2*alp*R*1e-3) # atmospheric loss vs range
Pr = EP*K*ocs/(16*np.pi*BQ**2*L**2) * (D**4/R**4) * (Nsys**2) * atmL

plt.figure(3,figsize = (9,6))
#plt.plot(R,Pr*1e9)
#plt.ylabel('Power Received (nW)')
plt.plot(R,10*np.log10(Pr)+30)
plt.ylabel('Power Received (dBm)')
plt.xlabel('Range (km)')
plt.title(f'LIDAR Power Detected Example:\nEP={EP*1e3}(mJ), PW={PW*1e9}(ns), '
          f'D={D*100}(cm), $\lambda$={L*1e6:.1f}(um), BQ={BQ}, OCS={ocs*convCm2:.1f}(cm$^{2}$)')
#plt.yscale('log')
plt.grid()

plt.figure(4,figsize = (9,6))
plt.plot(R,atmL)
plt.xlabel('Range (km)')
plt.ylabel('Atmospheric Loss')
plt.title(f'Atmospheric Loss with $\\alpha$={alp}')
plt.yscale('log')
plt.grid()