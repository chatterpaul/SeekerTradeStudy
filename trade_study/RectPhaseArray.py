# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 21:37:28 2021

@author: Chatterpaul
This show the basics of Rectangular Phase Array Antenna Patterns
# Reference: https://www.routledge.com/Radar-Systems-Analysis-and-Design-Using-MATLAB/Mahafza/p/book/9781439884959
# https://www.routledge.com/downloads/K13984/K13984_Downloads.zip :: rect_array.m
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def rect_array(Nxr,Nyr,dolxr,dolyr,theta0,phi0,winid,win,nbits):
    #%%%%%%%%%%%%%%%%%%%% ************************ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # This function computes the 3-D directive gain patterns for a plannar array
    # This function uses the fft2 to compute its output
    #%%%%%%%%%%%%%%%%%%%% ************  INPUTS ************ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Nxr ==> number of along x-axis; Nyr ==> number of elements along y-axis
    # dolxr ==> element spacing in x-direction; dolyr ==> element spacing in y-direction Both are in lambda units
    # theta0 ==> elevation steering angle in degrees, phi0 ==> azimuth steering angle in degrees
    # winid ==> window identifier; winid negative ==> no window ; winid positive ==> use window given by win
    # win ==> input window function (2-D window) MUST be of size (Nxr X Nyr)
    # nbits is the number of nbits used in phase quantization; nbits negative ==> NO quantization
    #%%%%%%%%%%%%%%%%%%%% *********** OUTPUTS ************* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # pattern ==> directive gain pattern
    #%%%%%%%%%%%%%%%%%%%% ************************ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    eps = 0.0001
    nx = np.arange(0,Nxr)
    ny = np.arange(0,Nyr)
    i = 1j

    # check that window size is the same as the array size
    nw = np.size(win,0);
    mw = np.size(win,1);
    if winid > 0:
        if nw != Nxr:
            print('STOP == Window size must be the same as the array')
            return
        if mw != Nyr:
            print('STOP == Window size must be the same as the array')
            return

    #if dol is > 0.5 then; choose dol = 0.5 and compute new N
    if (dolxr <=0.5):
       ratiox = 1
       dolx = dolxr
       Nx = Nxr
    else:
       ratiox = np.ceil(dolxr/.5)
       Nx = (Nxr -1 ) * ratiox + 1
       dolx = 0.5
    
    if (dolyr <=0.5):
       ratioy = 1
       doly = dolyr
       Ny = Nyr
    else:
       ratioy = np.ceil(dolyr/.5)
       Ny = (Nyr -1) * ratioy + 1
       doly = 0.5

    # choose proper size fft, for minimum value choose 256X256
    Nrx = 10 * Nx; 
    Nry = 10 * Ny;
    nfftx = 2**(np.ceil(np.log(Nrx)/np.log(2)));
    nffty = 2**(np.ceil(np.log(Nry)/np.log(2)));
    if nfftx < 256:
        nfftx = 256
    if nffty < 256:
        nffty = 256;

    # generate array of elements with or without window
    if winid < 0: 
        array = np.ones((Nxr,Nyr))
    else:
        array = win

    # convert steering angles (theta0, phi0) to radians
    theta0 = theta0 * np.pi / 180;
    phi0 = phi0 * np.pi / 180;
    # convert steering angles (theta0, phi0) to U-V sine-space
    u0 = np.sin(theta0) * np.cos(phi0);
    v0 = np.sin(theta0) * np.sin(phi0);
    print(f'Steering to UV={u0:0.3f},{v0:0.3f}')
    
    # Use formula thetal = (2*pi*n*dol) * sin(theta0) divided into 2^m levels
    # and rounded to the nearest qunatization level
    if nbits < 0:
        phasem = np.exp(i*2*np.pi*dolx*u0 * nx * ratiox)
        phasen = np.exp(i*2*np.pi*doly*v0 * ny * ratioy)
    else:
        levels = 2**nbits
        qlevels = 2.0*np.pi / levels # compute quantization levels
        sinthetaq = np.round(dolx * nx * u0 * levels * ratiox) * qlevels # vector of possible angles
        sinphiq = np.round(doly * ny * v0 * levels * ratioy) * qlevels # vector of possible angles
        phasem = np.exp(i*sinthetaq)
        phasen = np.exp(i*sinphiq)     
  
    # add the phase shift terms
    array = array * (np.outer(phasem,phasen));
  

    # determine if interpolation is needed (i.e N > Nr)
    w = np.zeros((int(Nxr*ratiox),int(Nyr*ratioy)))*i
    if (Nx > Nxr) or (Ny > Nyr):
        for xloop in np.arange(Nxr): #1 : Nxr
            temprow = array[xloop,:]
            w[int(xloop*ratiox),:int(Ny):int(ratioy)] = temprow
        array = w;
    else:
        w = array

    w = array;
    # Compute array pattern
    wpadd = np.pad(w,[(0,int(nfftx-np.size(w,0))),(0,int(nffty-np.size(w,1)))]);
    arrayfft = np.abs(np.fft.fftshift(np.fft.fft2(wpadd)))**2
    
    #compute [su,sv] matrix
    U = np.arange(-nfftx/2,nfftx/2)/(dolx*nfftx)
    #indexx = np.argwhere(np.abs(U) <= 1)
    #U = U[indexx]
    V = np.arange(-nffty/2,nffty/2)/(doly*nffty)
    #indexy = np.argwhere(np.abs(V) <= 1)
    #V = V[indexy]
    
    # Normalize to generate gain patern
    #rbar = np.sum(arrayfft[indexx,indexy]) / dolx/doly/4./nfftx/nffty;
    #arrayfft = arrayfft(indexx,indexy) /rbar;
    rbar = np.sum(arrayfft) / dolx/doly/4./nfftx/nffty
    arrayfft = arrayfft/rbar

    SV,SU = np.meshgrid(V,U);
    indx = (SU**2 + SV**2) > 1
    arrayfft[indx] = eps/10;

    pattern = 10*np.log10(arrayfft + eps)

    plt.figure(1)#,figsize = (9,6))
    #plt.contour(SV,SU,pattern)
    plt.imshow(pattern,extent=[-1,1,-1,1],origin = 'lower' )
    plt.plot(v0,u0,'m*')
    plt.colorbar()
    plt.xlabel('v')
    plt.ylabel('u')
    plt.title(f'Steered:u0={u0:0.3f}, v0={v0:0.3f}')

    fig = plt.figure(2,figsize = (9,6))
    ax = fig.gca(projection='3d')
    dist = np.round(pattern)
    dist_max = np.max(dist)
    my_col = cm.jet(dist/dist_max)
    ax.plot_surface(SV, SU, pattern, rstride=1, cstride=1, facecolors=my_col,
        linewidth=0, antialiased=False)
    plt.xlabel('v')
    plt.ylabel('u')
    ax.set_zlabel('gain (dB)')
    
    return pattern
#%% Main Script
#x_min = 16
dolxr = 0.8
dolyr = 0.8
#dolxr = 0.5
#dolyr = 0.5
#theta0 = 30
#phi0 = 30
theta0 = 0#40
phi0 = 0#-15
winid = 1
nbits = -1

Nxr = 16 #x_min
Nyr = 16 #1*x_min
if winid > 0:
    win = np.ones((Nxr,Nyr))
else:
    win = np.ones((1,1))

antPat = rect_array(Nxr,Nyr,dolxr,dolyr,theta0,phi0,winid,win,nbits)
antPatMax = np.max(antPat)
arrayGain = 10*np.log10(Nxr*Nyr)
bwx = 1.2/(dolxr*Nxr)
bwy = 1.2/(dolyr*Nyr)
bw = np.min([bwx,bwy])
antPatMaxCalc = 10*np.log10(4*np.pi/(bwx*bwy))
# need a factor of 2 before of both tx & rx?
#antPatMaxCalc = 10*np.log10(4*np.pi/(bwx*bwy*2)) 

print(f'Antenna Pattern Gain = {antPatMax:0.2f}(dB)')
print(f'Array Gain = {arrayGain:0.2f}(dB)')
print(f'Antenna Calculated Gain = {antPatMaxCalc:0.2f}(dB)')
print(f'Antenna Calculated 3dB-BW Az,El = {np.rad2deg(bwx):0.2f}, {np.rad2deg(bwy):0.2f}(deg)')

