# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 03:49:44 2021

@author: Chatterpaul
The Famous Range Range Equation
This show the basics of Rectangular Phase Array Antenna Patterns
# Reference: https://www.routledge.com/Radar-Systems-Analysis-and-Design-Using-MATLAB/Mahafza/p/book/9781439884959
# https://www.routledge.com/downloads/K13984/K13984_Downloads.zip :: radar_eq.m
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Pd_Swerling1 import pd_swerling1

def radar_eq(pt, freq, g, sigma, b, nf, loss, rng):
    #%% Monostatic Case
    # This function implements Eq. (2.22) of textbook
    #
    # Inputs:
    #   pt        == input peak power in Watts
    #   freq      == radar operating frequency in Hz
    #   g         == antenna gain in dB
    #   sigma     == radar cross section in meter squared
    #   b         == radar bandwidth in Hz
    #   nf        == noise Figure in dB
    #   loss      == total radar losses in dB
    #   rng     == range to target (single value or vector) in Km
    #    
    # Outputs:
    #   snr       == SNR in dB     
    # 
    c = 3.0e+8 # speed of light
    lambdaM = c / freq # wavelength
    p_peak = 10*np.log10(pt) # convert peak power to dB
    lambda_sqdb = 10*np.log10(lambdaM**2) # compute wavelength square in dB
    sigmadb = 10*np.log10(sigma) # convert sigma to dB
    four_pi_cub = 10*np.log10((4.0 * np.pi)**3) # (4pi)^3 in dB
    k_db = 10*np.log10(1.38e-23) # Boltzman's constant in dB
    to_db = 10*np.log10(290) # noise temp. in dB
    b_db = 10*np.log10(b) # bandwidth in dB
    range_pwr4_db = 40*np.log10(rng) # vector of target range^4 in dB
    # Implement Equation (2.22)
    num = p_peak + g + g + lambda_sqdb + sigmadb
    den = four_pi_cub + k_db + to_db + b_db + nf + loss + range_pwr4_db
    snr = num - den
    return snr

def radar_eq_bi(pt, freq, gtx, grx, sigma, b, nf, loss, rngtx, rngrx):
    # Bi-static or Semi-Active (non co-located Tx and Rx)
    # This function implements Eq. (2.22) of textbook
    #
    # Inputs:
    #   pt        == input peak power in Watts
    #   freq      == radar operating frequency in Hz
    #   gtx       == tx antenna gain in dB
    #   grx       == rx antenna gain in dB
    #   sigma     == bi-static radar cross section in meter squared
    #   b         == radar bandwidth in Hz
    #   nf        == noise Figure in dB
    #   loss      == total radar losses in dB
    #   rngtx     == tx antenna range to target (single value or vector) in Km
    #   rngrx     == rx antenna range to target (single value or vector) in Km
    #    
    # Outputs:
    #   snr       == SNR in dB     
    # 
    c = 3.0e+8 # speed of light
    lambdaM = c / freq # wavelength
    p_peak = 10*np.log10(pt) # convert peak power to dB
    lambda_sqdb = 20*np.log10(lambdaM) # compute wavelength square in dB
    sigmadb = 10*np.log10(sigma) # convert sigma to dB
    four_pi_cub = 30*np.log10(4.0 * np.pi) # (4pi)^3 in dB
    k_db = 10*np.log10(1.38e-23) # Boltzman's constant in dB
    to_db = 10*np.log10(290) # noise temp. in dB
    b_db = 10*np.log10(b) # bandwidth in dB
    tx_range_pwr2_db = 20*np.log10(rngtx) # vector of target range^2 in dB
    rx_range_pwr2_db = 20*np.log10(rngrx) # vector of target range^2 in dB
    # Implement Equation (2.22)
    num = p_peak + gtx + grx + lambda_sqdb + sigmadb
    den = four_pi_cub + k_db + to_db + b_db + nf + loss + \
        tx_range_pwr2_db + rx_range_pwr2_db
    snr = num - den
    return snr
#%% Main Script

'''
rng = np.linspace(100,20000,100)
snr_ms = radar_eq(pt=10, freq=30e9, g=24, sigma=.01, b=1e6, nf=6, loss=1, 
               rng=rng)


snr_bs = radar_eq_bi(pt=100, freq=30e9, gtx=40, grx=24, sigma=.01, b=1e6, nf=6, loss=1,
               rngtx = np.linspace(1e4,2e4,100), rngrx = rng)
'''

km2m = 1e3;
txdf = pd.read_csv('D:\work_project\STK\Interceptor_UAS\Facility-Radar-To-Aircraft-Target1_ScanEagle_AERAtt.csv')
txTime = txdf['Time (EpSec)']
txRng = txdf['Range (km)']*km2m
rxdf = pd.read_csv('D:\work_project\STK\Interceptor_UAS\Aircraft-Interceptor1-To-Aircraft-Target1_ScanEagle_AERAtt.csv')
#print(rxdf)
#rxdf.columns
rxTime = rxdf['Time (EpSec)']
rxRng = rxdf['Range (km)']*km2m


fig=plt.figure(1,figsize = (9,6))
plt.plot(txTime,txRng)
plt.plot(rxTime,rxRng)
plt.xlabel("Time (sec))")
plt.ylabel('Tx/Rx Range (km)')
plt.grid()
#plt.title(f'Swerling 1 Pd Model pfa={pfa} (probability of false alarm)')

rcs_dBsm = -20
rcs = 10**(rcs_dBsm/10)
tIdx = min(len(rxTime),len(txTime))
bw = 100e3
snr_ms = radar_eq(pt=10, freq=30e9, g=24, sigma=rcs, b=bw, nf=6, loss=1, 
               rng=rxRng)
snr_bs = radar_eq_bi(pt=100, freq=30e9, gtx=40, grx=24, sigma=rcs, b=bw, nf=6, loss=1,
               rngtx = txRng[:tIdx], rngrx = rxRng[:tIdx])

rng = rxdf['Range (km)'][:tIdx]


fig=plt.figure(2,figsize = (9,6))
plt.plot(rng,snr_bs)
plt.plot(rng,snr_ms)
plt.legend(('Mono-Static','Bi-Static'))
plt.xlabel('Rx Range (km)')
plt.ylabel('SNR (dB)')
plt.gca().invert_xaxis()
plt.grid()
#plt.title(f'Swerling 1 Pd Model pfa={pfa} (probability of false alarm)')


pfa = 1e-6
nPulse = 1
nPulseCpi = 64
nfa = nPulse*np.log(2)/pfa
#snrArray = snr_ms
snrArray = snr_bs
pdArray = np.zeros((len(snrArray),2))
for idx,snr in enumerate(snrArray):    
    snr = snrArray[idx]
    pd = pd_swerling1(nfa, nPulse, snr)
    pdCpi = pd_swerling1(nfa, nPulseCpi, snr)
    pdArray[idx,:] = [pd,pdCpi]

#fig=plt.figure(33,figsize = (9,6))
fig, ax = plt.subplots()
ax2 = ax.twinx()
ax.plot(rng,pdArray[:,0],'-.',color='blue')
ax.plot(rng,pdArray[:,1],'-',color='blue')
ax2.plot(rng,snr_bs,color='orange')
ax.legend(('Single Pulse',f'{nPulseCpi} Pulses'))
ax.set_xlabel("Rx Range (km)")
ax.set_ylabel("Pd",color='blue')
ax.set_ylim(0,1)
ax2.set_ylabel("SNR",color='orange')
ax.invert_xaxis()
ax.grid()
#plt.title(f'Swerling 1 Pd Model pfa={pfa} (probability of false alarm)')