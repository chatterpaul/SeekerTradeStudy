# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 02:23:35 2021

@author: Chatterpaul
# Probabilty of Detect using Swerling-1 model:
# From Wikipedia: https://en.wikipedia.org/wiki/Fluctuation_loss
# This applies to a target that is made up of many independent scatterers of
# roughly equal areas. As few as half a dozen scattering surfaces can produce 
# this distribution. This model is particularly useful for considering aircraft 
# shapes.
"""
import numpy as np
import scipy.special as sc
import matplotlib.pyplot as plt

def pd_swerling1(nfa, nPulse, snrbar):
    # This function is used to calculate the probability of detection
    # for Swerling 1 targets.
    #
    # Inputs
    #   nfa       == Marcum's false alarm number
    #   nPulse        == number of integrated pulses
    #   snrbar    == SNR 
    #
    # outputs
    #   pd        == probability of detection
    #format long
    snrbar = 10.0**(snrbar/10.0)
    eps = 0.00000001
    delta = eps
    #% Calculate the threshold Vt
    pfa =  nPulse * np.log(2) / nfa
    sqrtpfa = np.sqrt(-np.log10(pfa))
    sqrtnp = np.sqrt(nPulse)
    vt0 = nPulse - sqrtnp + 2.3 * sqrtpfa * (sqrtpfa + sqrtnp - 1.0)
    vt = vt0
    while (delta < (vt0/10000)):
       igf = sc.gammaincc(vt0,nPulse)
       num = 0.5**(nPulse/nfa) - igf
       deno = -np.exp(-vt0) * vt0**(nPulse-1) /np.math.factorial(nPulse-1)
       vt = vt0 - (num / (deno+eps))
       delta = np.abs(vt - vt0)
       vt0 = vt

    temp1 = 1.0 + nPulse * snrbar;
    temp2 = 1.0 / (nPulse *snrbar);
    temp = 1.0 + temp2;
    val1 = temp**(nPulse-1.);
    igf1 = sc.gammaincc(vt,nPulse-1);
    igf2 = sc.gammaincc(vt/temp,nPulse-1);
    pd = 1.0 - igf1 + val1 * igf2 * np.exp(-vt/temp1);
    
    if (nPulse == 1 or pd > 1 or snrbar < .01):
       temp = -vt / (1.0 + snrbar);
       pd = np.exp(temp);
       
    return pd

#%% Main Script
pfa = 1e-6
nPulse = 1
nPulseCpi = 64
nfa = nPulse*np.log(2)/pfa
snrArray = np.arange(-30,30,1)
pdArray = np.zeros((len(snrArray),2))
for idx,snr in enumerate(snrArray):    
    snr = snrArray[idx]
    pd = pd_swerling1(nfa, nPulse, snr)
    print(f"Pd = {pd:0.3f} for pfa={pfa}(probability of false alarm), number of pulse = {nPulse} and snr={snr}(dB)")
    pdCpi = pd_swerling1(nfa, nPulseCpi, snr)
    pdArray[idx,:] = [pd,pdCpi]

fig=plt.figure(33,figsize = (9,6))
plt.plot(snrArray,pdArray)
plt.legend(('Single Pulse',f'{nPulseCpi} Pulses'))
plt.xlabel("SNR (dB))")
plt.ylabel("Pd")
plt.grid()
plt.title(f'Swerling 1 Pd Model pfa={pfa} (probability of false alarm)')