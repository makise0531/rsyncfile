# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

def doppler(E1,E2,kB,T,M,mn):
    delta=sqrt((4.0*E1*kb*T)/(M/mn))
    S=1.0/(delta*np.qurt(np.pi))*np.exp(-(E2-E1)/delta)

    return S

def BW(E2,Es,ks,gs,gamma_ns,gamma_rs,Ep,kp,gp,gamma_np,gamma_rp,):
    gamma_s=gamma_ns+gamma_rs
    gamma_p=gamma_np+gamma_rp
    a_0s=1.0/(4.0*ks**2)*np.sqrt(Es/E2)*(gs*gamma_ns*gamma_rs)/((E2-Es)**2+gamma_s**2/4)
    a_0p=1.0/(4.0*kp**2)*np.sqrt(E2/Ep)*(gp*gamma_np*gamma_rp)/((E2-Ep)**2+gamma_p**2/4)
    return a_0s + a_0p

def gs_f(Is,Js):
    return 2.0*Js+1.0/2.0*(2.0*Is+1.0)

def gp_f(Ip,Jp):
    return 2.0*Jp+1.0/2.0*(2.0*Ip+1.0)

Is=1.0/2.0
Js=1.0
Es=27.5
ks=
gs=gs_f(I,J)

X = np.linspace(-np.pi, np.pi, 256, endpoint=True)
C, S = np.cos(X), np.sin(X)

plt.plot(X, C)
plt.plot(X, S)

plt.show()
