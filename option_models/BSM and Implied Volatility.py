import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss

## PARAMETERS
##
## lambda: mean-reversion speed
## vbar: long-term average volatility
## eta: volatility of vol process
## rho: correlation between stock and vol
## v0: initial volatility
## r: risk-free interest rate
## tau: time to maturity
## S0: initial share price
## K: strike price
##
## MODEL
## 
## dS_t = S_t r dt + S_t sqrt(V_t)dW_t^S
## dV_t = \lambda (\vbar - V_t)dt - eta sqrt(V_t)dW_t^V
## with d<W^S,W^V>_t = \rho dt

ONEYEAR = 250

def Moneyness (spot, strike, tau, intr): 
    K*exp(-intr*tau)/S

def BlackScholesCall (S0, strike, tau, intr, sigma, EPS=0.01):
    d1 =  (log(S0/strike) + (intr +0.5*sigma^2)*tau)/(sigma*sqrt(tau))
    d2 = d1 - sigma*sqrt(tau)
    
    if (tau < EPS):
        return(max(S0-strike,0))
    else:
        return(S0*pnorm(d1) - strike*exp(-intr*(tau))*pnorm(d2))
    

def ImpliedVolCall (S0, strike, tau, intr, price):
    f = ImpliedVolCall.BlackScholesCall(S0,strike,tau,intr,x) - price
    if (f(-1) * f(1) > 0):
        return(NA)
    '''uniroot(f,c(-1,1))$root''' #Finding the uniroot function in Python