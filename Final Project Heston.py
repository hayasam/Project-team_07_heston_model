import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss
from scipy.fftpack import fft

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
    
    def HestonCallClosedForm (lmbda, vbar, eta, rho, v0, intr, tau, S0, strike):
    PIntegrand = function(u, lmbda, vbar, eta, rho, v0, intr, tau, S0, strike, j) 
    F = S0*exp(intr/tau)
    x = log(F/strike)
    a = lmbda * vbar
            
    if(j == 1):
        b = lmbda - rho* eta
        alpha =  (u^2)/2 - (u/2 * one_i) + (one_i * u)      #Changed 1i to one_i
        beta = lmbda - (rho * eta) - (rho * eta * one_i * u)
    else:
        b = lmbda
        alpha = (u^2)/2 - (u/2 * one_i)
        beta = lmbda - (rho * eta * one_i * u)
            
            
        gamma = a^2/2
        d = sqrt(beta^2 - 4*alpha*gamma)
        rplus = (beta + d)/(2*gamma)
        rminus = (beta - d)/(2*gamma)
        g = rminus / rplus
            
        D = rminus * (1 - exp(-d*tau))/(1-g*exp(-d*tau))
        C = lmbda * (rminus * tau - 2/(eta^2) * log( (1-g*exp(-d*tau))/(1-g) ) )
        
        top = p(C*vbar + D*v0 + one_i*u*x)
        bottom = (i * u)
        Re(top/bottom)


def P(lmbda, vbar, eta, rho, v0, intr, tau, S0, strike, j):
    value = integrate(PIntegrand, lower = 0, upper = Inf,
                      subdivisions=1000,
                     )
                       
    ans =  0.5 + 1/pi * value
    return (ans)
    

    A = S0*P(lmbda, vbar, eta, rho, v0, r, tau, S0, strike, 1)
    B = strike*exp(-r*tau)*P(lmbda, vbar, eta, rho, v0, intr, tau, S0, strike, 0)
    A-B

    def HestonMonteCarlo(lmbda, vbar, eta, rho, v0, intr, thau, S0, strike, nSteps=2000, nPaths=3000, vneg=2):

        n = nSteps
        N = nPaths
        
        dt = tau / n
        
        negCount = 0
        
        S = rep(S0,N)
        v = rep(v0,N)
        
        for i in range(1, n):
           
                W1 = rnorm(N);
                W2 = rnorm(N);
                W2 = rho*W1 + sqrt(1 - rho^2)*W2;

                sqvdt = sqrt(v*dt)
                S = S*exp((intr-v/2)*dt + sqrt(v * dt) * W1)
                
                if ((vneg == 3) & (2*lmbda*vbar/(eta^2) <= 1)):
                    cat("Variance not guaranteed to be positive with choice of lambda, vbar, and eta\n")
                    cat("Defaulting to Reflection + Milstein method\n")
                    vneg = 2
                

                if (vneg == 0):
                    ## Absorbing condition
                    v = v + lmbda*(vbar - v)* dt + eta * sqvdt * W2
                    negCount = negCount + length(v[v < 0])
                    v[v < 0] = 0
                
                if (vneg == 1):
                    ## Reflecting condition
                    sqvdt = sqrt(v*dt)
                    v = v + lmbda*(vbar - v)* dt + eta * sqvdt * W2
                    negCount = negCount + length(v[v < 0])
                    v = ifelse(v<0, -v, v)
                
                if (vneg == 2):
                    ## Reflecting condition + Milstein
                    v = (sqrt(v) + eta/2*sqrt(dt)*W2)^2 - lmbda*(v-vbar)*dt - eta^2/4*dt
                    negCount = negCount + length(v[v < 0])
                    v = ifelse(v<0, -v, v)     
                
                if (vneg == 3):
                    ## Alfonsi - See Gatheral p.23
                    v = v -lmbda*(v-vbar)*dt +eta*sqrt(v*dt)*W2 - eta^2/2*dt      
                
            
        
        negCount = negCount / (n*N);

        ## Evaluate mean call value for each path
        V = exp(-intr*tau)*(S>K)*(S - strike); # Boundary condition for European call
        AV = mean(V);
        AVdev = 2 * sd(V) / sqrt(N);

        list(value=AV, lower = AV-AVdev, upper = AV+AVdev, zerohits = negCount)
    


def HestonSurface(lmbda, vbar, eta, rho, v0, intr, tau, S0, strike, N=5, min_tau = 1/ONEYEAR):
    LogStrikes = seq(-0.5, 0.5, length=N)
    Ks = rep(0.0,N)
    taus = seq(min.tau, tau, length=N)
    vols = matrix(0,N,N)

    TTM = Money = Vol = rep(0,N*N)
    
    
    class HestonPrice(strike, tau):
          def HestonCallClosedForm(lmbda, vbar, eta, rho, v0, intr, tau, S0, strike):
                n = 1
                for i in range(1, N):
                    for j in range(1, N):
                
                        Ks[i] = exp(r * taus[j]+LogStrikes[i]) * S0
                        price = HestonPrice(Ks[i],taus[j])
                        iv = ImpliedVolCall(S0, Ks[i], taus[j], intr, price)
                        TTM[n] = taus[j] * ONEYEAR # in days
                        Money[n] = Moneyness(S0,Ks[i],taus[j],intr)
                        Vol[n] = iv
                        n = n+1
    

    
#the following code is for plotting the Heston Surface

def PlotHestonSurface(lmbda=6.21, vbar=0.019, eta=0.61, rho=-0.7, v0=0.010201, intr=0.0319, tau=1.0, S0=100, strike=100, N=30, min_tau = 1/ONEYEAR):
                  Ks = seq(0.8*strike, 1.25 * strike, length=N)  
                  taus = seq(0.21, tau, length=N)
        
                  HestonPrice = Vectorize(function(k, t)(HestonCallClosedForm(lmbda, vbar, eta, rho, v0, r, t, S0, k))
        
                  IVHeston = Vectorize(function(k, t)(ImpliedVolCall(S0, k, t, r, HestonPrice(k, t)))
        
        z = outer(Ks, taus, IVHeston)
                                                                     
        nrz = nrow(z)
        ncz = ncol(z)
        nb.col = 256
        color = heat.colors(nb.col)
        facet =  -(z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz])
        facetcol = cut(facet, nb.col)
        
        x=Ks
        y=taus
        z
        theta = 40
        phi = 20
        expand = 0.5
        col=color[facetcol],
        xlab="Strikes"
        ylab="Time to maturity"
        zlab="Implied Volatility"
        ticktype="detailed"
        plt.plot(x,y,z)                            
        return(invisible(z))
    