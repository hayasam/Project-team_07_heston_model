
def FourierST1(option,params):  
   N=params.NAS  
   [S0,K,sigma,T,r,divi,american]=[option.S0,option.K,option.sigma,option.T,option.r,option.divi,option.american]  
  
   j=complex(0,1)  
   #create vector in the real space  
   x_min=-7.5  
   x_max=7.5  
   dx=(x_max-x_min)/(N-1)  
   x=linspace(x_min,x_max,N)  
     
   #create vector in the fourier space  
   w_max=np.pi/dx;  
   dw=2.0*w_max/(N);  
   w=np.concatenate((linspace(0,w_max,N/2+1),linspace(-w_max+dw,-dw,N/2-1)))  
  
   # Option payoff  
   s = S0*np.exp(x);  
   v_call = np.maximum(s-K,0)  
   v_put = np.maximum(K-s,0)  
     
   # FST method  
   char_exp_factor = np.exp((j*(r-0.5*sigma**2)*w - 0.5*sigma**2*(w**2)-r)*T)  
   VC = np.real(np.fft.ifft(np.fft.fft(v_call)*char_exp_factor))  
   VP = np.real(np.fft.ifft(np.fft.fft(v_put)*char_exp_factor))  
  
   #Interpolate option prices  
   tck=si.splrep(s,VC)  
   option.call.price=si.splev(S0,tck,der=0)     
   tck=si.splrep(s,VP)  
   option.put.price=si.splev(S0,tck,der=0)    
    


