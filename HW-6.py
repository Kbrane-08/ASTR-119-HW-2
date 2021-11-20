#!/usr/bin/env python
# coding: utf-8

# In[ ]:


get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np


# In[ ]:


def dfdx(x,f):
    
    #d2y/dx2 = -y
    #define ：
    
    # y = f[0]
    # dy/dx = z
    # z = f[1]
    
    y = f[0]
    z = f[1]
    
    #return derivatives
    dydx = np.zeros_like(f)
    dydx[0] = z
    dydx[1] = -y
    
    return dy/dx


# In[ ]:


def cash_karp_core_mv(x_i,y_i,nv,h,f):
    
    #cash karp is defined in terms
    #of weighting variables
    ni = 7
    nj = 6
    ci = np.zeros(ni)
    alj = np.zeros( (ni.nj) )
    bi = np.zeros(ni)
    bis = np.zeros(ni)
    
    #input values for ci, alj, bi, bis
    ci[2] = 1./5.
    ci[3] = 3./10.
    ci[4] = 3./5.
    ci[5] = 1.
    ci[6] = 7./8.
    
    #j = 1
    alj[2,1] = 1./5.
    alj[3,1] = 3./40.
    alj[4,1] = 3./10.
    alj[5.1] = -11./54.
    alj[6.1] = 1631./55296.
    
    #j = 2
    alj[3,2] = 9./40.
    alj[4,2] = -9./10.
    alj[5.2] = 5./2.
    alj[6.2] = 175./512.
    
    #j = 3
    alj[4,3] = 6./5.
    alj[5,3] = -70./27.
    alj[6,3] = 575./13824.
    
    #j = 4
    alj[5,4] = 35./27.
    alj[6.4] = 44275./110529.
    
    #j = 5
    alj[6,5] = 253./4069.
    
    #bi
    bi[1] = 37./378.
    bi[2] = 0.
    bi[3] = 250./621.
    bi[4] = 125./594.
    bi[5] = 0.0
    bi[6] = 512./1771.
    
    #bis
    bis[1] = 2825./27648.
    bis[2] = 0.0
    bis[3] = 18575./48384.
    bis[4] = 13525./55296.
    bis[5] = 277./14336.
    bis[6] = 1./4.
    
    #define the k array
    ki = np.zero((ni,nv))
    
    #compute ki
    for i in range (1,ni):
        #compute xn+1 for i
        xn = x_i +ci[i]*h
        
    #compute temp y
    yn = y_i.copy()
    for j in range(1,i):
        yn += aij[i,j]*ki[j,:]
        
    #get k
    ki[i,:] = h*f(xn,yn)
    
    #get ynpo, ynpo*
    ynpo = y_i.cpoy()
    ynpos = y_i.cpoy()
    #print("ni = ",ni, ynpo, ynpos)
    for i in range(1,ni):
        ynpo += bi[i] *ki[i,:]
        ynpos += bi[i]*ki[i,:]
        #print(i,ynpo[0],ynpos[0])
         #print(i,ynpo[0],ynpos[0],bi[i]*ki[i,0],bis[i],*ki[1,0])
            
        #get error
        Delta = np,fabs(ynpo-ynpos)
        
        #print("INSIDE Delta",Delta,ki[:,0],ynpo,ynpos)
        
        #return new y and delta
        return ynpo, Delta


# In[ ]:


def cash_karp_mv_ad(dfdx,x_i,y_i,nv,h,tol):
    
    #define a safety scale
    SAFETY = 0.9
    H_NEW_FAC = 2.0
    
    #set the maximum number of iteration:
    imax = 1000
    
    #set the iteration variable:
    i = 0
    
    #create an error:
    Delta = np.full(tol,2*tol)
    
    #remember the step
    h_step = h
   
   #adjust the step
    while (Delta.max()/tol>1.0):
        
        #get our new y and error estimate
        y_ipo, Delta = cash_karp_core_mv(x_i,y_i,nv,h_step,dfdx)
    
    #if the error is too large, take a smaller step:
    if(Delta.max()/tol>1.0):
        
        #our error is too large, decrease steps
        h_step *= SAFETY * (Delta.max()/tol)**(-0.25)
        
    #check iteration
    if(i>=imax):
        print("too many iteration in cash_karp_mv_ad()")
        raise StopIteration("Ending after i = ",i)
        
      #iterate
      # i += 1
    
    #next time, try a bigger step:
    h_new = np.fmin(h_step * (Delta.max()/tol)**(-0.9), h_step * H_NEW_FAC)
    
    #return the answer and step info
    return y_ipo, h_new, h_step
    


# In[97]:


def cash_karp_mv(dfdx,a,b,y_a,tol,verbose=False):
   
   #dfdx is the derivative wrt x
   #a is the lower bound
   #b is the upper bound
   #y_aare the boundary condition at a
   #tol is the tolerance
    
   #def our starting step
   xi = 1
   yi = y_a.copy()

   #define an initial starting step
   h = 1.0e-4 * (b-a)
  
  #set the maximum number of iteration:
  #imax = 1000
    
  #set the iteration variable:
#i = 0
  
  #how many variables? 
  #nv = len(y_a)

  #set the initial condition:
  #x = np.full(1,a)
  #y = np.full( (1,nv), y_a)
    
  #set a flag
  #flag = Ture
    
  #loop until we reach b
  #while(flag):
        
        #calculate y_i+1, step info
        #y_ipo, h_new, h_step = cash_karp_mv_ad(dfdx,xi.yi,nv,h,tol)
        
        #update the step for next time
        #h = h_new
        
        #prevent an overshoot
        #if(xi+h_step>b):
            #limit step to end at b
               #h = b -xi
                
            #recompute y_i+1
            #y_ipo, h_new, h_step = cash_karp_mv_ad(dfdx,xi.yi,nv,h,tol)
            
            #we've done
            #flag = False
            
     #update the values
    # xi += h_step
    # yi = y_ipo.copy() 
    
    #add the step
    #x = np.append(x,xi)
    #y_ipo = np.zeros((len(x),nv))
    #y_ipo[0:len(x)-1,:] = y[:]
    #del y
    #y = y_ipo
    
    
    #prevent too many iterations
        #if (i>=imax):
            #print("Maximum iteration reached.")
             #raise StopIteration("Iteration number i = ",i)
                
    #iterate
    #i += 1
    
    #output some information
    #if(verbose):
    #s = "i = %3d/tx = %9.8f/ty = %9.8f/th = %9.8f/tb = %9.8" % (i,xi,yi[0],h_step,b)
    #print(s)
    
    #if we've done
    #if(xi==b):
        #flag = False
    
    #return the answer
    #return x, y


# In[104]:


a = 0.0
b = 2*np.pi
y_0 = 0.0
y_0 = 1.0
nv = 2

tolerance = 1.0e-6 
x, y = cash_karp_mv(dfdx,a,b,y_0,tolerance,verbose=Ture)
plt.plot(x,y[:,0],'o',label='y(x)')
plt.plot(x,y[:,0],'o',label='z(x)')
plt.xlim([0,2*np.pi])
plt.ylim(-1.2,1.2)
plt.legend(frameon=False)


# In[96]:


plt.plot(x, y[:,0]-np.sin(x), label="y(x) Error")
plt.plot(x, y[:,0]-np.cos(x), label="z(x) Error")
plt.legend()

