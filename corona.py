#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 13:43:26 2018

@author:  (based on a code from Oregon State)


S0 - susceptable
I0 - infected 
R0 - recovered

mu - the probability of somone getting susceptable with a vaccine
    36% effective from time
        http://time.com/5415276/flu-shot-2018/
        https://www.cdc.gov/mmwr/volumes/67/wr/mm6706a2.htm?s_cid=mm6706a2_w
alpha - the probability of being infected with teh corona viris

beta - the probability of someone recovering teh corona virus

delta - the death rate of the corona virus
https://www.worldometers.info/coronavirus/coronavirus-death-rate/

"""
# This loads some pacakges that have arrays etc and the ODE integrators
import scipy, scipy.integrate
import pylab
    

# Parameters
alpha = .89 #######
beta = .3
delta = .1


# Initial condition
S0 = 58
I0 = 5
R0 = 0
#dead
#D0 = 0

Y0 = [ S0, I0, R0 ]

tMax = 1

# Time vector for solution
T = scipy.linspace(0, tMax, 10)


# This defines a function that is the right-hand side of the ODEs
# Warning!  Whitespace at the begining of a line is significant!
def rhs(Y, t, beta, gamma, mu):
    '''
    SIR model.
    
    This function gives the right-hand sides of the ODEs.
    '''
    
    # Convert vector to meaningful component vectors
    # Note: Indices start with index 0, not 1!
    S = Y[0]
    I = Y[1]
    R = Y[2]
    #D = Y[3]
    
    N = (5832710) 
    #number of people getting the flu shot
    
    # The right-hand sides
    dS =  alpha * I * S  - beta * I - delta * I 
    dI = -alpha * I * S 
    dR = beta * I
    
    
    # Convert meaningful component vectors into a single vector
    dY = [ dS, dI, dR ]

    return dY

# Integrate the ODE
# Warning!  The ODE solver over-writes the initial value.
# This can be a pain in the ass if you run the ODE multiple times.
# Also, 'args' passes parameters to right-hand-side function.
solution = scipy.integrate.odeint(rhs,
                                  Y0,
                                  T,
                                  args = (alpha, beta, delta))
        
S = solution[:, 0]
I = solution[:, 1]
R = solution[:, 2]

#N = S + I + R
#D = delta * I


# Make plots

# Load a plotting package
# PyLab is motivated by Matlab...


# I usually use PyLab for quick plots
# and the Python GnuPlot package for publication

pylab.figure()

pylab.plot(T, S,
           T, I ,
           T, R,)

pylab.xlabel('Time')
pylab.ylabel('Population')

pylab.legend([ 'Susceptible', 'Infective', 'Recovered'])

# Actually display the plot
pylab.show()