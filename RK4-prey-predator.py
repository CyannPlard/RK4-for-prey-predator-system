#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 2021

@author: cyannbuisson in collaboration with Valentin Poncet
"""
import numpy as np
import matplotlib.pyplot as plt

""" INTRODUCTION """

#Hi. This is a project about prey-predator system solved with the Runge-Kutta 4 method.
#You may use this code to solve another differential equations system. 
#All you have to do is :
    #changing the time parameters at the begininng of the RK4 function in the "Solver" part.
    #replacing the F function by the one you want to solve has explained in the "Function" part.
    #inserting the right initial conditions in the y0 table as explained in the "Plot" part.
    #change the plot display eventually
    
# The tyoe of problem we solve here are n differential equations of m order like :
        # dy/dx=g(y,x) with g a function determining this differential equation of 1st order
        # obviously there are then n variables to calculate (in the example above n=1)

""" RUNGE-KUTTA 4 SOLVER """



#F is the function resuming the n differential equation(s) (defined in the "function" part)
#y0 is a list containing the m*n initial conditions of the problem (explained in the "plot" part).

def RK4(F,y0):
    
    ### Time parameters of the resolution (arbitrary) ###
    x0=0
    xf=20
    h=0.1 #interval discretization, algorithm precision
    X=np.arange(x0,xf,h) #table of times/iterations used to solve the problem
    
    ### Iteration variables ###
    y=np.zeros((len(X),len(y0))) #calculation table
    y[0,:]=y0 #the initial conditions are put in the first line of the table
    
    ### Method ###
    for i in range (len(X)-1):
        k1=F(i,y[i,:])
        k2=F(i+h/2,y[i,:]+h/2*k1)
        k3=F(i+h/2,y[i,:]+h/2*k2)
        k4=F(i+h,y[i,:]+h*k3)
        y[i+1]=y[i,:]+h/6*(k1+2*k2+2*k3+k4)
    
    ### Preparing the return ###
    x = X.reshape((len(X),1))
    solution = np.append(x,y,1) #the abscissa and the solution are in the same table for plotting
    return solution

""" FUNCTION TO SOLVE """

# Entry : #time t or parameter x (here the function doesn't depend on it but it may !)
          #a list y (which is infact a line of the resolution table y in the RK4 function)
              #(here the list y has 2 values : the first y[0] is the number of prey at time t and the 2nd y[1] is the number of predator at time t)
              # y may contain the value at time t of each variable of each order involved in the equations, in the same order as in y0 (explained later)
# Return : #table containing the value of the functions determining the equation (like the value of g(x) in the example of the introduction)

def F(tx,y): #Lotka-Volterra prey-predator model 
    ### Model paramameters (arbitrary) ###
    a=2/3
    b=4/3
    c=1
    d=1
    ### Model equations (n=2, m=1) ###
    prey=y[0]*(a-b*y[1])
    predator=y[1]*(c*y[0]-d)
    return np.array([prey,predator])

""" Here an example of a second order equation (n=1, m=2)
We have the following equation:
d2y/dx2 = -y * dy/dx

We transform it in :
     dy/dx = z
     dz/dt = -yz

Then the corresponding F function in the program is :
    
def F(x,y):
    dy=y[1]
    dz=-y[0]*y[1]
    return np.array([dy,dz])
we can see that here the table y contain the variable y=[y(x),z(x)].

The initial values are :
y0[y(0), dy(0)/dx]

"""

""" PLOTING THE SOLUTIONS """

#y0 is a list containing the m*n initial conditions of the problem :
    #in order : the m conditions of the 1st equation
               #the m conditions of the 2nd equation
               # ...
               #the m conditions of the nth equation
    # the m conditions have to be put in the same order (for the orders) for each equation
y0=[0.9,0.9]
y = RK4(F, y0)

fig = plt.figure("Lotka-Volterra model")
ax1=fig.add_axes([0.1,0.1,0.8,0.35])
ax2=fig.add_axes([0.1,0.6,0.8,0.35])
ax1.plot(y[:,0],y[:,1],label="prey")
ax2.plot(y[:,0],y[:,2],label="predator")
#ax3.plot(y1[:,1],y1[:,2])