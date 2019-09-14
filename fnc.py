from numpy import *
from scipy import *



#**************************************************************
#**************************************************************
#This module define the function used for the DIRECT optimization


#Some test function

#Branin function variable range
min_arr=array([-5.0,0.0])
max_arr=array([10.0,15.0])

#Holder function variable range
#min_arr=array([-10.0,-10.0])
#max_arr=array([10.0,10.0])

#Camel function variable range
#min_arr=array([-3,-2])
#max_arr=array([3,2])


delta_arr=max_arr-min_arr

#get the absolute coordinates of the function minimum
def get_opt(val):
    return min_arr+val*delta_arr


#Defining the function
def real_fnc(arr): 
    
    #Branin fnc
    return 1.0*((arr[1]-(5.1/(4*(pi**2.0)))*(arr[0]**2.0)+\
    (5.0/pi)*arr[0]-6.0)**2.0)+10.0*(1.0-(1.0/(8.0*pi)))*cos(arr[0])+10.0 

    #Holder table fnc
    #return -abs(sin(arr[0])*cos(arr[1])*exp(abs(1.0-sqrt(dot(arr,arr))/pi)))

    #6-hump camel fnc
    #return((4.0-2.1*(arr[0]**2)+(arr[0]**4)/3.0)*(arr[0]**2)+arr[0]*arr[1]+(-4.0+4.0*(arr[1]**2))*(arr[1]**2))
    
    
#real_fnc using normalized coordinates as input.
def target_fnc(arr): 
    coord=min_arr+delta_arr*arr
    return real_fnc(coord)
    
