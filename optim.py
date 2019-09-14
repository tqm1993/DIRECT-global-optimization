from numpy import *
from scipy import *
from scipy import optimize
from fnc import *
import matplotlib.pyplot as plt
import time


#*****************************************************************
#*****************************************************************
#This module contain DIRECT global optimization code. The theory of the method can be found
#in Jones D.R, Perttunen C.D, Stuckman B.E, "Lipschitzian optimization without the Lipschitz constant"
#Journal of Optimization Theory and Application, 79(1), p.157 (1993)
#*****************************************************************
#*****************************************************************


#Assuming that we have a multi-dimensional functional, defined in a constrained vectorial space
#The DIRECT method relies on dividing such vectorial space into small multi-dimensional "rectangle".
#It compares the function value of each rectangle to converge towards the minima.

#This class rect_zone store information of all created rectangles, and also defines some operator acted on them.
class rect_zone():

    def __init__(self,dim,norm,center_coord,center_value): #in normalized coord
        self.center_coord=center_coord
        self.center_value=center_value
        self.dim=dim
        self.norm=norm
        self.dist=sqrt(dot(self.norm,self.norm))/2.0
        self.max_side,self.max_indices=self.all_max(self.norm)
        
    def all_max(self,arr):
        max_indices=array([],dtype=int)
        max_val=max(arr)
        for i in range (0,len(arr)):
            if (arr[i]==max_val):
                max_indices=append(max_indices,i)
            
        return (max_val,max_indices)

    def update_rect(self):
        self.dist=sqrt(dot(self.norm,self.norm))/2
        self.max_side,self.max_indices=self.all_max(self.norm)
  
    def divide_rect(self):
        divide_arr=ones(self.dim)
        def check_dum(i,target):
            if (i==1):
                print(self.norm)
                print(target[0].norm)
        add_list=[] #store new rectangle created
        dummy_w=array([]) #store new center value to be compared later
        lower_arr=array([])
        upper_arr=array([])

        for i in range (0,len(self.max_indices)):

            q=zeros((self.dim))
            q[self.max_indices[i]]=self.norm[self.max_indices[i]]/3.0 #split in third the dimension with max indices
            lower_val=target_fnc(self.center_coord-q)
            upper_val=target_fnc(self.center_coord+q)
            dummy_w=append(dummy_w,min(lower_val,upper_val)) #sorting condition
            lower_arr=append(lower_arr,lower_val)
            upper_arr=append(upper_arr,upper_val)

        sort_index=argsort(dummy_w)
        sorted_maxind=self.max_indices[sort_index]
        sorted_lower_arr=lower_arr[sort_index]
        sorted_upper_arr=upper_arr[sort_index]

        for i in range(0,len(sorted_maxind)):

            q=zeros((self.dim))
            q[sorted_maxind[i]]=self.norm[sorted_maxind[i]]/3.0
            divide_arr[sorted_maxind[i]]=3.0

            add_list.append(rect_zone(self.dim,self.norm/divide_arr,self.center_coord-q,sorted_lower_arr[i]))
            add_list.append(rect_zone(self.dim,self.norm/divide_arr,self.center_coord+q,sorted_upper_arr[i]))

        self.norm=self.norm/divide_arr
        self.update_rect()
        
        return add_list

#This class dist_fam regroups all rectangles which have the same dimension, then searching for the minimum and maximum values
#in each group
class dist_fam(): #regroup all member with same vertices distance

    def __init__(self,dist):
        self.dist=dist
        self.index=array([],dtype=int)
        self.value=array([])

    def search_mem(self,full_distlist,full_vallist): #to be optimized further
        for i in range(0,len(full_distlist)):
            if(full_distlist[i]==self.dist):
                self.index=append(self.index,i)
                self.value=append(self.value,full_vallist[i])
    
    
    #DONT search in output list !!!!! fix this part (add center_value in init)
    def search_min(self,output_list):
        dummy=array([])
        #for i in range (0,len(self.index)):
        #    dummy=append(dummy,output_list[self.index[i]].center_value) #fixed
        q=argmin(self.value)
        return (self.value[q],self.index[q]) #return both min value and corresponded index in output_list
    def search_max(self,output_list):
        dummy=array([])
        #for i in range (0,len(self.index)):
        #    dummy=append(dummy,output_list[self.index[i]].center_value)
        q=argmax(self.value)
        return (self.value[q],self.index[q]) #return both max value and corresponded index in output_list
            

#DIRECT method (re-check append method)

def pre_sampling(output_list): #sampling different families,to be further optimized,act on full output_list of rect_zone
    dummy_dist=array([]) #image of rect_zone class index
    dummy_cenval=array([])
    for i in range (0,len(output_list)):
        dummy_dist=append(dummy_dist,output_list[i].dist) #same index as output_list
        dummy_cenval=append(dummy_cenval,output_list[i].center_value)

    dist_list=unique(dummy_dist) #list of distance

    sorted_index=argsort(dist_list)
    fam_list=[]
    sorted_dist_list=dist_list[sorted_index] #sorted list of distance
    
    print ("dummy_dist: "+str(len(dummy_dist)))
    print ("sorteddistlist: "+str(len(sorted_dist_list)))
    for i in range (0,len(sorted_dist_list)):
        q=dist_fam(sorted_dist_list[i])
        q.search_mem(dummy_dist,dummy_cenval) #search index from distance list,image of output_list index
        fam_list.append(q)
    #print ("fam_list "+str(len(fam_list)))
    return fam_list #return a list of families


      
def search_potopt(output_list,min_val,epsilon=1.0e-4):
    part_time=time.time()
    fam_list=pre_sampling(output_list) #already sorted from least to most
    print ("presampling: "+str(time.time()-part_time))
    end_sig=False
    pot_optind=array([],dtype=int)
    #print ("length famlist: "+str(len(fam_list)))
    max_ls=array([])
    min_ls=array([])
    maxid_ls=array([],dtype=int)
    minid_ls=array([],dtype=int)
    for i in range (0,len(fam_list)):
        maxval,maxidval=fam_list[i].search_max(output_list)
        minval,minidval=fam_list[i].search_min(output_list)
        max_ls=append(max_ls,maxval)
        maxid_ls=append(maxid_ls,maxidval)
        min_ls=append(min_ls,minval)
        minid_ls=append(minid_ls,minidval)
        
    
    for i in range (0,len(fam_list)):
        #print ("this is: "+str(i))
        
        
        #print(i)
        #min_i,index_i=fam_list[i].search_min(output_list)
        min_i=min_ls[i]
        index_i=minid_ls[i]
        crit2=(min_i-(min_val-epsilon*abs(min_val)))/fam_list[i].dist
    
        #check search algo
        #start_time=time.time()
        for k in range(0,i):
            if (end_sig==False):
                #max_k,index_k=fam_list[k].search_max(output_list) #check search algo
                max_k=max_ls[k]
                #index_k=maxid_ls[k]
                if ((min_i<max_k)&(crit2<0)): # to check
                    end_sig=True
                    #print("lower")
            else:
                break
        for k in range (i+1,len(fam_list)):
            if (end_sig==False):
                #min_k,index_k=fam_list[k].search_min(output_list)
                min_k=min_ls[k]
                #index_k=minid_ls[k]
                crit=(min_k-min_i)/(fam_list[k].dist-fam_list[i].dist)
                if ((crit<0)|(crit<crit2)):
                    end_sig=True
                    #print("higher")
            else:
                break
        if (end_sig==False):
            pot_optind=append(pot_optind,index_i)   #check very carefully this function  
        end_sig=False
        #print ("analysis: "+str(time.time()-start_time))
    return pot_optind


def append_bar(output_list,add_list):
    dummy=output_list
    for i in range (0,len(add_list)):
        dummy.append(add_list[i])
    return dummy

#*****************************************************************************
#main DIRECT optimization process. The target function is defined in the fnc module

#Input:
#dim: dimension of the function
#T_max: Maximum number of iteration
#epsilon: Stopping criteria
#flag_conv: If True, the optimization stops once the minimum value converges. If False,
#the optimization stops only when T_max is reached

#Output:
#min_val: Minimum value of the function
#min_coord: Coordinates of the function minima
#coord_sampled: Coordinates of all rectangles created during the optimization

def direct_opt(dim,T_max,epsilon=1.0e-4,flag_conv=True): #T_max stands for maximum number of iteration
    process_time=time.time()
    print("starting the optimization process")
    end_sig=False
    coord_sampled=[]    
    dummy_val=array([])
    output_list=[]   
    norm=ones((dim))
    center_coord= 0.5*norm
    output_list.append(rect_zone(dim,norm,center_coord,target_fnc((center_coord))))
    old_min=0.0
    min_val=output_list[0].center_value #initialize the min value
    min_coord=output_list[0].center_coord
    t=0 #counter for iteration
    dummy_coord=[]
    dummy_norm=[]
    add_list=output_list[0].divide_rect()
    output_list=append_bar(output_list,add_list)
    
    for i in range (0,len(add_list)):
        if (dim<=3): #to avoid memory prob, esp for 18 variable optimization
            coord_sampled.append(add_list[i].center_coord)
        dummy_val=append(dummy_val,add_list[i].center_value)
        dummy_coord.append(add_list[i].center_coord)
        dummy_norm.append(add_list[i].norm)
    temp_arg=argmin(dummy_val)
    temp=dummy_val[temp_arg]
    temp_coord=add_list[temp_arg].center_coord 
    if (temp<min_val):
        min_val=temp
        min_coord=temp_coord
    dummy_val=array([])
    t+=1
    
    #main loop
    for t in range(1,T_max):
        if (end_sig==True):
            break
        else:
            start_time=time.time()
            print ("step "+str(t))
            potopt_ind=search_potopt(output_list,min_val,epsilon) 
            print (str(len(potopt_ind))+" is found in "+str(time.time()-start_time))
            
            for i in range (0,len(potopt_ind)): 
                if (end_sig==True):
                    break
                else:
                    part_time=time.time()
                    #have to choose :rather you parallelize the calculation of output_list or the main function itself
                    add_list=output_list[potopt_ind[i]].divide_rect()
                    output_list=append_bar(output_list,add_list)
                    for t in range (0,len(add_list)):
                        if (dim<=3): 
                            coord_sampled.append(add_list[t].center_coord) 
                        dummy_val=append(dummy_val,add_list[t].center_value)
                    temp_arg=argmin(dummy_val)
                    temp=dummy_val[temp_arg]
                    temp_coord=add_list[temp_arg].center_coord
                    
                    if (temp<min_val):
                        old_min=min_val
                        min_val=temp
                        min_coord=temp_coord
                    #end of separation
                        if (abs(min_val-old_min)<=(epsilon*old_min))&(flag_conv==True):
                            end_sig=True
                            print("converge at step: "+str(i))
                    dummy_val=array([])
                    print ("potential optimum number: "+str(i)+" in "+str(time.time()-part_time))
            print ("current min value: "+str(min_val))
            print ("step computing time: "+str(time.time()-start_time)+" s")
    print ("total CPU time: " +str(time.time()-process_time)+" s")
        
    return (min_val,min_coord,coord_sampled) #add norm of the min in return

    
#A simple plot of the center of all rectangle created. Only for two-dimensional function
#coord_sampled input extracted from direct_opt routine
#plot done in normalized coordinates
def plot_test(coord_sampled): #for 2d fnc only
    data=zeros((2,len(coord_sampled)))
    
    for i in range(0,len(coord_sampled)):
        data[0,i]=coord_sampled[i][0]
        data[1,i]=coord_sampled[i][1]
    
    plt.plot(data[0,:],data[1,:],'.b')
    plt.show()
    
    