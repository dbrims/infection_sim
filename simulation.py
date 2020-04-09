#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 15:22:22 2020

@author: dbrims
"""

import pandas as pd
import random
import numpy as np
import matplotlib.pyplot as plt
import time

#convert the inputed percent strings to decimal 
def make_pct(str):
    if '%' in str:
           pct=float(str[0:-1])/100
    else:
        pct=float(str)
        if pct>1:
            pct/=100

    return pct

#plot a color coated scatter plot of all the balls, IU need to figure
#how to annimate
def scatter_plot(df):
    col = df.infected.map({False:'b', True:'r'})
    df.plot.scatter(x='x coord', y='y coord',c=col, figsize=(10,7.5))

    return 
    
#shift all the balls each iteration based on their velocities, if 
#the balls go past the wall, asign them to the index of the wall    
def moving(df):

    df['x coord']=round(df['x coord']+df['Vx'],0)
    df['y coord']=round(df['y coord']+df['Vy'],0)
    df['x coord'].loc[df['x coord']<0]=0
    df['x coord'].loc[df['x coord']>100]=100
    df['y coord'].loc[df['y coord']<0]=0
    df['y coord'].loc[df['y coord']>100]=100

    return df
    
#if the ball hits the wall, it will bounce off based on 
#angle of incidence = angle of reflection (ie, the the sign of the 
#velocity normal to the wall is flipped)
def bounce_wall(df):

    df.loc[df['x coord']==0, 'Vx']=-1*df['Vx']
    df.loc[df['x coord']==100, 'Vx']=-1*df['Vx']
    df.loc[df['y coord']==0, 'Vy']=-1*df['Vy']
    df.loc[df['y coord']==100, 'Vy']=-1*df['Vy']

    
    return df

#I am ignoring the three body problem and only considering pairwise collisions     
#if two balls collide, we assume a perfectly elastic collision. we assume
# the mass of the balls are equivelent and magnitude of the velocity is 
#unchange.  the resulting x and y components of the velocity are 
#depend end on the initial velocity, angle of travel of the particle
#and collision angle of the two balls. the equation is given here:
#https://en.wikipedia.org/wiki/Elastic_collisiong 

    
def bounce_ball(df, infect_prob, speed):

    
    balls=collisions(df)

    
    for ball in balls:
        Vx1=df.loc[ball[0],'Vx']
        Vx2=df.loc[ball[1],'Vx']
        Vy1=df.loc[ball[0],'Vy']
        Vy2=df.loc[ball[1],'Vy']
        theta1=df.loc[ball[0],'Angle']
        theta2=df.loc[ball[1],'Angle']
        st1=df.loc[ball[0],'stationary']
        st2=df.loc[ball[1],'stationary']
    
        phi=get_phi(Vx1, Vx2, Vy1, Vy2)
        
        if st1==False and st2==False:
            Vx1_new=speed*(np.cos(theta2-phi)*np.cos(phi)+np.sin(theta1-phi)*np.cos(phi+np.pi/2))    
            Vx2_new=speed*(np.cos(theta1-phi)*np.cos(phi)+np.sin(theta2-phi)*np.cos(phi+np.pi/2))
            Vy1_new=speed*(np.cos(theta2-phi)*np.sin(phi)+np.sin(theta1-phi)*np.sin(phi+np.pi/2))   
            Vy2_new=speed*(np.cos(theta1-phi)*np.sin(phi)+np.sin(theta2-phi)*np.sin(phi+np.pi/2))
            
            df= update_V(df, ball, Vx1_new, Vx2_new, Vy1_new, Vy2_new)
            
           
        elif st1==True:
            Vx1_new=speed*(np.cos(theta2-phi)*np.cos(phi))
            Vx2_new=speed*(np.sin(theta2-phi)*np.cos(phi+np.pi/2)) 
            Vy1_new=speed*(np.cos(theta2-phi)*np.sin(phi))   
            Vy2_new=speed*(np.sin(theta2-phi)*np.sin(phi+np.pi/2))
            
            df= update_V(df, ball, Vx1_new, Vx2_new, Vy1_new, Vy2_new)
            
        elif st2==True:
            Vx1_new=speed*(np.sin(theta1-phi)*np.cos(phi+np.pi/2))  
            Vx2_new=speed*(np.cos(theta1-phi)*np.cos(phi))
            Vy1_new=speed*(np.sin(theta1-phi)*np.sin(phi+np.pi/2))  
            Vy2_new=speed*(np.cos(theta1-phi)*np.sin(phi))
            
            df= update_V(df, ball, Vx1_new, Vx2_new, Vy1_new, Vy2_new)          
    
    return df

#funcion to update the data frame with the new velocity terms coming from 
# the bounce_ball function
def update_V(df, ball, Vx1_new, Vx2_new, Vy1_new, Vy2_new):

    df.loc[ball[0],'Vx']=Vx1_new
    df.loc[ball[0],'Vy']=Vy1_new
    df.loc[ball[1],'Vx']=Vx2_new
    df.loc[ball[1],'Vy']=Vy2_new

    return df
    
#function which finds which balls are colliding, based on their being
# at the same grid coordinate and returns those balls as a list
def collisions(df):
    
    coords= list(zip(df['x coord'], df['y coord']))

    l=len(coords)
    balls=[]
    for x in range(l):
        for y in range(x,l):
            if x!=y:
                if coords[x]==coords[y]:
                    balls.append([x,y])
                    
    return balls
    
#calculating the collision angle between two colliding balls
def get_phi(vx1, vx2, vy1, vy2):
    dx =vx1-vx2
    dy = vy1-vy2
    if dx ==0:
        phi = np.pi/2.
    else:
        phi = np.arctan(dy/dx) 
    return phi
    
#if an infected ball stricks another ball, based on the probability of
# infection we assign the colliding ball either infected or not    
def infecting(df, infect_prob,radius):  
    if radius==0:
        balls=collisions(df)
    else:
        balls=infect_radius(df,radius)  
          
    for ball in balls:
        if df.loc[ball[0],'infected']==True:      
            if random.uniform(0,1)<=infect_prob:
                df.loc[ball[1],'infected']=True           
        if df.loc[ball[1],'infected']==True:        
            if random.uniform(0,1)<=infect_prob:
                df.loc[ball[0],'infected']=True               
    return df

def infect_radius(df, radius):
        
    coords= list(zip(df['x coord'], df['y coord']))

    l=len(coords)
    close_balls=[]
    for x in range(l):
        for y in range(x,l):
            if x!=y:
                if ((coords[x][0]-coords[y][0])**2+(coords[x][1]-coords[y][1])**2)<=radius**2:
                    close_balls.append([x,y])
    return close_balls
    

    
#if we want to input initial condition manually, these are the inputs    
#density=int(input('what is the density of balls you want (1-1000)?:'))
#stationary_str=input('what percent of balls do you want stationary?:')
#infect_prob_str=input('what is the prbability of infection?:')
#speed=int(input('what speed are the particles moving?:'))
#radius=float(input('what radius of infection do you want?:'))
    
#stationary=make_pct(stationary_str)
#infect_prob=make_pct(infect_prob_str)

#we want to do a sensitity analysis, itterating through parameters
# to see what gives the greated effect on infection spread
# I should do at least 5 trials, and average the results but not sure
#how to create that dataframe
    
#for trials in range(5):
cycle=0
for infectivity in range(10,100,10): 
    if cycle==0:
        density_time = pd.DataFrame(columns = ['Density',f'{infectivity}'+'% infectivity'])
    else:
        density_time[f'{infectivity}'+'% infectivity']=0   
#

#for static in range(0,50,10):
#    if cycle==0:
#        density_time = pd.DataFrame(columns = ['Density',f'{static}'+'% static'])
        
#    else:
 #       density_time[f'{static}'+'% static']=0
    

    step=0
    for density in range(100,1000,100):
        speed=5
#        infect_prob=1
        radius=5
        stationary=0
 #       stationary=(static)/100
        infect_prob=(infectivity)/100
    
        location=[]
        velocity=[]
        infected=[]

#building a data frame with the initial state of each ball
#balls are randomly assinged non-overlapping x and y postions within the box
        for ball in range(density):
            infected.append(False) 
            xy=[random.randrange(1,99),random.randrange(1,99)]
            test=False
            while test==False:
                if xy in location:
                    xy=[random.randint(1,99),random.randint(1,99)]
                else:
                    location.append(xy)
                    test=True
#initial velocities are generated and stationary balls are assigned                
            if ball<(density-stationary*density):        
                sign=random.randint(0,1)
                Vx=random.uniform(-1,1)
                if sign==0:
                    Vy=np.sqrt(speed**2-Vx**2)
                    velocity.append([Vx,Vy, False])
                elif sign==1:
                    Vy=np.sqrt(speed**2-Vx**2)*-1
                    velocity.append([Vx,Vy, False])
            else:
                velocity.append([0,0, True])
    
#based on the initial conditions, the the angle in radians of travel
# of the ball is computed and the intial infected ball is assigned
        balls_p= pd.DataFrame(location, columns=['x coord', 'y coord'])
        balls_v= pd.DataFrame(velocity, columns=['Vx','Vy', 'stationary'])
        balls_df=pd.concat([balls_p, balls_v],axis=1)
        balls_df['Angle']=np.arctan(balls_df['Vy']/balls_df['Vx'])
        balls_df['infected']=infected
        balls_df.loc[0,'infected']=True

#the time it takes for 2/3 of the balls to becom infected is determined
#by measuring the itterations required to achieve that level of infection
#iteration involves moving each ball, determing if a ball has collided with
#a wall, a ball has collided with another ball, and if an infected ball has 
#collided with another ball, spreading the disease
        i=0
        while (balls_df.infected[balls_df.infected==True].count()/(balls_df.infected.count()))<.67:
        #scatter_plot(balls_df)
        #time.sleep(2)
            balls_df=moving(balls_df)
            balls_df=bounce_wall(balls_df)
            balls_df=bounce_ball(balls_df,infect_prob,speed)
            balls_df=infecting(balls_df, infect_prob, radius)
            i+=1
            print(balls_df.infected[balls_df.infected==True].count()/(balls_df.infected.count()))

#create a pandas df of the simulation results. which is then plotted
        if cycle==0:
            
            den_list = [density, i]
#           den_list=[density, i]
            den_list
            length = len(density_time)
            density_time.loc[length] = den_list
#        else:
#            density_time.loc[density,f'{static}'+'% static']=i
#            step+=1
        else:
            density_time.loc[density,f'{infectivity}'+'% infectivity']=i
            step+=1
        print(density_time)
    if cycle==0:
        density_time.set_index('Density', inplace=True)
    cycle+=1
    density_time.plot()
    
#file='trial'+f'{trials}'+'.csv'
#density_time.to_csv(rfile, index = False, header=True)
    
#for trial in range(5):
#    if trial==0:
#        file='trial'+f'{trial}'+'.csv'
#        density_time.read_csv(file, index = False, header=True)
#    else:
#        file='trial'+f'{trial}'+'.csv'
#        df.read_csv(file, index = False, header=True)
#    density_time=density_time.add(df)
#density_time=density_time.div(5)



            
density_time.plot()

                

        


    

            
    


