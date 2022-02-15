'''
Created on 19. july 2017

@author: ELP
'''

import os,sys
import netCDF4
from netCDF4 import Dataset
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import Tkinter as tk

from Tkinter import *
from tkFileDialog import askopenfilename # python 2
root = tk.Tk()
root.withdraw()
#from tkinter import * python3 
#from tkinter.tkFiledialog import askopenfilename # for python 3
 
import time

#plt.style.use('ggplot')
plt.style.use('bmh')

#[u'dark_background', u'bmh', u'grayscale', u'ggplot', u'fivethirtyeight']
#3demo('default')
#demo('seaborn')

with open('wat_hardan2.txt', 'r') as f:
    # important to specify delimiter right 
    reader = csv.reader(f, delimiter=',')
    r = []
    for row in reader:
        # if you don't know which delimiter is used 
        # print one row to view it 
        #print (row)
        #break
        r.append(row)        
r1 = np.transpose(np.array(r[1:])) 
#print (r1[0])
st,nh4  = r1[0], r1[1]    
#nh4 = r1[1]  
po4 = r1[2]
ptot = r1[3]
ntot = r1[4]
si = r1[5]  
no3 = r1[6]  
o2 = r1[7]  
doc = r1[8]  
pressure = r1[9]  
fe3 = r1[10]  
mn2 = r1[11]  
dop = r1[12]  
don = r1[13] 

tic = [2190.56, 2021.71, 988.85, 999.57]
tic_depth = [470,200,0, 10]
alk = [2301.98,2180.317, 1059.57, 1058.218]

start_list = [0,7,12,17,23,28] 
end_list = [7,12,17,23,28,33]
st_list = [1,2,3,4,5,6]

#Cruise data from sediments 
depth200 = [0.50, 1.20, 2.50 , 4.00, 7.5]
depth199 = [0.50, 1.50, 2.50 , 4.00, 7.5]

si200 = [285.7,357.1,328.6,300.0,314.3]
no3200 = [35.71428571, 16.42857143, 9.64, 8.92,2.857142857 ]
po4200 = [9.68, 9.68, 9.68, 8.39, 6.45]
tic200 = [2916.7, 3841.7, 3683.3, 3675.0, 4033.3]
alk200 = [2860, 3460, 3410, 3410, 3850]
nh4200 = [121.4285714,100.,100.,85.71428571,71.42857143]

si199 = [189.3, 171.4, 214.3, 175.0, 221.4]
no3199 = [3.928, 4.071, 5.571, 2.071, 2.714] 
po4199 = [41.94,30.97, 21.94,24.52, 10.32]
alk199 = [3070, 3050, 3590, 4040, 3800]
tic199 = [3141.7, 3641.7, 4191.7, 4133.3, 3916.7]
nh4199 = [100.,114.2857143,145.7142857,107.1428571,92.85714286]


ask_filename = askopenfilename(
    initialdir= os.getcwd(),
    filetypes =(("NetCDF file", "*.nc"),
                ("All Files","*.*")),
                title = "Choose a file.")
  
fname = ask_filename # "BROM_Hfjord_out.nc"
#os.path.split(ask_filename)[1] 

fh =  Dataset(fname)

kz =  np.array(fh.variables['Kz'][:,:,0])
depth_brom = np.array(fh.variables['z'][:])
sed = np.array(fh.variables['z'][15]) # depth of the SWI
sed_depth_brom = np.array(fh.variables['z'][:])-sed
sed_depth_brom = sed_depth_brom *100
si_brom = np.array(fh.variables['Si'][:,:,:]) 
nh4_brom = np.array(fh.variables['NH4'][:,:,:])
po4_brom = np.array(fh.variables['PO4'][:,:,:])
no3_brom = np.array(fh.variables['NO3'][:,:,:])
o2_brom = np.array(fh.variables['O2'][:,:,:])
alk_brom = np.array(fh.variables['Alk'][:,:,:])
tic_brom = np.array(fh.variables['DIC'][:,:,:])
time_brom = np.array(fh.variables["time"][:])

len_time = len(time_brom)
i_brom = np.array(fh.variables['i'][:]) 
len_i = len(i_brom)

#def create_figure():
#figure1 = plt.figure(figsize=(5,7), dpi=100)
#gs = gridspec.GridSpec(2,1) #yx
#ax = figure1.add_subplot(gs[0])
#ax1 = figure1.add_subplot(gs[1])



# plot cruise data water  
def plot_water(var,axis):
    if var == alk or var == tic: 
        axis.scatter(var,tic_depth,color = '#825572',
                     zorder=10,edgecolors='k')
        pass
    else: 
        for m in range(0,6):        
            axis.plot(var[start_list[m]:end_list[m]], 
                    pressure[start_list[m]:end_list[m]],
                    'o-', color = '#9f6b8c')  #,
                    #'o-',alpha = 0.7)        
            #axis.annotate(st_list[m],
            #(var[start_list[m]],pressure[start_list[m]]),
            #              textcoords='offset points', 
            #xytext=(3, 3))
            
            #axis.scatter(var[33],pressure[33],
            #             c = '#8b0000',s = 300,
            #             zorder = 10)        
                 
def plot_sed(var,var_depth,axis): 
    if var == None :
        pass
    else:    
        axis.plot(var,var_depth,'o-')
        
def plot_brom(var_brom,axis):    
    for day in range(0,len_time,20):  # plot every 20'th day  
        for n in range(0,len_i,2): # plot every 2d column
                #for m in range(0,100):# time_brom-1:
            axis.plot(var_brom[day,:15,n],depth_brom[:15],
                      color = '#a1a7af',alpha = 0.8,
                      linewidth = 0.2) 
                   
def plot_brom_sed(var_brom,axis):  
    for day in range(0,len_time,20):  
        for n in range(0,len_i,2):
            axis.plot(var_brom[day,13:,n],sed_depth_brom[13:],
                      color = '#a1a7af',alpha = 0.3,
                      linewidth = 0.2,zorder = 1 )          

def plot_all(var_brom,var_water,var_sed1,var_sed2,title):    
    #ax1.cla()
    figure1 = plt.figure(figsize=(5,7), dpi=100)
    gs = gridspec.GridSpec(2,1) #yx

    ax = figure1.add_subplot(gs[0])
    ax1 = figure1.add_subplot(gs[1])
    ax.set_ylim(310,0)
    ax1.set_ylim(5,-5) 
    
    ax.axhspan(310, 0,color='#cce6f5',
                 alpha = 0.4, label = "water")     
    ax1.axhspan(0, -5,color='#cce6f5',
                 alpha = 0.4, label = "water")    
    ax1.axhspan(5, 0,color='#dcc196',
                 alpha = 0.4, label = "sediment")        
     
    #plt.gcf().clear()
    plot_brom(var_brom,ax)
    plot_brom_sed(var_brom,ax1)            
    plot_water(var_water,ax)
    plot_sed(var_sed1,depth200,ax1)
    plot_sed(var_sed2,depth199,ax1)
    ax.set_title(title)
    script_dir = os.path.dirname(__file__)
    
    results_dir = os.path.join(script_dir, 'Results/')
    
    if not os.path.isdir(results_dir):
       os.makedirs(results_dir)
       
    plt.savefig(results_dir+title+'.png')
    #plt.show()
    plt.clf()
    
# plt.savefig("/home/dir/subdir/filename.png")   
plot_all(alk_brom,alk,alk199,alk200,'Alkalinity')  
plot_all(si_brom,si,si199,si200,'Si')
plot_all(no3_brom,no3,no3199,no3200,'NO3')
plot_all(po4_brom,po4,po4199,po4200,'PO4')
plot_all(o2_brom,o2,None,None,'O2')
  
plot_all(tic_brom,tic,tic199,tic200,'TIC')
plot_all(nh4_brom,nh4,nh4199,nh4200,'NH4')
