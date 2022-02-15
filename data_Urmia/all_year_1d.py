#!/usr/bin/python
# -*- coding: utf-8 -*-
# this ↑ comment is important to have 
# at the very first line 
# to define using unicode 


'''
Created on 30. jun. 2017

@author: ELP
'''


import matplotlib.pyplot as plt
from PyQt5 import QtGui,QtWidgets
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter

import readdata
import read_field


def plot(self): 
               
    plt.clf()        
    try:
        index = str(self.qlistwidget.currentItem().text())
    except AttributeError: 
        print ("Choose the variable to print ")       
        messagebox = QtWidgets.QMessageBox.about(self, "Retry",
                                             'Choose variable,please') 
        return None 
     
    start = self.numday_box.value() 
    stop = self.numday_stop_box.value() 
    #index = str(self.time_prof_box.currentText())
    #print ('test all year', index) 
    data_units = self.fh.variables[index].units                
    self.figure.patch.set_facecolor('white') 
    gs = gridspec.GridSpec(3,1) 
    gs.update(left=0.3, right=0.7,top = 0.94,bottom = 0.04,
               wspace=0.2,hspace=0.3) 
    
    self.ax00 = self.figure.add_subplot(gs[0]) # water         
    self.ax10 = self.figure.add_subplot(gs[1]) # bbl
    self.ax20 = self.figure.add_subplot(gs[2]) # sediment 
    
    for axis in (self.ax00,self.ax10,self.ax20):
        axis.yaxis.grid(True,'minor')
        axis.xaxis.grid(True,'major')                
        axis.yaxis.grid(True,'major')    
                     
    numcol = self.numcol_2d.value() # 
    
    # read chosen variable 
    z = np.array(self.fh.variables[index])
    z = np.array(z[:,:,numcol]) 
    
    self.ax00.set_title(index +', ' + data_units) 
    
    #Label y axis        
    self.ax00.set_ylabel('h, m', 
                    fontsize= self.font_txt) 
    self.ax10.set_ylabel('h, m', 
                    fontsize= self.font_txt)   
    self.ax20.set_ylabel('h, cm', 
                    fontsize= self.font_txt)
    
    self.ax00.set_ylim(self.y1max,0) 
    self.ax00.axhspan(self.y1max,0,color='#dbf0fd',
                 alpha = 0.7,label = "water" )
     
    self.ax10.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    self.ax10.set_ylim(self.y2max, self.y1max)   
    self.ax10.axhspan(self.y2max, self.y1max,color='#c5d8e3',
                 alpha = 0.4, label = "bbl"  )                
    
    
    self.ax20.set_ylim(self.ysedmax, self.ysedmin) 

    self.ax20.axhspan(self.ysedmin,0,
                 color='#c5d8e3',alpha = 0.4,
                 label = "bbl"  )        
    self.ax20.axhspan(self.ysedmax,0,
                 color='#b08b52',alpha = 0.4,
                 label = "sediment"  )
    
    
    for n in range(start,stop,10):#365
        """if (n>0 and n <60) or (n>=335 and n<365) : #"winter"
        #if n >= 0 and n<=60 or n >= 335 and n <365 : #"winter"                               
            ax00.plot(z[n][0:self.ny2max],
                  self.depth[0:self.ny2max],
                  self.wint,alpha = self.a_w, 
                  linewidth = self.linewidth , zorder = 10) 
         
            ax10.plot(z[n][0:self.ny2max],
                  self.depth[0:self.ny2max],
                  self.wint,alpha = self.a_w, 
                  linewidth = self.linewidth , zorder = 10) 
        
            ax20.plot(z[n][self.nysedmin-1:],
                  self.depth_sed[self.nysedmin-1:],
                  self.wint, alpha = self.a_w,
                  linewidth = self.linewidth, zorder = 10) """  
        #else: 
        self.ax00.plot(z[n][0:self.ny2max],
              self.depth[0:self.ny2max],
              self.spr_aut,alpha = self.a_w, 
              linewidth = self.linewidth , zorder = 8) 

        self.ax10.plot(z[n][0:self.ny2max],
              self.depth[0:self.ny2max],
              self.spr_aut,alpha = self.a_w, 
              linewidth = self.linewidth , zorder = 8) 
    
        self.ax20.plot(z[n][self.nysedmin-1:],
              self.depth_sed[self.nysedmin-1:],
              self.spr_aut, alpha = self.a_w,
              linewidth = self.linewidth, zorder = 8)                          
        #ax20.scatter(z[n][self.nysedmin-1:],
        #      self.depth_sed[self.nysedmin-1:]) 
   
    if self.fielddata_checkbox.isChecked():
        #print ('is checked')
        read_field.read(self) 
        
 

                      
    self.canvas.draw()     



'''def all_year_charts(self): 
    #messagebox = QtGui.QMessageBox.about(self, "Next time",
    #                                     'it does not work yet =(')           
    plt.clf()
    gs = gridspec.GridSpec(3,3) 
    gs.update(left=0.06, right=0.93,top = 0.94,bottom = 0.04,
               wspace=0.2,hspace=0.1)   
    self.figure.patch.set_facecolor('white') 
    #self.figure.patch.set_facecolor(self.background) 
    #Set the background color  
    self.ax00 = self.figure.add_subplot(gs[0]) # water         
    self.ax10 = self.figure.add_subplot(gs[1]) # water
    self.ax20 = self.figure.add_subplot(gs[2]) # water 

    self.ax01 = self.figure.add_subplot(gs[3])          
    self.ax11 = self.figure.add_subplot(gs[4])
    self.ax21 = self.figure.add_subplot(gs[5])

    ax02 = self.figure.add_subplot(gs[6])    
    ax12 = self.figure.add_subplot(gs[7])
    ax22 = self.figure.add_subplot(gs[8])

    ax00.set_ylabel('Depth (m)',fontsize= self.font_txt) #Label y axis
    ax01.set_ylabel('Depth (m)',fontsize= self.font_txt)   
    ax02.set_ylabel('Depth (cm)',fontsize= self.font_txt) 
                                 
    for n in range(1,len(self.vars)):
        if (self.all_year_1d_box.currentIndex() == n) :
            
            varname1 = self.vars[n][0] 
            varname2 = self.vars[n][1] 
            varname3 = self.vars[n][2] 
            #print (n)
            z123 = readdata.read_all_year_var(self,
                        self.fname,varname1,varname2,varname3)

            z0 = np.array(z123[0])
            z1 = np.array(z123[1])
            z2 = np.array(z123[2])
            
            ax00.set_title(str(self.titles_all_year[n][0]), 
            fontsize=self.xlabel_fontsize, fontweight='bold') 
            
            ax10.set_title(str(self.titles_all_year[n][1]), 
            fontsize=self.xlabel_fontsize, fontweight='bold') 
            
            ax20.set_title(str(self.titles_all_year[n][2]), 
            fontsize=self.xlabel_fontsize, fontweight='bold')                                 
            self.num_var = n  

    for axis in (ax00,ax10,ax20,ax01,ax11,ax21,ax02,ax12,ax22):
        #water          
        axis.yaxis.grid(True,'minor')
        axis.xaxis.grid(True,'major')                
        axis.yaxis.grid(True,'major') 
                
    ax00.set_ylim(self.y1max,0)   
    ax10.set_ylim(self.y1max,0)  
    ax20.set_ylim(self.y1max,0) 
    
    ax01.set_ylim(self.y2max, self.y2min)   
    ax11.set_ylim(self.y2max, self.y2min)  
    ax21.set_ylim(self.y2max, self.y2min) 

    ax02.set_ylim(self.ysedmax, self.ysedmin)   
    ax12.set_ylim(self.ysedmax, self.ysedmin)  
    ax22.set_ylim(self.ysedmax, self.ysedmin) 
    #
    #n0 = self.varmax(self,z0,1) #[0:self.y2max_fill_water,:].max() 
    start = 0
    stop = 365 
    #### to change""!!!!
    
    
    watmin0 = readdata.varmin(self,z0,"wattime",start,stop) 
    watmin1 = readdata.varmin(self,z1,"wattime",start,stop) 
    watmin2 = readdata.varmin(self,z2,"wattime",start,stop)          

    watmax0 = readdata.varmax(self,z0,"wattime",start,stop) 
    watmax1 = readdata.varmax(self,z1,"wattime",start,stop)
    watmax2 = readdata.varmax(self,z2,"wattime",start,stop)  
             
    sed_min0 = readdata.varmin(self,z0,"sedtime",start,stop) 
    sed_min1 = readdata.varmin(self,z1,"sedtime",start,stop) 
    sed_min2 = readdata.varmin(self,z2,"sedtime",start,stop)    

    sed_max0 = readdata.varmax(self,z0,"sedtime",start,stop) 
    sed_max1 = readdata.varmax(self,z1,"sedtime",start,stop)         
    sed_max2 = readdata.varmax(self,z2,"sedtime",start,stop)         
    
    if self.num_var == 5: #pH 
        watmax1 = 9
        watmin1 = 6.5
    elif self.num_var == 2: #po4, so4
        watmax0 = 3  
        #watmax1 = 7000.          
        #watmin1 = 4000.            
    else:
        pass
                
    
    self.m0ticks = readdata.ticks(watmin0,watmax0)
    self.m1ticks = readdata.ticks(watmin1,watmax1)
    self.m2ticks = readdata.ticks(watmin2,watmax2)  
    
    self.sed_m0ticks = readdata.ticks(sed_min0,sed_max0)
    self.sed_m1ticks = readdata.ticks(sed_min1,sed_max1)
    self.sed_m2ticks = readdata.ticks(sed_min2,sed_max2)                 
    #for axis in (ax00,ax10,ax20):             
    
    ax00.set_xlim(watmin0,watmax0)   
    ax01.set_xlim(watmin0,watmax0)         
    ax02.set_xlim(sed_min0,sed_max0)
    
    ax10.set_xlim(watmin1,watmax1)   
    ax11.set_xlim(watmin1,watmax1)         
    ax12.set_xlim(sed_min1,sed_max1)         
     
    ax20.set_xlim(watmin2,watmax2)   
    ax21.set_xlim(watmin2,watmax2)         
    ax22.set_xlim(sed_min2,sed_max2) 
            
    ax10.set_xlim(watmin1,watmax1)   
    ax11.set_xlim(watmin1,watmax1)         
    #ax12.set_xlim(sed_min1,sed_max1) 
                 
    ax20.set_xlim(watmin2,watmax2)   
    ax21.set_xlim(watmin2,watmax2)         
    #ax22.set_xlim(sed_min2,sed_max2)                  
    #water

                 
    ax00.fill_between(
                    self.m0ticks, self.y1max, 0,
                    facecolor= self.wat_col1, alpha=0.1 ) #self.a_w
    ax01.fill_between(
                    self.m0ticks, self.y2min_fill_bbl ,self.y2min,
                    facecolor= self.wat_col1, alpha=0.1 ) #self.a_w    
    ax01.fill_between(self.m0ticks, self.y2max, self.y2min_fill_bbl,
                           facecolor= self.bbl_col1, alpha=self.a_bbl) 
        
    ax02.fill_between(self.sed_m0ticks,self.ysedmin_fill_sed,-10,
                           facecolor= self.bbl_col1, alpha=self.a_bbl)          
    ax02.fill_between(self.sed_m0ticks, self.ysedmax, self.ysedmin_fill_sed,
                           facecolor= self.sed_col1, alpha=self.a_s)          
    
        #axis.fill_between(self.xticks, self.y2max, self.y2min_fill_bbl,
        #                   facecolor= self.bbl_color, alpha=self.alpha_bbl)        
        
    ax10.fill_between(
                    self.m1ticks, self.y1max, 0,
                    facecolor= self.wat_col1, alpha=0.1 ) #self.a_w
    
    ax11.fill_between(
                    self.m1ticks, self.y2min_fill_bbl ,self.y2min,
                    facecolor= self.wat_col1, alpha=0.1 ) #self.a_w    
    ax11.fill_between(self.m1ticks, self.y2max, self.y2min_fill_bbl,
                           facecolor= self.bbl_col1, alpha=self.a_bbl)      
    ax12.fill_between(self.sed_m1ticks,self.ysedmin_fill_sed,-10,
                           facecolor= self.bbl_col1, alpha=self.a_bbl) 
    ax12.fill_between(self.sed_m1ticks, self.ysedmax, self.ysedmin_fill_sed,
                          facecolor= self.sed_col1, alpha=self.a_s)        
    ax20.fill_between(
                    self.m2ticks, self.y1max, 0,
                    facecolor= self.wat_col1, alpha=0.1 ) #self.a_w
    
    ax21.fill_between(
                    self.m2ticks, self.y2min_fill_bbl ,self.y2min,
                    facecolor= self.wat_col1, alpha=0.1 ) #self.a_w    
    ax21.fill_between(self.m2ticks, self.y2max, self.y2min_fill_bbl,
                           facecolor= self.bbl_col1, alpha=self.a_bbl)     
    ax22.fill_between(self.sed_m2ticks,self.ysedmin_fill_sed,-10,
                           facecolor= self.bbl_col1, alpha=self.a_bbl)                         
    ax22.fill_between(self.sed_m2ticks, self.ysedmax, self.ysedmin_fill_sed,
                           facecolor= self.sed_col1, alpha=self.a_s) 
                        

    
            
    for n in range(0,3): #365
        if n >= 0 and n<=60 or n >= 335 and n <365 : #"winter" 
            linewidth = self.linewidth
                              
            ax00.plot(z0[n],self.depth,self.wint,alpha = 
                      self.a_w, linewidth = linewidth , zorder = 10) 
            ax10.plot(z1[n],self.depth,self.wint,alpha = 
                      self.a_w, linewidth = linewidth , zorder = 10)
            ax20.plot(z2[n],self.depth,self.wint,alpha = 
                      self.a_w, linewidth = linewidth, zorder = 10 )  
            
            ax01.plot(z0[n],self.depth,self.wint,alpha = 
                      self.a_w, linewidth = linewidth, zorder = 10 ) 
            ax11.plot(z1[n],self.depth,self.wint,alpha = 
                      self.a_w, linewidth = linewidth , zorder = 10)
            ax21.plot(z2[n],self.depth,self.wint,alpha = 
                      self.a_w, linewidth = linewidth, zorder = 10 ) 

            ax02.plot(z0[n],self.depth_sed,self.wint,alpha = 
                      self.a_w, linewidth = linewidth, zorder = 10 ) 
            ax12.plot(z1[n],self.depth_sed,self.wint,alpha = 
                      self.a_w, linewidth = linewidth, zorder = 10 )
            ax22.plot(z2[n],self.depth_sed,self.wint,alpha = 
                      self.a_w, linewidth = linewidth, zorder = 10 ) 
        elif n >= 150 and n < 249: #"summer"
            ax00.plot(z0[n],self.depth,self.summ,alpha = 
                      self.a_s, linewidth = linewidth, zorder = 10 ) 
            ax10.plot(z1[n],self.depth,self.summ,alpha = 
                      self.a_s, linewidth = linewidth, zorder = 10 )
            ax20.plot(z2[n],self.depth,self.summ,alpha = 
                      self.a_s, linewidth = linewidth, zorder = 10 )  
            
            ax01.plot(z0[n],self.depth,self.summ,alpha = 
                      self.a_s, linewidth = linewidth, zorder = 10 ) 
            ax11.plot(z1[n],self.depth,self.summ,alpha = 
                      self.a_s, linewidth = linewidth, zorder = 10 )
            ax21.plot(z2[n],self.depth,self.summ,alpha = 
                      self.a_s, linewidth = linewidth, zorder = 10 ) 

            ax02.plot(z0[n],self.depth_sed,self.summ,alpha = 
                      self.a_s, linewidth = linewidth, zorder = 10 ) 
            ax12.plot(z1[n],self.depth_sed,self.summ,alpha = 
                      self.a_s, linewidth = linewidth, zorder = 10 )
            ax22.plot(z2[n],self.depth_sed,self.summ,alpha = 
                      self.a_s, linewidth = linewidth, zorder = 10 ) 
        else : #"autumn and spring"
            ax00.plot(z0[n],self.depth,self.spr_aut,alpha = 
                      self.a_aut, linewidth = linewidth, zorder = 10 ) 
            ax10.plot(z1[n],self.depth,self.spr_aut,alpha = 
                      self.a_aut, linewidth = linewidth, zorder = 10 )
            ax20.plot(z2[n],self.depth,self.spr_aut,alpha = 
                      self.a_aut, linewidth = linewidth, zorder = 10 )  
            
            ax01.plot(z0[n],self.depth,self.spr_aut,alpha = 
                      self.a_aut, linewidth = linewidth, zorder = 10 ) 
            ax11.plot(z1[n],self.depth,self.spr_aut,alpha = 
                      self.a_aut, linewidth = linewidth, zorder = 10 )
            ax21.plot(z2[n],self.depth,self.spr_aut,alpha = 
                      self.a_aut, linewidth = linewidth, zorder = 10 ) 

            ax02.plot(z0[n],self.depth_sed,self.spr_aut,
                      alpha = self.a_aut, zorder = 10) 
            ax12.plot(z1[n],self.depth_sed,self.spr_aut,
                      alpha = self.a_aut, zorder = 10)
            ax22.plot(z2[n],self.depth_sed,self.spr_aut,
                      alpha = self.a_aut, zorder = 10)      
 
    self.canvas.draw()  ''' 
    
      
