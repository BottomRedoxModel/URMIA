#!/usr/bin/python
# -*- coding: utf-8 -*-
# this ↑ comment is important to have 
# at the very first line 
# to define using unicode 

#import matplotlib


from matplotlib.figure import Figure
import os,sys
import numpy as np

from PyQt4 import QtGui, QtCore
from PyQt4.QtGui import QSpinBox,QLabel,QComboBox

from matplotlib import rc
from matplotlib import style
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas)
from matplotlib.backends.backend_qt4agg import (
    NavigationToolbar2QT as NavigationToolbar)
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.gridspec as gridspec

import readdata


#for i in xrange(5):
class Window(QtGui.QDialog):
    def __init__(self, parent=None, ):
        super(Window, self).__init__(parent)

        # open file system to choose needed nc file 
        fname = str(QtGui.QFileDialog.getOpenFileName(self,
        'Open netcdf ', os.getcwd(), "netcdf (*.nc);; all (*)"))   

        # call submodule to read the nc file 
        readdata.readdata_brom(self,fname) 
        readdata.colors(self)        
        readdata.calculate_ywat(self)
        readdata.calculate_ybbl(self) 
        readdata.depth_sed(self)   #calc depths in cm 
        # kz is determined at midpoints, should be 
        # calculated separately *self.depth2
        readdata.depth_sed2(self)  #calc depths in cm for kz         
        readdata.y2max_fill_water(self)                
        readdata.calculate_ysed(self)
        readdata.y_coords(self)  
                      
        num_years = float(len(self.time)/365.)    

        qapp = QtGui.QApplication(sys.argv)
        self.qwidget = QtGui.QWidget()
        self.qwidget.setGeometry(QtCore.QRect(30, 30, 1500, 800)) # w,h
        
        #QHBoxLayout lines up widgets horizontally
        #QHBoxLayout(QWidget *parent)       
        qlayout = QtGui.QHBoxLayout(self.qwidget)
        self.qwidget.setLayout(qlayout)
        
        
        #QScrollArea provides a scrolling view onto another widget.
        self.qscroll = QtGui.QScrollArea(self.qwidget)
        self.qscroll.setGeometry(QtCore.QRect(30, 30,1500,800))
        #self.qscroll.setFrameStyle(QtGui.QFrame.NoFrame)
        qlayout.addWidget(self.qscroll)
        
        self.qscrollContents = QtGui.QWidget()
        
        #QVBoxLayout lines up widgets vertically
        #QVBoxLayout(QWidget *parent)
        self.qscrollLayout = QtGui.QVBoxLayout(self.qscrollContents)
        self.qscrollLayout.setGeometry(QtCore.QRect(10,10, 1500, 700))
        
        self.qscroll.setWidget(self.qscrollContents)
        self.qscroll.setWidgetResizable(True)

        
        self.qfigWidget = QtGui.QWidget(self.qscrollContents)
        print num_years 
        
        xx = 2.0*num_years
        print  xx
        self.figure = Figure((xx, 1.0), dpi=100,facecolor='white')
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setParent(self.qfigWidget)
        toolbar = NavigationToolbar(self.canvas, self.qfigWidget)
        #axes = self.figure.add_subplot(111)
        #axes.plot([1,2,3,4])
        
    
        # Create widget 1                    
        self.time_prof_box = QtGui.QComboBox(self.qscrollContents)
        # add items to Combobox
        for i in self.var_names_profile:
            self.time_prof_box.addItem(str(i))
        # Define connection between clicking the button and 
        # calling the function to plot figures        
        self.time_prof_box.currentIndexChanged.connect(
            self.time_profile) 
        
        
        
        
        
        
        
        # place plot components in a layout
        plotLayout = QtGui.QVBoxLayout()
        plotLayout.addWidget(toolbar)
        plotLayout.addWidget(self.time_prof_box)
        plotLayout.addWidget(self.canvas)
        self.qfigWidget.setLayout(plotLayout)
        
        # prevent the canvas to shrink beyond a point
        # original size looks like a good minimum size
        self.canvas.setMinimumSize(self.canvas.size())    
        self.qscrollLayout.addWidget(self.qfigWidget)
        
        self.qwidget.show()
        exit(qapp.exec_()) 
        
    def time_profile(self):

        plt.clf()

        #self.figure = Figure((7.0, 5.0), dpi=100,facecolor='white')
        for n in range(1,len(self.var_names_profile)):
            if (self.time_prof_box.currentIndex() == n) :
                z = self.vars[n] #self.po4
                title = self.time_prof_box.currentText() 
                self.num_var = n  
                                 
        gs = gridspec.GridSpec(2, 1) 
        gs.update(left=0.06, right=0.92,top = 0.91,bottom = 0.07,
                   wspace=0.2,hspace=0.1)         
        x = self.time #np.arange(6)

        #print num_years
        
        #self.xsize = 8.27 
        #self.figure = plt.figure(figsize=(5 ,self.xsize))       
              
        # change size of figure       
 
        xsize = int(500.*len(self.time)/365.)
        
                     
        #self.qfigWidget.setGeometry(QtCore.QRect(30, 30, 15000 , 700)) #w,h
        #self.qwidget.setGeometry(QtCore.QRect(10,10, 500 , 500)) #w,h   
        #self.figure.set_size_inches(30., 6.) #w,h      
        #self.qscrollLayout.setGeometry(QtCore.QRect(10,10, 900, 900))
        #self.qscroll.setGeometry(QtCore.QRect(10, 10, 1200, 700))  #gray area 
                        
        y = self.depth #np.arange(5)
        y_sed = self.depth_sed
        
        z1 = np.array(z).flatten()
        zz = z1.reshape((len(x),len(y))).T
        
        watmin = readdata.varmin(self,zz,0) #0 - water 
        watmax = readdata.varmax(self,zz,0)
        
        sed_min = readdata.varmin(self,zz,1)
        sed_max = readdata.varmax(self,zz,1)  

        #set constant min and max for pH
        if self.num_var == 34: 
            watmax = 9
            watmin = 6.8
            sed_max = 10
            sed_min = 6 

        X,Y = np.meshgrid(x,y)
        X_sed,Y_sed = np.meshgrid(x,y_sed)
        ax = self.figure.add_subplot(gs[0])
        ax2 = self.figure.add_subplot(gs[1])        
        #self.figure, (ax, ax2) = plt.subplots(2, 1, sharex=True)
        #f.set_size_inches(11.69,8.27)

        ax.set_title(title)
        ax.set_ylim(self.y1max-1,0)       
        ax2.set_ylim(self.ysedmax,self.y3min) #ysedmin
        xlen = len(x)
        ax.set_xlim(0,xlen)


        ax.set_ylabel('Depth (m)',fontsize= self.font_txt) #Label y axis 
        ax2.set_ylabel('Depth (cm)',fontsize= self.font_txt) 
        ax2.set_xlabel('Number of day',fontsize= self.font_txt) 
            
      
        
        self.num = 50.            
        wat_levs = np.linspace(watmin,watmax,num= self.num)
        sed_levs = np.linspace(sed_min,sed_max,
                             num = self.num)
                
        int_wat_levs = []
        int_sed_levs= []
                
        for n in wat_levs:
            n = readdata.int_value(self,n,watmin,watmax)
            int_wat_levs.append(n)            
        for n in sed_levs:
            n = readdata.int_value(self,n,sed_min,
                                   sed_max)
            int_sed_levs.append(n)            
                      
        #define color maps 
        cmap = plt.cm.jet #gnuplot#jet#gist_rainbow
        cmap1 = plt.cm.rainbow  

        ## contourf() draw contour lines and filled contours
        # levels = A list of floating point numbers indicating 
        # the level curves to draw, in increasing order    
        # If None, the first value of Z will correspond to the lower
        # left corner, location (0,0).  
        # If image, the rc value for image.origin will be used.    
        CS = ax.contourf(X,Y, zz, levels= int_wat_levs,
                              cmap=cmap)

        CS1 = ax2.contourf(X_sed,Y_sed, zz, levels = int_sed_levs,
                              cmap=cmap1) #, origin='lower'
       
        # Add an axes at position rect [left, bottom, width, height]
        cax = self.figure.add_axes([0.94, 0.53, 0.02, 0.35])                
        cax1 = self.figure.add_axes([0.94, 0.1, 0.02, 0.35])
        
        wat_ticks = readdata.ticks(watmin,watmax)        
        sed_ticks = readdata.ticks(sed_min,sed_max)
        
        cb = plt.colorbar(CS,cax = cax,ticks = wat_ticks)
        cb_sed = plt.colorbar(CS1,cax = cax1 )
        #cb.set_label('Water')
        

        #cb.set_ticks(wat_ticks)
        cb_sed.set_ticks(sed_ticks)  
        ax2.axhline(0, color='white', linestyle = '--',linewidth = 1 )
        for n in range(0, len(self.time)) :
            if np.remainder(n, 365) == 0 :
                ax2.axvline(n, color='white', linestyle = '--',linewidth = 1)
                ax.axvline(n, color='white', linestyle = '--',linewidth = 1)                
            else :
                pass
                
        #self.qwidget.show()
        self.canvas.draw() 
              
        
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    app.setStyle("plastique")
    main = Window()#None, QtCore.Qt.WindowSystemMenuHint | 
                    #     QtCore.Qt.WindowTitleHint)#QtGui.QWidget()
    #main.setWindowFlags(
    #    main.windowFlags() | QtCore.Qt.FramelessWindowHint)

    
    #setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
    #main = QtGui.QDialog
    main.setStyleSheet("background-color:#d8c9c2;")
    # int posx, int posy, int w, int h
    #main.setGeometry(50, 50, 900,500)
    main.show()
    #PySide.QtCore.Qt.WindowFlags
    sys.exit(app.exec_()) 
        
        