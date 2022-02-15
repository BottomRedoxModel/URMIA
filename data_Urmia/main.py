#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Created on 14. des. 2016

@author: E.Protsenko
'''

import os,sys
import numpy as np
from netCDF4 import Dataset 
from PyQt5 import QtGui, QtCore
from PyQt5 import QtWidgets

from matplotlib import rc

from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas)

from matplotlib.backends.backend_qt5agg import (
    NavigationToolbar2QT as NavigationToolbar)


import matplotlib.pyplot as plt

import readdata
import time_plot 
import dist_plot
import fluxes_plot
import all_year_1d
import help_dialog

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
class Window(QtWidgets.QDialog):
    
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)

        # function to display the names of the window flags        
        # Qt.Window Indicates that the widget is a window, 
        # usually with a window system frame and a title bar
        # ! it is not possible to unset this flag if the widget 
        # does not have a parent.
        
        self.setWindowFlags(QtCore.Qt.Window)   
        
        self.setWindowIcon(QtGui.QIcon('bromlogo2.png'))       
        self.figure = plt.figure(figsize=(5.69 , 6.27),
                        facecolor='None',edgecolor='None') 
        self.figure.patch.set_alpha(0)    
            
        # open file system to choose needed nc file 
        
        self.fname ,_  = (QtWidgets.QFileDialog.getOpenFileName(self,
        'Open netcdf ', os.getcwd(), "netcdf (*.nc);; all (*)"))  #str
        
        totitle = os.path.split(self.fname)[1]
        self.setWindowTitle("BROM Pictures ("+str(totitle)+')') 
        #self.totitle = totitle[16:-3]           
        readdata.readdata_brom(self,self.fname)    
         
        # Add group Boxes - boxes of widgets
        createOptionsGroup(self)
        createTimeGroup(self)
        createDistGroup(self) 
                
        # Create widgets
        self.label_choose_var = QtWidgets.QLabel('Choose variable:')                   
        self.qlistwidget = QtWidgets.QListWidget() 

        self.qlistwidget.setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection)
        
        self.all_year_box = QtWidgets.QComboBox()                                      
        self.dist_prof_button = QtWidgets.QPushButton()           
        #self.injlines_checkbox = QtWidgets.QCheckBox(
        #    'Draw inject lines')                     
        self.time_prof_last_year =  QtWidgets.QPushButton()    
        self.time_prof_all =  QtWidgets.QPushButton()  
        
        #self.time_prof_box = QtWidgets.QComboBox()                                   
        self.all_year_button =  QtWidgets.QPushButton()                                  
        self.fick_box = QtWidgets.QPushButton() 
        self.help_button = QtWidgets.QPushButton(' ')
        
        self.fh =  Dataset(self.fname)        
        #def testvar(self):                    
        readdata.read_num_col(self,self.fname)
            #max_num_col = self.testvar.shape[0]
            #return (self.testvar)
            
        # Add group Boxes - boxes of widgets       
        createOptionsGroup(self)
        createTimeGroup(self)
        createDistGroup(self)        
        createCmapLimitsGroup(self)
                              
        if 'i' in self.names_vars:
            self.dist = np.array(self.fh.variables['i'])                          
            self.numcol_2d.setRange(0, int(self.testvar.shape[0]-1))               
            self.numday_box.setRange(0, self.lentime-1)              
            self.numday_stop_box.setRange(0, self.lentime-1)             
            self.numday_stop_box.setValue(self.lentime-1)
        
        self.fh.close()
        
        # Add titles to buttons 
        self.time_prof_all.setText('Time: all year')
        self.fick_box.setText('Fluxes SWI')
        self.all_year_button.setText('1D plot')
        self.time_prof_last_year.setText('Time: last year')               
        self.dist_prof_button.setText('Show Dist Profile')  
           
        ### Define connection between clicking the button and 
        ### calling the function to plot figures         
                                 
        self.time_prof_last_year.released.connect(self.call_print_lyr)
        self.all_year_button.released.connect(self.call_all_year)
        self.time_prof_all.released.connect(self.call_print_allyr)        
        self.fick_box.released.connect(self.call_fluxes)        
        self.dist_prof_button.released.connect(self.call_print_dist)                             
        self.help_button.released.connect(self.call_help)
                  
        self.canvas = FigureCanvas(self.figure)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                        QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)  
        self.canvas.setSizePolicy(sizePolicy)         
        self.toolbar = NavigationToolbar(self.canvas, self)
        #self.canvas.setMinimumSize(self.canvas.size())
        
        ## The QGridLayout class lays out widgets in a grid          
        self.grid = QtWidgets.QGridLayout(self)
        
        readdata.widget_layout(self)        
        readdata.readdata2_brom(self,self.fname)   
                 
        if 'Kz'  in self.names_vars :
            readdata.calculate_ywat(self)
            readdata.calculate_ybbl(self)   
            readdata.y2max_fill_water(self)
            readdata.depth_sed(self)
            readdata.calculate_ysed(self)
            readdata.calculate_ysed(self)
            readdata.calc_nysedmin(self)  
            readdata.y_coords(self)        
        else: 
            self.sediment = False
           
        readdata.colors(self)
        readdata.set_widget_styles(self) 
        
        self.num = 50. 

        self.qlistwidget.addItems(self.sorted_names)
        self.qlistwidget.setFixedSize(
        self.qlistwidget.sizeHintForColumn(0)+ 2 * self.qlistwidget.frameWidth()+50,
              self.canvas.height())
        
    def call_all_year(self): 
          
        all_year_1d.plot(self)
        
    def call_fluxes(self):    
        fluxes_plot.fluxes(self)
        
    def call_print_dist(self): 
        
        dist_plot.dist_profile(self)
            
    def call_print_lyr(self): #last year
        stop = len(self.time)
        start = stop - 365
        time_plot.time_profile(self,start,stop)
        #self.time_profile(start,stop) 
            
    def call_print_allyr(self):  
        start = self.numday_box.value() 
        stop = self.numday_stop_box.value()  
        time_plot.time_profile(self,start,stop)  
                       
    '''def save_figure(self): 
        #does not work 
        printer = QtWidgets.QPrinter(QtWidgets.QPrinter.HighResolution)
        printer.setPageSize(QtWidgets.QPrinter.A9)
        printer.setColorMode(QtWidgets.QPrinter.Color)
        printer.setOutputFormat(QtWidgets.QPrinter.PdfFormat)
        printer.setOutputFileName(self.edit.text())
        self.render(printer)
        #plt.savefig('pdf_fig.pdf',format = 'pdf')    
        #self.figure.savefig('pic.png', format='png')'''
   
   
    def call_help(self):
        help_dialog.show(self) 

        
def createDistGroup(self):  
        
    self.dist_groupBox = QtWidgets.QGroupBox("Distance axis")  
         
    self.dist_grid = QtWidgets.QGridLayout(self.dist_groupBox)
    
    self.col_label = QtWidgets.QLabel('Column: ')
    self.numcol_2d = QtWidgets.QSpinBox() 
    readdata.read_num_col(self,self.fname)
    try:
        self.label_maxcol = QtWidgets.QLabel(
        'max\ncolumn: '+ str(self.testvar.shape[0]-1)) 
    except AttributeError: 
        pass
    #max_col = readdata.read_num_col(self,self.fname)
    #self.label_maxcol_n =
    # QtWidgets.QLabel('max\ncolumn: ') #+ str(testvar.shape[0]-1))      
    
    self.dist_grid.addWidget(self.col_label,0,0,1,1) 
    self.dist_grid.addWidget(self.numcol_2d,1,0,1,1) 
    try:
        self.dist_grid.addWidget(self.label_maxcol,1,1,1,1) 
    except AttributeError: 
        pass    
       
def createTimeGroup(self):  
     
    self.last_year_button = QtWidgets.QPushButton('last year')    
    self.time_groupBox = QtWidgets.QGroupBox(" Time axis")        

    self.numday_start_label = QtWidgets.QLabel('start: ') 
    self.numday_stop_label = QtWidgets.QLabel('stop: ')  
    self.maxday_label = QtWidgets.QLabel('max day: ')
           

         
    self.numday_box = QtWidgets.QSpinBox()     #start
    self.numday_stop_box = QtWidgets.QSpinBox()  
    self.value_maxday = QtWidgets.QLabel(str(self.lentime-1))  
          
    self.cmap_water_box = QtWidgets.QComboBox() 
    self.cmap_water_label = QtWidgets.QLabel('cmap water: ') 
    self.cmap_water_box.addItems(readdata.cmap_list(self))  
    self.cmap_sed_box = QtWidgets.QComboBox()
    self.cmap_sed_box.addItems(readdata.cmap_list(self))
    self.cmap_sed_label = QtWidgets.QLabel('cmap sed: ')
    
         
    self.time_grid = QtWidgets.QGridLayout(self.time_groupBox)   

    #line 0
    self.time_grid.addWidget(self.cmap_water_label,0,0,1,1)
    self.time_grid.addWidget( self.cmap_sed_label,0,1,1,1) 
        
    #line 1
    self.time_grid.addWidget(self.cmap_water_box,1,0,1,1)
    self.time_grid.addWidget( self.cmap_sed_box,1,1,1,1)    
    
    #line 2 
    self.time_grid.addWidget(self.numday_start_label,2,0,1,1)
    self.time_grid.addWidget( self.numday_stop_label,2,1,1,1) 
    self.time_grid.addWidget(      self.maxday_label,2,2,1,1)      
      
    #self.time_grid.addWidget(self.numday_start_label,1,0,1,1)
    
    #line 3             
    self.time_grid.addWidget(     self.numday_box,3,0,1,1) 
    self.time_grid.addWidget(self.numday_stop_box,3,1,1,1)      
    self.time_grid.addWidget(   self.value_maxday,3,2,1,1)   
   
def createOptionsGroup(self):
        self.options_groupBox = QtWidgets.QGroupBox(" Properties ")  
        
        self.scale_all_axes = QtWidgets.QCheckBox(
            "Scale:all columns, all time")                 
        self.yearlines_checkbox = QtWidgets.QCheckBox(
            'Draw year lines')           
        self.datescale_checkbox = QtWidgets.QCheckBox(
            'Format time axis')         
        self.fielddata_checkbox = QtWidgets.QCheckBox(
            'Add field data (1D)') 
        self.interpolate_checkbox = QtWidgets.QCheckBox(
            'Interpolate')       
                     
        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(self.scale_all_axes)
        vbox.addWidget(self.yearlines_checkbox)
        vbox.addWidget(self.datescale_checkbox)
        vbox.addWidget(self.fielddata_checkbox) 
        vbox.addWidget(self.interpolate_checkbox)       
        vbox.addStretch(1)
        self.options_groupBox.setLayout(vbox)     

def createCmapLimitsGroup(self):
        
        self.cmap_groupBox = QtWidgets.QGroupBox("colour map limits ")  
        self.change_limits_checkbox = QtWidgets.QCheckBox('Change limits')
        self.label_maxwater = QtWidgets.QLabel('cmap max water: ')
        
        self.label_minwater = QtWidgets.QLabel('min water: ')   
        self.label_maxsed = QtWidgets.QLabel('cmap max sediment: ')
        self.label_minsed = QtWidgets.QLabel('min sediment: ')  
        self.box_minwater = QtWidgets.QSpinBox()
        self.box_maxwater = QtWidgets.QSpinBox()
        
        self.box_minsed = QtWidgets.QSpinBox()
        self.box_maxsed = QtWidgets.QSpinBox()    
                
        self.box_minwater.setMaximum(1000000000)   
        self.box_minsed.setMaximum(1000000000)           
        self.box_maxwater.setMaximum(1000000000)   
        self.box_maxsed.setMaximum(1000000000)  
        cmap_grid = QtWidgets.QGridLayout(self.cmap_groupBox) 
        
        cmap_grid.addWidget(self.change_limits_checkbox,0,0,1,1) 
        
        cmap_grid.addWidget(self.label_minwater,1,0,1,1)
        cmap_grid.addWidget(self.label_maxwater,1,1,1,1)
        cmap_grid.addWidget(self.box_minwater,2,0,1,1)
        cmap_grid.addWidget(self.box_maxwater,2,1,1,1)  
                 
        cmap_grid.addWidget(self.label_minsed,3,0,1,1)        
        cmap_grid.addWidget(self.label_maxsed,3,1,1,1)
        cmap_grid.addWidget(self.box_minsed,4,0,1,1)
        cmap_grid.addWidget(self.box_maxsed,4,1,1,1)       
                                
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle("plastique")
    main = Window()
    main.setStyleSheet("background-color:#dceaed;")
    main.show()
    sys.exit(app.exec_()) 
