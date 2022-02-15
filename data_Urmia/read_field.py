import csv 
import numpy as np
from PyQt5 import QtGui

def read(self):    
    index = str(self.qlistwidget.currentItem().text())    
    import csv 
    #import numpy as np

    
    #index = str(self.qlistwidget.currentItem().text()) 
    try:
    #    index = str(self.qlistwidget.currentItem().text())
    #except AttributeError: 
    #    print ("Choose the variable to print ")    
        with open('field_data\BBL_18-08-2009.txt','r') as f:
                # important to specify delimiter right 
            reader = csv.reader(f,delimiter = '\t',dialect='excel')
            r = []    
            for row in reader:
                # if you don't know which delimiter is used 
                # print one row to view it 
                #print (row)
                #break
                r.append(row)  
                
        r1 = np.transpose(np.array(r[:]))  
        r1 = r1.astype(float)
                          
        sed_depth_1808  = r1[1]
        depth_1808 = sed_depth_1808/100 + 9
        pH_1808 = r1[2]    
        H2S_1808 = r1[3]    
        Alk_1808 = r1[4]   
        PO4_1808 = r1[5]     
        SO4_1808 = r1[6]      
        DON_1808 = r1[7]      
        NH4_1808 = r1[8]      
        NO3_1808 = r1[9]     
        NO2_1808 = r1[10]     
        Fe_1808 = r1[11]      
        Mn_1808 = r1[12]      
        Ni_1808 = r1[13]      
        Hgtot_1808 = r1[14]     
        MeHg_1808 = r1[15]  
    
        f.close() 
           
        with open('field_data\BBL_15-07-2009.txt','r') as f:
            reader = csv.reader(f,delimiter = '\t',dialect='excel')
            r = []    
            for row in reader:
                r.append(row)  
                
        r1 = np.transpose(np.array(r[:]))  
        r1 = r1.astype(float)
                          
        depth_1507  = r1[0]
        sed_depth_1507  = r1[1]
        #depth_1507 = sed_depth_1507/100 + 9
        pH_1507 = r1[2]    
        H2S_1507 = r1[3]    
        Alk_1507 = r1[4]   
        PO4_1507 = r1[5]     
        SO4_1507 = r1[6]      
        DON_1507 = r1[7]      
        NH4_1507 = r1[8]      
        NO3_1507 = r1[9]     
        NO2_1507 = r1[10]     
        Fe_1507 = r1[11]      
        Mn_1507 = r1[12]      
        Ni_1507 = r1[13]      
        Hgtot_1507 = r1[14]     
        MeHg_1507 = r1[15]        
        f.close()
          
        with open('field_data/bbl_porewater_2009.txt','r') as f:
            reader = csv.reader(f,delimiter = '\t',dialect='excel')
            r = []    
            for row in reader:
                r.append(row)  
                
        r1 = np.transpose(np.array(r[:]))  
        r1 = r1.astype(float)
        
        depth_2009  = r1[0]
        sed_depth_2009  = r1[1]
        pH_2009 = r1[2]    
        H2S_2009 = r1[3]    
        Alk_2009 = r1[4]   
        PO4_2009 = r1[5]     
        SO4_2009 = r1[6]      
        DON_2009 = r1[7]      
        NH4_2009 = r1[8] 
        Fe_2009 = r1[9]  
        Mn_2009 = r1[10]  
        Ni_2009 = r1[11]  
        Hgtot_2009 = r1[12]  
        MeHg_2009 = r1[13]  
        f.close()
    
        with open('field_data/bbl_porewater_2010.txt','r') as f:
            reader = csv.reader(f,delimiter = '\t',dialect='excel')
            r = []    
            for row in reader:
                #print (row)
                #break
                r.append(row)  
                
        r1 = np.transpose(np.array(r[:]))  
        r1 = r1.astype(float)
        #np.array(r1[7]).astype(np.double)
        depth_2010  = np.array(r1[0]).astype(float)
        sed_depth_2010  = np.array(r1[1]).astype(float)
        pH_2010 = np.array(r1[2]).astype(np.float)    
        H2S_2010 = np.array(r1[3]).astype(float)    
        Alk_2010 = np.array(r1[4]).astype(float)   
        PO4_2010 = np.array(r1[5]).astype(float)    
        DON_2010 = np.array(r1[6]).astype(float)      
        NH4_2010 = np.array(r1[7]).astype(float)
        NO3_2010 = np.array(r1[8]).astype(float)
        SO4_2010 = np.array(r1[9]).astype(float)     
        Fe_2010 = np.array(r1[10]).astype(float)  
        Mn_2010 = np.array(r1[11]).astype(float)  
        Ni_2010 = np.array(r1[12]).astype(float)  
     
        f.close()
                
        if index == 'Alk':
            self.ax00.plot(Alk_1808,depth_1808,'ro-',zorder = 10)
            self.ax00.plot(Alk_1507,depth_1507,'ro-',zorder = 10)
    
                                
            self.ax10.plot(Alk_1808,depth_1808,'ro-',zorder = 10)
            self.ax10.plot(Alk_1507,depth_1507,'ro-',zorder = 10)       
            self.ax10.plot(Alk_2009,depth_2009,'ro-',zorder = 10)   
            self.ax10.plot(Alk_2010,depth_2010,'ro-',zorder = 10)  
            
            self.ax20.plot(Alk_2009,sed_depth_2009,'ro-',zorder = 10) 
            self.ax20.plot(Alk_2010,sed_depth_2010,'ro-',zorder = 10) 
                         
        elif index == 'pH':
                    
            self.ax00.plot(pH_1808,depth_1808,'ro-',zorder = 10)
            self.ax00.plot(pH_1507,depth_1507,'ro-',zorder = 10)
                         
            self.ax10.plot(pH_1808,depth_1808,'ro-',zorder = 10)
            self.ax10.plot(pH_1507,depth_1507,'ro-',zorder = 10) 
            self.ax10.plot(pH_2009,depth_2009,'ro-',zorder = 10)   
            self.ax10.plot(pH_2010,depth_2010,'ro-',zorder = 10)
                    
            self.ax20.plot(pH_2009,sed_depth_2009,'ro-',zorder = 10) 
            self.ax20.plot(pH_2010,sed_depth_2010,'ro-',zorder = 10) 
                    
        elif index == 'H2S':
                      
            self.ax10.plot(H2S_1808,depth_1808,'ro-',zorder = 10)
            self.ax10.plot(H2S_1507,depth_1507,'ro-',zorder = 10)         
            self.ax10.plot(H2S_2009,depth_2009,'ro-',zorder = 10)   
            self.ax10.plot(H2S_2010,depth_2010,'ro-',zorder = 10)
                    
            self.ax20.plot(H2S_2009,sed_depth_2009,'ro-',zorder = 10) 
            self.ax20.plot(H2S_2010,sed_depth_2010,'ro-',zorder = 10)       
            
        elif index == 'PO4':
                    
            self.ax00.plot(PO4_1808,depth_1808,'ro-',zorder = 10)
            self.ax00.plot(PO4_1507,depth_1507,'ro-',zorder = 10)
                         
            self.ax10.plot(PO4_1808,depth_1808,'ro-',zorder = 10)
            self.ax10.plot(PO4_1507,depth_1507,'ro-',zorder = 10) 
            self.ax10.plot(PO4_2009,depth_2009,'ro-',zorder = 10)   
            self.ax10.plot(PO4_2010,depth_2010,'ro-',zorder = 10)
                    
            self.ax20.plot(PO4_2009,sed_depth_2009,'ro-',zorder = 10) 
            self.ax20.plot(PO4_2010,sed_depth_2010,'ro-',zorder = 10)         
            
        elif index == 'DON':
                    
            self.ax00.plot(DON_1808,depth_1808,'ro-',zorder = 10)
            self.ax00.plot(DON_1507,depth_1507,'ro-',zorder = 10)
                         
            self.ax10.plot(DON_1808,depth_1808,'ro-',zorder = 10)
            self.ax10.plot(DON_1507,depth_1507,'ro-',zorder = 10) 
            self.ax10.plot(DON_2009,depth_2009,'ro-',zorder = 10)   
            self.ax10.plot(DON_2010,depth_2010,'ro-',zorder = 10)
                    
            self.ax20.plot(DON_2009,sed_depth_2009,'ro-',zorder = 10) 
            self.ax20.plot(DON_2010,sed_depth_2010,'ro-',zorder = 10)          
            
        elif index == 'NH4':
                    
            self.ax00.plot(NH4_1808,depth_1808,'ro-',zorder = 10)
            self.ax00.plot(NH4_1507,depth_1507,'ro-',zorder = 10)
                         
            self.ax10.plot(NH4_1808,depth_1808,'ro-',zorder = 10)
            self.ax10.plot(NH4_1507,depth_1507,'ro-',zorder = 10) 
            self.ax10.plot(NH4_2009,depth_2009,'ro-',zorder = 10)   
            self.ax10.plot(NH4_2010,depth_2010,'ro-',zorder = 10)
                    
            self.ax20.plot(NH4_2009,sed_depth_2009,'ro-',zorder = 10) 
            self.ax20.plot(NH4_2010,sed_depth_2010,'ro-',zorder = 10)         
            
        elif index == 'NO3':
                    
            self.ax00.plot(NO3_1808,depth_1808,'ro-',zorder = 10)
            self.ax00.plot(NO3_1507,depth_1507,'ro-',zorder = 10)
                         
            self.ax10.plot(NO3_1808,depth_1808,'ro-',zorder = 10)
            self.ax10.plot(NO3_1507,depth_1507,'ro-',zorder = 10)   
            self.ax10.plot(NO3_2010,depth_2010,'ro-',zorder = 10)
    
            self.ax20.plot(NO3_2010,sed_depth_2010,'ro-',zorder = 10) 
            
        elif index == 'SO4':
                    
            self.ax00.plot(SO4_1808,depth_1808,'ro-',zorder = 10)
            self.ax00.plot(SO4_1507,depth_1507,'ro-',zorder = 10)
                         
            self.ax10.plot(SO4_1808,depth_1808,'ro-',zorder = 10)
            self.ax10.plot(SO4_1507,depth_1507,'ro-',zorder = 10) 
            self.ax10.plot(SO4_2009,depth_2009,'ro-',zorder = 10)   
            self.ax10.plot(SO4_2010,depth_2010,'ro-',zorder = 10)
                    
            self.ax20.plot(SO4_2009,sed_depth_2009,'ro-',zorder = 10) 
            self.ax20.plot(SO4_2010,sed_depth_2010,'ro-',zorder = 10) 
                    
        elif index == 'Fe2':
                    
            self.ax00.plot(Fe_1808,depth_1808,'ro-',zorder = 10)
            self.ax00.plot(Fe_1507,depth_1507,'ro-',zorder = 10)
                         
            self.ax10.plot(Fe_1808,depth_1808,'ro-',zorder = 10)
            self.ax10.plot(Fe_1507,depth_1507,'ro-',zorder = 10) 
            self.ax10.plot(Fe_2009,depth_2009,'ro-',zorder = 10)   
            self.ax10.plot(Fe_2010,depth_2010,'ro-',zorder = 10)
                    
            self.ax20.plot(Fe_2009,sed_depth_2009,'ro-',zorder = 10) 
            self.ax20.plot(Fe_2010,sed_depth_2010,'ro-',zorder = 10)     
            
        elif index == 'Mn2':
                    
            self.ax00.plot(Mn_1808,depth_1808,'ro-',zorder = 10)
            self.ax00.plot(Mn_1507,depth_1507,'ro-',zorder = 10)
                         
            self.ax10.plot(Mn_1808,depth_1808,'ro-',zorder = 10)
            self.ax10.plot(Mn_1507,depth_1507,'ro-',zorder = 10) 
            self.ax10.plot(Mn_2009,depth_2009,'ro-',zorder = 10)   
            self.ax10.plot(Mn_2010,depth_2010,'ro-',zorder = 10)
                    
            self.ax20.plot(Mn_2009,sed_depth_2009,'ro-',zorder = 10) 
            self.ax20.plot(Mn_2010,sed_depth_2010,'ro-',zorder = 10)   
            
        elif index == 'Ni':
                    
            self.ax00.plot(Ni_1808,depth_1808,'ro-',zorder = 10)
            self.ax00.plot(Ni_1507,depth_1507,'ro-',zorder = 10)
                         
            self.ax10.plot(Ni_1808,depth_1808,'ro-',zorder = 10)
            self.ax10.plot(Ni_1507,depth_1507,'ro-',zorder = 10) 
            self.ax10.plot(Ni_2009,depth_2009,'ro-',zorder = 10)   
            self.ax10.plot(Ni_2010,depth_2010,'ro-',zorder = 10)
                    
            self.ax20.plot(Ni_2009,sed_depth_2009,'ro-',zorder = 10) 
            self.ax20.plot(Ni_2010,sed_depth_2010,'ro-',zorder = 10)         
                  
        elif index == 'Hg0' or index == 'Hg2':
                    
            self.ax00.plot(Hgtot_1808,depth_1808,'ro-',zorder = 10)
            self.ax00.plot(Hgtot_1507,depth_1507,'ro-',zorder = 10)
                         
            self.ax10.plot(Hgtot_1808,depth_1808,'ro-',zorder = 10)
            self.ax10.plot(Hgtot_1507,depth_1507,'ro-',zorder = 10) 
            self.ax10.plot(Hgtot_2009,depth_2009,'ro-',zorder = 10)   
                    
            self.ax20.plot(Hgtot_2009,sed_depth_2009,'ro-',zorder = 10) 
            
        elif index == 'MeHg':
                    
            self.ax00.plot(MeHg_1808,depth_1808,'ro-',zorder = 10)
            self.ax00.plot(MeHg_1507,depth_1507,'ro-',zorder = 10)
                         
            self.ax10.plot(MeHg_1808,depth_1808,'ro-',zorder = 10)
            self.ax10.plot(MeHg_1507,depth_1507,'ro-',zorder = 10) 
            self.ax10.plot(MeHg_2009,depth_2009,'ro-',zorder = 10)   
                    
            self.ax20.plot(MeHg_2009,sed_depth_2009,'ro-',zorder = 10) 
          
    except IOError:
        QtGui.QMessageBox.about(self, "Retry",
                                 'You do not have field data \n'  
                 'or the files are in wrong format')
        pass                          