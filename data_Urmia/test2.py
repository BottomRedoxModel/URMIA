
from netCDF4 import Dataset

fh = Dataset('water_day.nc')

names_vars = [] 
for names,vars in fh.variables.items():
    names_vars.append(names) 
    
if 'z' in names_vars:
    print ('y')
else: 
    print ('no')         
   
#    