! This file is part of Bottom RedOx Model (BROM, v.1.1).
! BROM is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the BROM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Shamil Yakubov, 
!                     Elizaveta Protsenko, Phil Wallhead
!----------------------------------------------------------------------- 
    
    
    
    module ids
    
    use io_ascii, only: find_index
    use fabm_types, only: attribute_length
    
    implicit none
    
    !Indices of state variables that are needed in brom-transport and subroutines (i.e. boundary conditions description)
    integer :: id_O2, id_Mn2, id_Mn3, id_Mn4, id_H2S, id_Fe2, id_Fe3, id_FeS,   &
               id_MnS, id_DON, id_PON, id_NH4, id_NO2, id_NO3, id_S0, id_S2O3,  &
               id_SO4, id_DIC, id_Alk, id_pCO2, id_PO4, id_Si, id_Sipart, id_Phy, id_Het, &
               id_Baae, id_Bhae, id_Baan, id_Bhan, id_Hplus, id_CaCO3, id_FeS2, id_MnCO3, &
               id_BaSO4, id_Cl, id_Na, id_NaCl, id_V_air, id_V_wat, id_V_sed, &
               id_dV_ch,id_dV_sink, &
               id_hz_m, id_z_m, id_area,id_Na2CaSO42, id_K2Ca2MgSO44, id_CaMgCO32, id_CaSO4, &
               id_KCl, id_Mg, id_K, id_Ca_u, id_MgSO4 
    
    
    contains
    

!======================================================================================================================= 
    subroutine get_ids(par_name)
    
    !Input variables
    character(len=attribute_length),dimension(:)  :: par_name

    id_O2  = find_index(par_name, 'B_BIO_O2')      
    id_DON = find_index(par_name, 'B_BIO_DON')         
    id_PON = find_index(par_name, 'B_BIO_PON')       
    id_Phy = find_index(par_name, 'B_BIO_Phy')        
    id_Het = find_index(par_name, 'B_BIO_Het')          
    id_Baae = find_index(par_name, 'B_BACT_Baae')           
    id_Bhae = find_index(par_name, 'B_BACT_Bhae')          
    id_Baan = find_index(par_name, 'B_BACT_Baan')          
    id_Bhan = find_index(par_name, 'B_BACT_Bhan')    
    id_NO3 = find_index(par_name, 'B_NUT_NO3')             
    id_NO2 = find_index(par_name, 'B_NUT_NO2')            
    id_NH4 = find_index(par_name, 'B_NUT_NH4')         
    id_PO4 = find_index(par_name, 'B_NUT_PO4')      
    id_Si = find_index(par_name, 'B_NUT_Si')      
    id_Sipart = find_index(par_name, 'B_Si_Sipart')          
    id_Mn4 = find_index(par_name, 'B_Mn_Mn4')           
    id_Mn2 = find_index(par_name, 'B_Mn_Mn2')          
    id_Mn3 = find_index(par_name, 'B_Mn_Mn3')          
    id_MnS = find_index(par_name, 'B_Mn_MnS')
    id_MnCO3 = find_index(par_name, 'B_Mn_MnCO3')       
    id_Fe3 = find_index(par_name, 'B_Fe_Fe3')          
    id_Fe2 = find_index(par_name, 'B_Fe_Fe2')          
    id_FeS = find_index(par_name, 'B_Fe_FeS')          
    id_FeS2 = find_index(par_name, 'B_Fe_FeS2')
    id_SO4 = find_index(par_name, 'B_S_SO4')          
    id_S2O3 = find_index(par_name, 'B_S_S2O3')          
    id_S0 = find_index(par_name, 'B_S_S0')           
    id_H2S = find_index(par_name, 'B_S_H2S')         
    id_DIC = find_index(par_name, 'B_C_DIC')          
    id_Alk = find_index(par_name, 'B_C_Alk')  
    id_Hplus = find_index(par_name, 'B_pH_Hplus')          
    id_CaCO3 = find_index(par_name, 'B_Ca_CaCO3')
    id_Na = find_index(par_name, 'B_halite_Na')
    id_Cl = find_index(par_name, 'B_halite_Cl')
    id_NaCl = find_index(par_name, 'B_halite_NaCl')
    id_V_air = find_index(par_name, 'B_volumes_V_air')
    id_V_wat = find_index(par_name, 'B_volumes_V_wat')
    id_V_sed = find_index(par_name, 'B_volumes_V_sed')
    id_dV_ch = find_index(par_name, 'B_volumes_dV_ch')
    id_dV_sink = find_index(par_name, 'B_volumes_dV_sink')
    id_hz_m = find_index(par_name, 'B_volumes_hz_m')
    id_z_m = find_index(par_name, 'B_volumes_z_m')
    id_area = find_index(par_name, 'B_volumes_area')
    id_Na2CaSO42 = find_index(par_name, 'B_minerals_Na2CaSO42')
    id_K2Ca2MgSO44 = find_index(par_name, 'B_minerals_K2Ca2MgSO44')
    id_CaMgCO32 = find_index(par_name, 'B_minerals_CaMgCO32')
    id_CaSO4 = find_index(par_name, 'B_minerals_CaSO4')
    id_KCl = find_index(par_name, 'B_minerals_KCl')
    id_Mg = find_index(par_name, 'B_minerals_Mg')
    id_K = find_index(par_name, 'B_minerals_K')
    id_Ca_u = find_index(par_name, 'B_minerals_Ca_u')
    id_MgSO4 = find_index(par_name, 'B_minerals_MgSO4')

    end subroutine get_ids
!======================================================================================================================= 
    
    
    end module ids
