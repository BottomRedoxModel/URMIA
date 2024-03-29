!-----------------------------------------------------------------------
! BROM2 is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the BROM2 distribution.
!-----------------------------------------------------------------------

#include "fabm_driver.h"

module fabm_niva_brom_bio
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public :: type_niva_brom_bio
    !variables allocated here
    type(type_state_variable_id):: id_Phy,id_Het
    type(type_state_variable_id):: id_O2,id_PON,id_DON
    !state dependencies
    type(type_state_variable_id):: id_NO2,id_NO3,id_NH4,id_PO4
    type(type_state_variable_id):: id_Baae,id_Baan,id_Bhae,id_Bhan
    type(type_state_variable_id):: id_DIC,id_H2S,id_Si,id_Sipart,id_Alk

    type(type_diagnostic_variable_id):: id_DcPM_O2,id_DcDM_O2
    type(type_diagnostic_variable_id):: id_MortHet,id_Grazing,id_RespHet
    type(type_diagnostic_variable_id):: id_GrazBhae,id_GrazBhan
    type(type_diagnostic_variable_id):: id_GrazBaae,id_GrazBaan
    type(type_diagnostic_variable_id):: id_GrazPhy,id_GrazPOP,id_GrazBact
    type(type_diagnostic_variable_id):: id_MortPhy,id_ExcrPhy,id_LimNH4
    type(type_diagnostic_variable_id):: id_LimN,id_GrowthPhy
    type(type_diagnostic_variable_id):: id_LimT,id_LimP,id_LimNO3,id_LimSi
    type(type_diagnostic_variable_id):: id_LimLight, id_N_fixation
    type(type_diagnostic_variable_id):: id_eO2mO2,id_osat

    type(type_dependency_id):: id_temp,id_salt,id_par,id_pres
    type(type_dependency_id):: id_Hplus
    type(type_horizontal_dependency_id):: id_windspeed

    !Model parameters
    !specific rates of biogeochemical processes
    real(rk):: K_DON_ox,K_PON_ox,K_PON_DON,Tda,beta_da,K_omox_o2
    !----Phy  ----------!
    real(rk):: K_phy_gro,k_Erlov,Iopt
    real(rk):: K_phy_mrt,K_phy_exc,LatLight
    integer :: phy_t_dependence ! select dependence on T: (1) Old; (2) for Arctic; (3) ERSEM
    !----Het -----------!
    real(rk):: K_het_phy_gro,K_het_phy_lim,K_het_pom_gro,K_het_pom_lim,K_het_bac_gro
    real(rk):: K_het_res,K_het_mrt,Uz,Hz,limGrazBac
    !---- O2--------!
    !Upper boundary, for oxygen flux calculations
    !SY - pvel?
    real(rk):: pvel = 5._rk ! wind speed [m/s]
    real(rk):: a0 = 31.25_rk !oxygen saturation [uM]
    real(rk):: a1 = 14.603_rk !oxygen saturation [-]
    real(rk):: a2 = 0.4025_rk !oxygen saturation [1/degC]
    !---- N, P, Si--!
    real(rk):: K_nox_lim,K_nh4_lim,K_psi,K_nfix,K_po4_lim,K_si_lim
    !---- Sinking---!
    real(rk):: Wsed,Wphy,Whet
    !---- Stoichiometric coefficients ----!
    real(rk):: r_n_p, r_o_n, r_c_n, r_si_n
  contains
    procedure :: initialize
    procedure :: do
    procedure :: do_surface
  end type
contains
  !
  !
  !
  subroutine initialize(self,configunit)
    class (type_niva_brom_bio), intent(inout), target :: self
    integer,                    intent(in)            :: configunit

    call self%get_parameter(&
         self%LatLight,'LatLight','degree','Latitude',default=50.0_rk)
    call self%get_parameter(&
         self%K_DON_ox,'K_DON_ox','[1/day]',&
         'Specific rate of oxidation of DON with O2',&
         default=0.01_rk)
    call self%get_parameter(&
         self%K_PON_ox,'K_PON_ox','[1/day]',&
         'Specific rate of oxidation of PON with O2',&
         default=0.002_rk)
    call self%get_parameter(&
         self%K_PON_DON, 'K_PON_DON', '[1/day]',&
         'Specific rate of Autolysis of PON to DON',&
         default=0.1_rk)
    call self%get_parameter(&
         self%beta_da,'beta_da','[1/day]',&
         'Temperature control coefficient for OM decay',&
         default=20.0_rk)
    call self%get_parameter(&
         self%K_omox_o2,'K_omox_o2','[uM]',&
         'half sat. of o2 for OM mineralization',&
         default=1.0_rk)
    call self%get_parameter(&
         self%Tda,'Tda','[1/day]',&
         'Temperature control coefficient for OM decay',&
         default=13.0_rk)
    !----Phy----------!
    call self%get_parameter(&
         self%K_phy_gro,'K_phy_gro','1/d','Maximum specific growth rate',&
         default=2.0_rk)
    call self%get_parameter(&
         self%k_Erlov,'k_Erlov', '1/m','Extinction coefficient',&
         default=0.05_rk)
    call self%get_parameter(&
         self%Iopt,'Iopt','Watts/m**2/h','Optimal irradiance',&
         default=25.0_rk)
    call self%get_parameter(&
         self%K_phy_mrt,'K_phy_mrt','1/d','Specific rate of mortality',&
         default=0.10_rk)
    call self%get_parameter(&
         self%K_phy_exc,'K_phy_exc','1/d','Specific rate of excretion',&
         default=0.01_rk)
    call self%get_parameter(&
         self%phy_t_dependence,'phy_t_dependence','-','T dependence fro Phy growth',&
         default=1)
    !----Het----------!
    call self%get_parameter(&
         self%K_het_phy_gro,'K_het_phy_gro','1/d',&
         'Max.spec. rate of grazing of Het on Phy',&
         default=1.0_rk)
    call self%get_parameter(&
         self%K_het_phy_lim,'K_het_phy_lim','nd',&
         'Half-sat.const.for grazing of Het on Phy for Phy/Het ratio',&
         default=1.1_rk)
    call self%get_parameter(&
         self%K_het_pom_gro,'K_het_pom_gro','mmol/m**3',&
         'Max.spec.rate of grazing of Het on POM',&
         default=0.70_rk)
    call self%get_parameter(&
         self%K_het_bac_gro,'K_het_bac_gro','mmol/m**3',&
         'Max.spec.rate of grazing of Het on POM',&
         default=0.70_rk) 
    call self%get_parameter(&
         self%K_het_pom_lim,'K_het_pom_lim','nd',&
         'Half-sat.const.for grazing of Het on POM for POM/Het ratio',&
         default=0.2_rk)
    call self%get_parameter(&
         self%K_het_res,'K_het_res','1/d',&
         'Specific respiration rate',&
         default=0.02_rk)
    call self%get_parameter(&
         self%K_het_mrt,'K_het_mrt','1/d',&
         'Maximum specific rate of mortality of Het',&
         default=0.05_rk)
    call self%get_parameter(&
         self%Uz,'Uz','nd',&
         'Food absorbency for Het',&
         default=0.5_rk)
    call self%get_parameter(&
         self%Hz,'Hz','nd',&
         'Ratio betw. diss. and part. excretes of Het',&
         default=0.5_rk)
    call self%get_parameter(&
         self%limGrazBac,'limGrazBac','mmol/m**3',&
         'Limiting parameter for bacteria grazing by Het',&
         default=2._rk)
    !----N---------------
    call self%get_parameter(&
         self%K_psi,'K_psi','[nd]',&
         'Strength of NH4 inhibition of NO3 uptake constant',&
         default=1.46_rk)
    call self%get_parameter(&
         self%K_nox_lim,'K_nox_lim','[mmol/m**3]',&
         'Half-sat.const.for uptake of NO3+NO2',&
         default=0.15_rk)
    call self%get_parameter(&
         self%K_nh4_lim,'K_nh4_lim','[mmol/m**3]',&
         'Half-sat.const.for uptake of NH4',&
         default=0.02_rk)
    call self%get_parameter(&
         self%K_nfix,'K_nfix','[1/d]',&
         'Max. specific rate of mitrogen fixation',&
         default=10._rk)
    !----P---------------
    call self%get_parameter(&
         self%K_po4_lim,'K_po4_lim','[mmol/m**3]',&
         'Half-sat. constant for uptake of PO4 by Phy',&
         default=0.02_rk)
    !----Si-------------
    call self%get_parameter(&
         self%K_si_lim,'K_si_lim','[mmol/m**3]',&
         'Half-sat. constant for uptake of Si by Phy',&
         default=0.02_rk)
    !----Sinking--------
    call self%get_parameter(&
         self%Wsed,'Wsed','[1/day]',&
         'Rate of sinking of detritus (POP, PON)',&
         default=5.00_rk)
    call self%get_parameter(&
         self%Wphy,'Wphy','[m/day]',&
         'Rate of sinking of Phy',&
         default=0.10_rk)
    call self%get_parameter(&
         self%Whet,'Whet','[m/day]',&
         'Rate of sinking of Het',&
         default=1.00_rk)
    !----Stoichiometric coefficients----!
    call self%get_parameter(self%r_n_p,'r_n_p','[-]','N[uM]/P[uM]',&
                            default=16.0_rk)
    call self%get_parameter(self%r_o_n,'r_o_n','[-]','O[uM]/N[uM]',&
                            default=6.625_rk)
    call self%get_parameter(self%r_c_n,'r_c_n','[-]','C[uM]/N[uM]',&
                            default=8.0_rk)
    call self%get_parameter(self%r_si_n,'r_si_n','[-]','Si[uM]/N[uM]',&
                            default=2.0_rk)
    !Register state variables
    call self%register_state_variable(&
         self%id_Phy,'Phy','mmol/m**3','Phy',&
         minimum=0.0001_rk,initial_value=0.0001_rk,&
         vertical_movement=-self%Wphy/86400._rk)
    call self%register_state_variable(&
         self%id_Het,'Het','mmol/m**3','Het',minimum=0.0_rk,&
         vertical_movement=-self%Whet/86400._rk)
    call self%register_state_variable(&
         self%id_PON,'PON','mmol/m**3','PON',minimum=0.0_rk,&
         vertical_movement=-self%Wsed/86400._rk)
    call self%register_state_variable(&
         self%id_DON,'DON','mmol/m**3','DON',minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_O2,'O2','mmol/m**3','O2',minimum=0.0_rk)
    !Register state dependencies
    call self%register_state_dependency(&
         self%id_PO4,'PO4','mmol/m**3','PO4')
    call self%register_state_dependency(&
         self%id_Si,'Si','mmol/m**3','Si')
    call self%register_state_dependency(&
         self%id_Sipart,'Sipart','mmol/m**3','Si particulate')
    call self%register_state_dependency(&
         self%id_NO3,'NO3','mmol/m**3','NO3')
    call self%register_state_dependency(&
         self%id_NH4,'NH4','mmol/m**3','NH4')
    call self%register_state_dependency(&
         self%id_NO2,'NO2','mmol/m**3','NO2')
    call self%register_state_dependency(&
         self%id_DIC,'DIC','mmol/m**3','DIC')
    call self%register_state_dependency(&
         self%id_H2S,'H2S','mmol/m**3','H2S')
    call self%register_state_dependency(self%id_Alk,&
         standard_variables%alkalinity_expressed_as_mole_equivalent)
    call self%register_state_dependency(&
         self%id_Baae,'Baae','mmol/m**3','aerobic autotrophic bacteria')
    call self%register_state_dependency(&
         self%id_Bhae,'Bhae','mmol/m**3','aerobic heterotrophic bacteria')
    call self%register_state_dependency(&
         self%id_Baan,'Baan','mmol/m**3','anaerobic aurotrophic bacteria')
    call self%register_state_dependency(&
         self%id_Bhan,'Bhan','mmol/m**3','anaerobic heterotrophic bacteria')
    !diagnostic dependency
    call self%register_dependency(&
         self%id_Hplus,'Hplus','mmol/m**3','H+ hydrogen', required=.false.)
    !Register diagnostic variables
    call self%register_diagnostic_variable(&
         self%id_DcPM_O2,'DcPM_O2','mmol/m**3',&
         'POM with O2 mineralization',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDM_O2,'DcDM_O2','mmol/m**3',&
         'DOM with O2 mineralization',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_MortHet,'MortHet','mmol/m**3','Mortality of Het',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Grazing,'Grazing','mmol/m**3','Grazing of Het',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_RespHet,'RespHet','mmol/m**3','Respiration rate of Het',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrazBhae,'GrazBhae','mmol/m**3','GrazBhae',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrazBhan,'GrazBhan','mmol/m**3','GrazBhan',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrazBaae,'GrazBaae','mmol/m**3','GrazBaae',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrazBaan,'GrazBaan','mmol/m**3','GrazBaan',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrazPhy,'GrazPhy','mmol/m**3','GrazPhy',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrazPOP,'GrazPOP','mmol/m**3','GrazPOP',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrazBact,'GrazBact','mmol/m**3','GrazBact',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_MortPhy,'MortPhy','mmol/m**3','MortPhy',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_ExcrPhy,'ExcrPhy','mmol/m**3','ExcrPhy',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimNH4,'LimNH4','mmol/m**3','LimNH4',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimN,'LimN','mmol/m**3','LimN',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrowthPhy,'GrowthPhy','mmol/m**3','GrowthPhy',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimT,'LimT','mmol/m**3','LimT',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimP,'LimP','mmol/m**3','LimP',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimNO3,'LimNO3','mmol/m**3','LimNO3',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimSi,'LimSi','mmol/m**3','LimSi',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimLight,'LimLight','mmol/m**3','LimLight',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_N_fixation,'N_fixation','mmol/m**3/d','N_fixation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_eO2mO2,'eO2mO2',&
         '1','relative oxygen saturation', &
         standard_variable=standard_variables%fractional_saturation_of_oxygen)
    call self%register_diagnostic_variable(&
         self%id_osat,'osat',&
         'mmol O_2/m^3','oxygen saturation concentration')
    !Register environmental dependencies
    call self%register_dependency(&
         self%id_par,&
         standard_variables%downwelling_photosynthetic_radiative_flux)
    call self%register_dependency(&
         self%id_temp,standard_variables%temperature)
    call self%register_dependency(&
         self%id_pres,standard_variables%pressure)
    call self%register_dependency(&
         self%id_salt,standard_variables%practical_salinity)
    call self%register_dependency(&
         self%id_windspeed,standard_variables%wind_speed)
    !Specify that are rates computed in this module are per day
    !(default: per second)
    self%dt = 86400._rk
  end subroutine initialize
  !
  !
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_bio),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    real(rk):: temp,salt,pres,Iz
    real(rk):: NH4,NO2,NO3,PO4,Phy,Het,H2S,O2,Baae,Baan,Bhae,Bhan
    real(rk):: PON,DON,Si,Sipart,Alk,Hplus
    real(rk):: LimLight,LimT,LimP,LimNO3,LimNH4,LimN,LimSi
    real(rk):: GrowthPhy,MortPhy,ExcrPhy,dAlk,N_fixation
    real(rk):: GrazPhy,GrazPOP,GrazBaae,GrazBaan,GrazBhae
    real(rk):: GrazBhan,GrazBact,Grazing,RespHet,MortHet
    real(rk):: Autolysis,DcDM_O2,DcPM_O2,Dc_OM_total
    real(rk):: OSAT
    integer :: phy_t_dependence ! select dependence on T: (1) Old; (2) for Arctic; (3) ERSEM
    !increments
    real(rk):: d_NO2,d_NO3,d_PO4,d_Si,d_DIC,d_O2,d_NH4
    real(rk):: d_Sipart,d_Phy,d_Het,d_Baae,d_Baan,d_Bhae,d_Bhan
    real(rk):: d_DON,d_PON

    ! Enter spatial loops (if any)
    _LOOP_BEGIN_
      ! Retrieve current environmental conditions.
      _GET_(self%id_par,Iz) ! local photosynthetically active radiation
      _GET_(self%id_temp,temp) ! temperature
      _GET_(self%id_salt,salt) ! temperature
      _GET_(self%id_pres,pres) ! pressure in dbar
      ! Retrieve current (local) state variable values.
      !diagnostic
      if (_AVAILABLE_(self%id_Hplus)) then
      _GET_(self%id_Hplus,Hplus)
      else
         Hplus = 1.0e8
      end if
      _GET_(self%id_Hplus,Hplus)
      !state variables
      _GET_(self%id_NO2,NO2)
      _GET_(self%id_NO3,NO3)
      _GET_(self%id_PO4,PO4)
      _GET_(self%id_Si,Si)
      _GET_(self%id_Alk,Alk)
      !gases
      _GET_(self%id_NH4,NH4)
      _GET_(self%id_H2S,H2S)
      _GET_(self%id_O2,O2)
      !solids
      _GET_(self%id_Phy,Phy)
      _GET_(self%id_Het,Het)
      _GET_(self%id_Baae,Baae)
      _GET_(self%id_Baan,Baan)
      _GET_(self%id_Bhae,Bhae)
      _GET_(self%id_Bhan,Bhan)
      _GET_(self%id_PON,PON)
      _GET_(self%id_DON,DON)
      _GET_(self%id_Sipart,Sipart)

      !PON and DON (Savchuk, Wulff,1996)
      Autolysis = self%K_PON_DON*PON
      !(CH2O)106(NH3)16H3PO4+106O2->106CO2+106H2O+16NH3+H3PO4
      DcDM_O2 = self%K_DON_ox*DON*O2/(O2+self%K_omox_o2) &
               *(1._rk+self%beta_da*yy(self%tda,temp))
      DcPM_O2 = self%K_PON_ox*PON*O2/(O2+self%K_omox_o2) &
               *(1._rk+self%beta_da*yy(self%tda,temp))

      !Phy
      !Influence of the Irradiance on photosynthesis
      LimLight = Iz/self%Iopt*exp(1._rk-Iz/self%Iopt)
      !Influence of Temperature on photosynthesis
      LimT = f_t(temp,self%phy_t_dependence) !, phy_t_dependence)
      !dependence of photosynthesis on P
      LimP = yy(self%K_po4_lim*self%r_n_p,PO4/max(Phy,1.e-10_rk))
      !dependence of photosynthesis on Si
      LimSi = yy(self%K_si_lim/self%r_si_n,Si/max(Phy,1.e-10_rk))
      !dependence of photosynthesis on NO3+NO2
      LimNO3 = yy(self%K_nox_lim,(NO3+NO2)/max(Phy,1.e-10_rk))*&
               exp(-self%K_psi*(NH4/max(Phy,1.e-10_rk)))
      !dependence of photosynthesis on NH4
      LimNH4 = yy(self%K_nh4_lim,NH4/max(Phy,1.e-10_rk))*&
               (1._rk-exp(-self%K_psi*(NH4/max(Phy,1.e-10_rk))))
      !dependence of photosynthesis on N
      LimN = min(1._rk,LimNO3+LimNH4)
      LimN = max(0.0001,LimN)
      !Grouth of Phy (gross primary production in uM N)
      GrowthPhy = self%K_phy_gro*LimLight*LimT*min(LimP,LimN,LimSi)*Phy
      !Rate of mortality of phy
      MortPhy = max(0.99_rk,(self%K_phy_mrt+(&
                0.5_rk-0.5_rk*tanh(O2-60._rk))*&
                0.45_rk+(0.5_rk-0.5_rk*tanh(O2-20._rk))*0.45_rk))*Phy
      !Excretion of phy
      ExcrPhy = self%K_phy_exc*Phy

      !Het
      !Grazing of Het on phy
      GrazPhy = self%K_het_phy_gro*Het*&
                yy(self%K_het_phy_lim,Phy/(Het+0.0001_rk))
      !Grazing of Het on detritus
      GrazPOP = self%K_het_pom_gro*Het*&
                yy(self%K_het_pom_lim,PON/(Het+0.0001_rk))
      !Grazing of Het on  bacteria
      GrazBaae = 1.0_rk*self%K_het_bac_gro*Het*&
                 yy(self%limGrazBac,Baae/(Het+0.0001_rk))
      GrazBaan = 0.5_rk*self%K_het_bac_gro*Het*&
                 yy(self%limGrazBac,Baan/(Het+0.0001_rk))
      GrazBhae = 1.0_rk*self%K_het_bac_gro*Het*&
                 yy(self%limGrazBac, Bhae/(Het+0.0001_rk))
      GrazBhan = 1.3_rk*self%K_het_bac_gro*Het*&
                 yy(self%limGrazBac, Bhan/(Het+0.0001_rk))
      GrazBact =GrazBaae+GrazBaan+GrazBhae+GrazBhan
      !Total grazing of Het
      Grazing = GrazPhy+GrazPOP+GrazBact
      !Respiration of Het
      RespHet = self%K_het_res*Het*(0.5_rk+0.5_rk*tanh(O2-20._rk))
      MortHet = (0.25_rk+(0.5_rk-0.5_rk*tanh(O2-20._rk))*0.3_rk+&
                (0.5_rk+0.4_rk*tanh(H2S-10._rk))*0.45_rk)*Het

      !Nitrogen fixation described as appearence of NH4 available for
      !phytoplankton: N2 -> NH4 :
      N_fixation = self%K_nfix*LimP*&
                   1._rk/(1._rk+((NO3+NO2+NH4)/max(PO4,0.000001_rk)*16._rk)**4._rk)*GrowthPhy

      !Summariazed OM mineralization
      Dc_OM_total = DcPM_O2+DcDM_O2

      !components of temporal derivarives calculated in this module:
      !Changes in alkalinity
      dAlk = &
             !the nutrient-H+-compensation principle.
             !Formulated by Wolf-Gladrow et al., 2007 :
             !"an increase of alkalinity by 1 mole when nitrate or
             !nitrite is the N source,
             +1._rk*GrowthPhy*(LimNO3/LimN)&
             !decrease of H+ to compensate NO3 consumption
             !and a decrease of alkalinity by 1 mole when ammonia is used"
             -1._rk*GrowthPhy*(LimNH4/LimN)+N_fixation !&
!             + Dc_OM_total
      _SET_ODE_(self%id_Alk,dAlk)
      d_NO2 = (-GrowthPhy*(LimNO3/LimN)*&
               (NO2/(0.00001_rk+NO2+NO3)))
      _SET_ODE_(self%id_NO2,d_NO2)
      d_NO3 = (-GrowthPhy*(LimNO3/LimN)*&
               ((NO3+0.00001_rk)/(0.00001_rk+NO2+NO3)))
      _SET_ODE_(self%id_NO3,d_NO3)
      d_PO4 = ((Dc_OM_total-GrowthPhy+RespHet)/self%r_n_p)
      _SET_ODE_(self%id_PO4,d_PO4)
      d_Si = ((-GrowthPhy+ExcrPhy)*self%r_si_n)
      _SET_ODE_(self%id_Si,d_Si)
      d_DON = (Autolysis-DcDM_O2+ExcrPhy+Grazing*(1._rk-self%Uz)*self%Hz)
      _SET_ODE_(self%id_DON,d_DON)
      d_DIC = ((Dc_OM_total-GrowthPhy+RespHet)*self%r_c_n)
      _SET_ODE_(self%id_DIC,d_DIC)
      !gases
      d_O2 = ((-DcDM_O2-DcPM_O2+GrowthPhy-RespHet)*self%r_o_n)
      _SET_ODE_(self%id_O2,d_O2)
      d_NH4 = (Dc_OM_total+RespHet-GrowthPhy*(LimNH4/LimN)+N_fixation)
      _SET_ODE_(self%id_NH4,d_NH4)
      !solids
      d_Sipart = ((MortPhy+GrazPhy)*self%r_si_n) !+ExcrPhy
      _SET_ODE_(self%id_Sipart,d_Sipart)
      d_Phy = (GrowthPhy-MortPhy-ExcrPhy-GrazPhy)
      _SET_ODE_(self%id_Phy,d_Phy)
      d_Het = (self%Uz*Grazing-MortHet-RespHet)
      _SET_ODE_(self%id_Het,d_Het)
      d_Baae = -GrazBaae
      _SET_ODE_(self%id_Baae,d_Baae)
      d_Baan = -GrazBaan
      _SET_ODE_(self%id_Baan,d_Baan)
      d_Bhae = -GrazBhae
      _SET_ODE_(self%id_Bhae,d_Bhae)
      d_Bhan = -GrazBhan
      _SET_ODE_(self%id_Bhan,d_Bhan)
      d_PON = (-Autolysis-DcPM_O2+MortPhy+MortHet+Grazing*&
               (1._rk-self%Uz)*(1._rk-self%Hz)-GrazPOP)
      _SET_ODE_(self%id_PON,d_PON)

      OSAT = oxygen_saturation_concentration(temp,salt)
      
      _SET_DIAGNOSTIC_(self%id_osat,OSAT)
      _SET_DIAGNOSTIC_(self%id_eO2mO2,max(0.0_rk,O2/OSAT))
      _SET_DIAGNOSTIC_(self%id_DcPM_O2,DcPM_O2)
      _SET_DIAGNOSTIC_(self%id_DcDM_O2,DcDM_O2)
      _SET_DIAGNOSTIC_(self%id_MortHet,MortHet)
      _SET_DIAGNOSTIC_(self%id_Grazing,Grazing)
      _SET_DIAGNOSTIC_(self%id_RespHet,RespHet)
      _SET_DIAGNOSTIC_(self%id_GrazBhae,GrazBhae)
      _SET_DIAGNOSTIC_(self%id_GrazBhan,GrazBhan)
      _SET_DIAGNOSTIC_(self%id_GrazBaae,GrazBaae)
      _SET_DIAGNOSTIC_(self%id_GrazBaan,GrazBaan)
      _SET_DIAGNOSTIC_(self%id_GrazPhy,GrazPhy)
      _SET_DIAGNOSTIC_(self%id_GrazPOP,GrazPOP)
      _SET_DIAGNOSTIC_(self%id_GrazBact,GrazBact)
      _SET_DIAGNOSTIC_(self%id_MortPhy,MortPhy)
      _SET_DIAGNOSTIC_(self%id_ExcrPhy,ExcrPhy)
      _SET_DIAGNOSTIC_(self%id_LimNH4,LimNH4)
      _SET_DIAGNOSTIC_(self%id_LimN,LimN)
      _SET_DIAGNOSTIC_(self%id_GrowthPhy,GrowthPhy)
      _SET_DIAGNOSTIC_(self%id_LimT,LimT)
      _SET_DIAGNOSTIC_(self%id_LimP,LimP)
      _SET_DIAGNOSTIC_(self%id_LimNO3,LimNO3)
      _SET_DIAGNOSTIC_(self%id_LimSi,LimSi)
      _SET_DIAGNOSTIC_(self%id_LimLight,LimLight)
      _SET_DIAGNOSTIC_(self%id_N_fixation,N_fixation)
    _LOOP_END_
  end subroutine do
  !
  !O2 saturation
  !
  subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
    class (type_niva_brom_bio),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_SURFACE_
    real(rk)                   :: O2, temp, salt, windspeed
    real(rk)                   :: Ox, Oa, TempT, Obe, Q_O2

    _HORIZONTAL_LOOP_BEGIN_
      _GET_(self%id_O2,O2)
      _GET_(self%id_temp,temp) ! temperature
      _GET_(self%id_salt,salt) ! salinity
      _GET_HORIZONTAL_(self%id_windspeed,windspeed)

      Ox = 1953.4_rk-128._rk*temp+3.9918_rk*temp*temp-&
           0.050091_rk*temp*temp*temp !(Wanninkoff, 1992)
      if (Ox>0._rk) then
        Oa = 0.028_rk*(windspeed**3._rk)*sqrt(400._rk/Ox)
      else
        Oa = 0._rk
      end if

      !Calculation of O2 saturation Obe according to UNESCO, 1986
      TempT = (temp+273.15_rk)/100._rk
      Obe = exp(-173.4292_rk+249.6339_rk/TempT+143.3483_rk*&
            log(TempT)-21.8492_rk*TempT+salt*(-0.033096_rk+&
            0.014259_rk*TempT-0.0017_rk*TempT*TempT)) !Osat
      Obe = Obe*1000._rk/22.4_rk !convert from ml/l into uM
      Q_O2 = windspeed*(Obe-O2) !After (Burchard et al., 2005)

      _SET_SURFACE_EXCHANGE_(self%id_O2,Q_O2)
    _HORIZONTAL_LOOP_END_
  end subroutine
  !
  !Phy temperature limiter
  !
  elemental real(rk) function f_t(temperature,phy_t_dependence ) !, phy_t_dependence)
    real(rk),intent(in):: temperature
    integer,intent(in) :: phy_t_dependence
!    case (1) ! Old
    real(rk):: bm
    real(rk):: cm
!    case (2) ! for Arctic
    real(rk):: t_0 !reference temperature
    real(rk):: temp_aug_rate !temperature augmentation rate
!    case (3) ! ERSEM
    real(rk):: q10 !Coefficient for uptake rate dependence on t
    real(rk):: t_upt_min !Low t limit for uptake rate dependence on t
    real(rk):: t_upt_max !High t limit for uptake rate dependence on t

!    select dependence on T: (1) Old; (2) for Arctic; (3) ERSEM
    select case (phy_t_dependence)
     case (1) ! Old
        bm = 0.12_rk
        cm = 1.4_rk
     case (2) ! for Arctic   !(Moore et al.,2002;Jin et al.,2008)
        t_0           = 0._rk
        temp_aug_rate = 0.0663_rk
     case (3) ! ERSEM
        q10       = 2.0_rk
        t_upt_min = 10.0_rk
        t_upt_max = 32.0_rk
    end select
!   Some others:
 !  LimT     = 0.5(1+tanh((t-tmin)/smin)) (1-0.5(1+th((t-tmax)/smax))) !Smin= 15  Smax= 15  Tmin=  10 Tmax= 35   (Deb et al., .09)
 !  LimT     = exp(self%bm*temp-self%cm))        !Dependence on Temperature (used in (Ya,So, 2011) for Arctic)  
 !  LimT     = 1./(1.+exp(10.-temp))             !Dependence on Temperature (ERGOM for cya)
 !  LimT     = 1.-temp*temp/(temp*temp +12.*12.) !Dependence on Temperature (ERGOM for dia)
 !  LimT     = 2.**((temp- 10.)/10.) -2**((temp-32.)/3.) !(ERSEM)
 !  LimT     =q10*(T-20)/10 !Q10=1.88 (Gregoire, 2000)       

    select case (phy_t_dependence)
     case (1) ! Old
      f_t = exp(bm*temperature-cm) !Influence of T photosynthesis
     case (2) ! for Arctic
      f_t = exp(temp_aug_rate*(temperature-t_0))
     case (3) ! ERSEM
      f_t = q10**((temperature-t_upt_min)/10.)-q10**((temperature-t_upt_max)/3.)
    end select

  end function f_t
  !
  !adapted from ersem
  !
  function oxygen_saturation_concentration(ETW,X1X) result(OSAT)
    real(rk),                      intent(in) :: ETW,X1X
    real(rk)                                  :: OSAT

    real(rk),parameter :: A1 = -173.4292_rk
    real(rk),parameter :: A2 = 249.6339_rk
    real(rk),parameter :: A3 = 143.3483_rk
    real(rk),parameter :: A4 = -21.8492_rk
    real(rk),parameter :: B1 = -0.033096_rk
    real(rk),parameter :: B2 = 0.014259_rk
    real(rk),parameter :: B3 = -0.0017_rk
    real(rk),parameter :: R = 8.3145_rk
    real(rk),parameter :: P = 101325_rk
    real(rk),parameter :: T = 273.15_rk

    ! volume of an ideal gas at standard temp (25C) and pressure (1 atm)
    real(rk),parameter :: VIDEAL = (R * 298.15_rk / P) *1000._rk

    real(rk)           :: ABT

    ! calc absolute temperature
    ABT = ETW + T

    ! calc theoretical oxygen saturation for temp + salinity
    ! From WEISS 1970 DEEP SEA RES 17, 721-735.
    ! units of ln(ml(STP)/l)
    OSAT = A1 + A2 * (100._rk/ABT) + A3 * log(ABT/100._rk) &
            + A4 * (ABT/100._rk) &
            + X1X * ( B1 + B2 * (ABT/100._rk) + B3 * ((ABT/100._rk)**2))

    ! convert units to ml(STP)/l then to mMol/m3
    OSAT = exp( OSAT )
    OSAT = OSAT * 1000._rk / VIDEAL
  end function
  !
  ! Original author(s): Hans Burchard, Karsten Bolding
  ! DESCRIPTION:
  ! This is a squared Michaelis-Menten type of limiter
  !
  real(rk) function yy(a,x)
    real(rk),intent(in):: a,x

    yy=x**2._rk/(a**2._rk+x**2._rk)
  end function yy
end module fabm_niva_brom_bio
