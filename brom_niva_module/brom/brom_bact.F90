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

module fabm_niva_brom_bact
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public :: type_niva_brom_bact
    !variables allocated here
    type(type_state_variable_id):: id_Baae,id_Bhae,id_Baan,id_Bhan
    !state dependencies
    type(type_state_variable_id):: id_O2,id_NH4,id_H2S
    type(type_state_variable_id):: id_PON,id_DON,id_DIC,id_PO4
    type(type_state_variable_id):: id_Alk

    type(type_diagnostic_variable_id):: id_ChemBaae,id_ChemBaan
    type(type_diagnostic_variable_id):: id_HetBhan,id_HetBhae
    type(type_diagnostic_variable_id):: id_MortBhan,id_MortBaan
    type(type_diagnostic_variable_id):: id_MortBaae,id_MortBhae

    !processes we need to calculate bact
    type(type_dependency_id):: id_Nitrif1,id_Nitrif2
    type(type_dependency_id):: id_mn_ox1,id_s0_ox,id_anammox
    type(type_dependency_id):: id_DcPM_O2,id_DcDM_O2
    type(type_dependency_id):: id_mn_rd1,id_mn_rd2,id_fe_rd
    type(type_dependency_id):: id_hs_no3,id_hs_ox
    type(type_dependency_id):: id_s2o3_ox,id_fe_ox1
    type(type_dependency_id):: id_DcPM_NOX,id_DcDM_NOX
    type(type_dependency_id):: id_DcPM_SO4,id_DcDM_SO4
    type(type_dependency_id):: id_DcDM_Mn4,id_DcPM_Mn4
    type(type_dependency_id):: id_DcPM_Fe,id_DcDM_Fe
    type(type_dependency_id):: id_DcPM_ch4

    !Model parameters
    !sinking
    real(rk):: Wbact
    !specific rates of biogeochemical processes
    !----Bacteria-!
    real(rk):: K_Baae_gro,K_Baae_mrt,K_Baae_mrt_h2s,limBaae
    real(rk):: K_Bhae_gro,K_Bhae_mrt,K_Bhae_mrt_h2s,limBhae
    real(rk):: K_Baan_gro,K_Baan_mrt,limBaan
    real(rk):: K_Bhan_gro,K_Bhan_mrt,K_Bhan_mrt_o2,limBhan
    !---- Stoichiometric coefficients ----!
    real(rk):: r_n_p,r_c_n
  contains
    procedure :: initialize
    procedure :: do
  end type
contains
  !
  !
  !
  subroutine initialize(self,configunit)
    class (type_niva_brom_bact), intent(inout), target :: self
    integer,                      intent(in)            :: configunit

    !-----Model parameters------
    !Sinking
    call self%get_parameter(&
         self%Wbact,'Wbact','[1/day]',&
         'Rate of sinking of bacteria (Bhae,Baae,Bhan,Baan)',&
         default=0.4_rk)
    !Specific rates of biogeochemical processes
    !----Bacteria-!
    call self%get_parameter(&
         self%K_Baae_gro, 'K_Baae_gro', '[1/day]',&
         'Baae maximum specific growth rate',&
         default=0.019_rk)
    call self%get_parameter(&
         self%K_Baae_mrt, 'K_Baae_mrt', '[1/day]',&
         'Baae specific rate of mortality',&
         default=0.005_rk)
    call self%get_parameter(&
         self%K_Baae_mrt_h2s, 'K_Baae_mrt_h2s', '[1/day]',&
         'Baae increased specific rate of mortality due to H2S',&
         default=0.899_rk)
    call self%get_parameter(&
         self%limBaae, 'limBaae', '[1/day]',&
         'Limiting parameter for nutrient consumprion by Baae',&
         default=2.0_rk)
    call self%get_parameter(&
         self%K_Bhae_gro, 'K_Bhae_gro', '[1/day]',&
         'Bhae maximum specific growth rate',&
         default=0.5_rk)
    call self%get_parameter(&
         self%K_Bhae_mrt, 'K_Bhae_mrt', '[1/day]',&
         'Bhae specific rate of mortality',&
         default=0.01_rk)
    call self%get_parameter(&
         self%K_Bhae_mrt_h2s, 'K_Bhae_mrt_h2s', '[1/day]',&
         'Bhae increased specific rate of mortality due to H2S',&
         default=0.799_rk)
    call self%get_parameter(&
         self%limBhae, 'limBhae', '[1/day]',&
         'Limiting parameter for OM consumprion by Bhae',&
         default=5.0_rk)
    call self%get_parameter(&
         self%K_Baan_gro, 'K_Baan_gro', '[1/day]',&
         'Baan maximum specific growth rate',&
         default=0.012_rk)
    call self%get_parameter(&
         self%K_Baan_mrt, 'K_Baan_mrt', '[1/day]',&
         'Baan specific rate of mortality',&
         default=0.012_rk)
    call self%get_parameter(&
         self%limBaan, 'limBaan', '[1/day]',&
         'Limiting parameter for nutrient consumprion by Baan',&
         default=2.0_rk)
    call self%get_parameter(&
         self%K_Bhan_gro, 'K_Bhan_gro', '[1/day]',&
         'Bhan maximum specific growth rate',&
         default=0.2_rk)
    call self%get_parameter(&
         self%K_Bhan_mrt, 'K_Bhan_mrt', '[1/day]',&
         'Bhan specific rate of mortality',&
         default=0.007_rk)
    call self%get_parameter(&
         self%K_Bhan_mrt_o2, 'K_Bhan_mrt_o2', '[1/day]',&
         'Bhan increased specific rate of mortality due to O2',&
         default=0.899_rk)
    call self%get_parameter(&
         self%limBhan, 'limBhan', '[1/day]',&
         'Limiting parameter for OM consumprion by Bhan',&
         default=2.0_rk)
    !----Stoichiometric coefficients----!
    call self%get_parameter(&
         self%r_n_p,   'r_n_p',  '[-]',&
         'N[uM]/P[uM]',&
         default=16.0_rk)
    call self%get_parameter(&
         self%r_c_n,   'r_c_n',  '[-]',&
         'C[uM]/N[uM]',&
         default=8.0_rk)

    !Register state variables
    call self%register_state_variable(&
         self%id_Baae, 'Baae', 'mmol/m**3','Aerobic Autotrophic Bacteria',&
         minimum=0.0_rk,vertical_movement=-self%Wbact/86400._rk)
    call self%register_state_variable(&
         self%id_Bhae, 'Bhae', 'mmol/m**3','Aerobic Heterotrophic&
         Bacteria',&
         minimum=0.0_rk,vertical_movement=-self%Wbact/86400._rk)
    call self%register_state_variable(&
         self%id_Baan, 'Baan', 'mmol/m**3','Anaerobic Autotrophic&
         Bacteria',&
         minimum=0.0_rk,vertical_movement=-self%Wbact/86400._rk)
    call self%register_state_variable(&
         self%id_Bhan, 'Bhan', 'mmol/m**3','Anaerobic Heterotrophic&
         Bacteria',&
         minimum=0.0_rk,vertical_movement=-self%Wbact/86400._rk)

    !Register dependencies
    call self%register_dependency(&
         self%id_Nitrif1,'Nitrif1','mmol/m**3',&
         'Nitrification 1 stage')
    call self%register_dependency(&
         self%id_Nitrif2,'Nitrif2','mmol/m**3',&
         'Nitrification 2 stage')
    call self%register_dependency&
         (self%id_anammox ,'anammox','mmol/m**3',&
         'Anammox')
    call self%register_dependency&
         (self%id_mn_ox1,'mn_ox1','mmol/m**3',&
         'Mn(II) with O2 oxidation')
    call self%register_dependency&
         (self% id_fe_ox1,'fe_ox1','mmol/m**3',&
         'Fe(II) with O2 oxidation')
    call self%register_dependency&
         (self%id_s0_ox,'s0_ox','mmol/m**3',&
         'S0 with O2 oxidation')
    call self%register_dependency&
         (self%id_DcPM_O2,'DcPM_O2','mmol/m**3',&
         'POM with O2 mineralization')
    call self%register_dependency&
         (self%id_DcDM_O2,'DcDM_O2','mmol/m**3',&
         'DOM with O2 mineralization')
    call self%register_dependency&
         (self%id_mn_rd1,'mn_rd1','mmol/m**3',&
         'Mn(IV) with H2S reduction')
    call self%register_dependency&
         (self%id_mn_rd2,'mn_rd2','mmol/m**3',&
         'Mn(III) with H2S reduction')
    call self%register_dependency&
         (self%id_fe_rd,'fe_rd','mmol/m**3',&
         'Fe (III) with H2S reduction')
    call self%register_dependency&
         (self%id_hs_ox,'hs_ox','mmol/m**3',&
         'H2S with O2 oxidation')
    call self%register_dependency&
         (self%id_hs_no3,'hs_no3','mmol/m**3',&
         'H2S with NO3 oxidation')
    call self%register_dependency&
         (self%id_s2o3_ox,'s2o3_ox','mmol/m**3',&
         'Specific rate of oxidation of S2O3 with O2')
    call self%register_dependency&
         (self%id_DcPM_NOX,'DcPM_NOX','mmol/m**3',&
         'POM denitrification (1+2 stage)')
    call self%register_dependency&
         (self%id_DcDM_NOX,'DcDM_NOX','mmol/m**3',&
         'DOM denitrification (1+2 stage)')
    call self%register_dependency&
         (self%id_DcPM_SO4,'DcPM_SO4','mmol/m**3',&
         'POM sulfatereduction (1+2 stage)')
    call self%register_dependency&
         (self%id_DcDM_SO4,'DcDM_SO4','mmol/m**3',&
         'DOM sulfatereduction (1+2 stage)')
    call self%register_dependency&
         (self%id_DcPM_Mn4,'DcPM_Mn4','mmol/m**3',&
         'POM with Mn(IV) mineralization')
    call self%register_dependency&
         (self%id_DcDM_Mn4,'DcDM_Mn4','mmol/m**3',&
         'DOM with Mn(IV) mineralization')
    call self%register_dependency&
         (self%id_DcPM_Fe,'DcPM_Fe','mmol/m**3',&
         'POM with Fe(III) mineralization')
    call self%register_dependency&
         (self%id_DcDM_Fe,'DcDM_Fe','mmol/m**3',&
         'DOM with Fe(III) mineralization')
    call self%register_dependency&
         (self%id_DcPM_ch4,'DcPM_ch4','mmol/m**3',&
         'CH4 production from PON and DON')

    !Register state dependencies
    call self%register_state_dependency(&
         self%id_NH4,'NH4','mmol/m**3',&
         'ammonium')
    call self%register_state_dependency(&
         self%id_PON,'PON','mmol/m**3',&
         'particulate organic nitrogen')
    call self%register_state_dependency(&
         self%id_DON,'DON','mmol/m**3',&
         'dissolved organic nitrogen')
    call self%register_state_dependency(&
         self%id_po4,'PO4','mmol/m**3',&
         'phosphate',required=.false.)
    call self%register_state_dependency(&
         self%id_O2, 'O2', 'mmol/m**3',&
         'dissolved oxygen')
    call self%register_state_dependency(&
         self%id_H2S,'H2S','mmol/m**3',&
         'hydrogen sulfide')
    call self%register_state_dependency(&
         self%id_DIC,'DIC','mmol/m**3',&
         'total dissolved inorganic carbon',required=.false.)
    call self%register_state_dependency(self%id_Alk,&
         standard_variables%alkalinity_expressed_as_mole_equivalent)

    !Register diagnostic variables
    call self%register_diagnostic_variable(&
         self%id_ChemBaae,'ChemBaae','mmol/m**3',&
         'Growth of  Aerobic Autotrophic Bacteria',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_HetBhae,'HetBhae','mmol/m**3',&
         'Growth of  Aerobic Heterotrophic Bacteria',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_HetBhan,'HetBhan','mmol/m**3',&
         'Growth of Anaerobic Heterotrophic Bacteria ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_ChemBaan,'ChemBaan','mmol/m**3',&
         'Growth of Anaerobic Autotrophic Bacteria ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_MortBhan,'MortBhan','mmol/m**3',&
         'Mortality of Anaerobic Heterotrophic Bacteria ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_MortBaan,'MortBaan','mmol/m**3',&
         'Mortality of Anaerobic Autotrophic Bacteria ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_MortBaae,'MortBaae','mmol/m**3',&
         'Mortality of Aerobic Autotrophic Bacteria ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_MortBhae,'MortBhae','mmol/m**3',&
         'Mortality of Aerobic Heterotrophic Bacteria ',&
         output=output_time_step_integrated)

    !Specify that are rates computed in this module are per day
    !(default: per second)
    self%dt = 86400._rk
  end subroutine initialize
  !
  !
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_bact),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    !state variables
    real(rk):: Baae,Bhae,Baan,Bhan

    !state dependincies
    real(rk):: O2,NH4,H2S
    real(rk):: PON,DON,DIC,PO4

    !diagnostic variables
    real(rk):: ChemBaae,ChemBaan
    real(rk):: HetBhan,HetBhae
    real(rk):: MortBhan,MortBaan
    real(rk):: MortBaae,MortBhae

    !diagnostic dependencies
    real(rk):: Nitrif1,Nitrif2
    real(rk):: mn_ox1,s0_ox,anammox
    real(rk):: DcPM_O2,DcDM_O2
    real(rk):: mn_rd1,mn_rd2,fe_rd,hs_ox,hs_no3
    real(rk):: s2o3_ox,fe_ox1
    real(rk):: DcPM_NOX,DcDM_NOX
    real(rk):: DcPM_SO4,DcDM_SO4
    real(rk):: DcDM_Mn4,DcPM_Mn4
    real(rk):: DcPM_Fe,DcDM_Fe
    real(rk):: DcPM_ch4

    !increments
    real(rk):: d_Baae,d_Baan,d_Bhae,d_Bhan,d_Alk

    _LOOP_BEGIN_
      _GET_(self%id_Baae,Baae)
      _GET_(self%id_Bhae,Bhae)
      _GET_(self%id_Baan,Baan)
      _GET_(self%id_Bhan,Bhan)

      _GET_(self%id_Nitrif1,Nitrif1)
      _GET_(self%id_Nitrif2,Nitrif2)
      _GET_(self%id_mn_ox1,mn_ox1)
      _GET_(self%id_fe_ox1,fe_ox1)
      _GET_(self%id_s0_ox,s0_ox)
      _GET_(self%id_anammox,anammox)
      _GET_(self%id_DcPM_O2,DcPM_O2)
      _GET_(self%id_DcDM_O2,DcDM_O2)
      _GET_(self%id_mn_rd1,mn_rd1)
      _GET_(self%id_mn_rd2,mn_rd2)
      _GET_(self%id_fe_rd,fe_rd)
      _GET_(self%id_hs_ox,hs_ox)
      _GET_(self%id_hs_no3,hs_no3)
      _GET_(self%id_s2o3_ox,s2o3_ox)
      _GET_(self%id_DcPM_NOX,DcPM_NOX)
      _GET_(self%id_DcDM_NOX,DcDM_NOX)
      _GET_(self%id_DcPM_SO4,DcPM_SO4)
      _GET_(self%id_DcDM_SO4,DcDM_SO4)
      _GET_(self%id_DcDM_Mn4,DcDM_Mn4)
      _GET_(self%id_DcPM_Mn4,DcPM_Mn4)
      _GET_(self%id_DcPM_Fe,DcPM_Fe)
      _GET_(self%id_DcDM_Fe,DcDM_Fe)
      _GET_(self%id_DcPM_ch4,DcPM_ch4)

      _GET_(self%id_NH4,NH4)
      _GET_(self%id_PO4,PO4)
      _GET_(self%id_O2,O2)
      _GET_(self%id_PON,PON)
      _GET_(self%id_DON,DON)
      _GET_(self%id_DIC,DIC)
      _GET_(self%id_H2S,H2S)

      !Bacteria
      !OXIC CONDITIONS
      !aerobic autotrophs
      ChemBaae = (Nitrif1+Nitrif2+mn_ox1+fe_ox1+s2o3_ox+s0_ox+Anammox)&
                *self%K_Baae_gro*Baae*min(yy(self%limBaae,NH4&
                /(Baae+0.0001_rk)),yy(self%limBaae,PO4/(Baae+0.0001_rk)))
      MortBaae = (self%K_Baae_mrt+self%K_Baae_mrt_h2s&
                *(0.5_rk*(1._rk-tanh(1._rk-H2S))))*Baae*Baae
      !aerobic heterotroph
      HetBhae = (DcPM_O2+DcDM_O2)&
               *self%K_Bhae_gro*Bhae*yy(self%limBhae,DON/(Bhae+0.0001_rk))
      MortBhae = (self%K_Bhae_mrt+ self%K_Bhae_mrt_h2s&
               *(0.5_rk*(1._rk-tanh(1._rk-H2S))))*Bhae
      !ANOXIC CONDITIONS
      !anaerobic autotrophs
      ChemBaan = (mn_rd1+mn_rd2+fe_rd+hs_ox+hs_no3)&
                *self%K_Baan_gro*Baan*min(yy(self%limBaan,NH4&
                /(Baan+0.0001_rk)),yy(self%limBaan,PO4/(Baan+0.0001_rk)))
      MortBaan = self%K_Baan_mrt*Baan*Baan
      !anaerobic heterotroph
      HetBhan = (DcPM_NOX+DcDM_NOX+DcPM_SO4+DcDM_SO4+DcDM_Mn4+DcPM_Mn4 &
               +DcPM_Fe+DcDM_Fe+DcPM_ch4)&
               *self%K_Bhan_gro*Bhan*yy(self%limBhan,DON/(Bhan+0.0001_rk))
      MortBhan = (self%K_Bhan_mrt+ self%K_Bhan_mrt_o2&
                *(0.5_rk+0.5_rk*(tanh(1._rk-O2))))*Bhan

      !Alkalinity changes due to redox reactions:
      d_Alk = -ChemBaae-ChemBaan !+/- NH3
      _SET_ODE_(self%id_Alk,d_Alk)
      !Bacteria
      d_Baae = ChemBaae-MortBaae
      _SET_ODE_(self%id_Baae,d_Baae)
      d_Baan = ChemBaan-MortBaan
      _SET_ODE_(self%id_Baan,d_Baan)
      d_Bhae = HetBhae-MortBhae
      _SET_ODE_(self%id_Bhae,d_Bhae)
      d_Bhan = HetBhan-MortBhan
      _SET_ODE_(self%id_Bhan,d_Bhan)

      _SET_DIAGNOSTIC_(self%id_ChemBaae,ChemBaae)
      _SET_DIAGNOSTIC_(self%id_HetBhae,HetBhae)
      _SET_DIAGNOSTIC_(self%id_ChemBaan,ChemBaan)
      _SET_DIAGNOSTIC_(self%id_HetBhan,HetBhan)
      _SET_DIAGNOSTIC_(self%id_MortBaan,MortBaan)
      _SET_DIAGNOSTIC_(self%id_MortBhan,MortBhan)
      _SET_DIAGNOSTIC_(self%id_MortBaae,MortBaae)
      _SET_DIAGNOSTIC_(self%id_MortBhae,MortBhae)

      _SET_ODE_(self%id_NH4,-ChemBaae-ChemBaan)
      _SET_ODE_(self%id_DIC,(-ChemBaae-ChemBaan)*self%r_c_n)
      _SET_ODE_(self%id_PO4,(-ChemBaae-ChemBaan)/self%r_n_p)
      _SET_ODE_(self%id_DON,-HetBhae-HetBhan)
      _SET_ODE_(self%id_PON,MortBaae+MortBaan+MortBhae+MortBhan)
    _LOOP_END_
  end subroutine do
  !
  !Original author(s): Hans Burchard, Karsten Bolding
  !This is a squared Michaelis-Menten type of limiter
  !
  real(rk) function yy(a,x)
    real(rk),intent(in):: a,x

    yy=x**2/(a**2+x**2)
  end function yy
end module fabm_niva_brom_bact
