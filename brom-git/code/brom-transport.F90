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

    module brom_transport

    use fabm
    use fabm_config
    use fabm_types, only: attribute_length, rk
    use fabm_standard_variables
    use io_netcdf
    use io_ascii
    use mtridiagonal, only: init_tridiagonal,clean_tridiagonal
    use ids         !Provides access to variable indices id_O2 etc.


    implicit none
    private
    public init_brom_transport, do_brom_transport, clear_brom_transport


    !FABM model with all data and procedures related to biogeochemistry
    type (type_model) :: model
    type (type_horizontal_standard_variable),parameter :: &
        ice_thickness = type_horizontal_standard_variable(name=' ice_thickness',units='m')
    type (type_bulk_standard_variable),parameter :: &
        volume_of_cell = type_bulk_standard_variable(name=' volume_of_cell')


    !Solution parameters to be read from brom.yaml (see brom.yaml for details)
    integer   :: i_min, i_water, i_max     !x-axis related 
    integer   :: k_min, k_wat_bbl, k_points_below_water !z-axis related
    integer   :: max_hz_inc !, max_porosity ! URMIA: maximum increase of a layer before spliting, maximum porosity
    real(rk)  :: dt, water_layer_thickness
    integer   :: i_day, year, days_in_yr, freq_turb, freq_sed, last_day   !time related
    integer   :: diff_method, kz_bbl_type, bioturb_across_SWI  !vertical diffusivity related
    integer   :: h_adv, h_relax, h_turb  !horizontal transport  (advection and relaxation) switches
    real(rk)  :: K_O2s
    integer   :: input_type, use_Eair, use_hice, port_initial_state, ncoutfile_type !I/O related
    character(len=64) :: icfile_name, outfile_name, ncoutfile_name
    character :: hmix_file
    integer   :: k_bbl_sed, k_max          !z-axis related
    integer   :: k_air_wat                !z-axis related Urmia
    integer   :: par_max                   !no. BROM variables

    !Forcings to be provided to FABM: These must have the POINTER attribute
    real(rk), pointer, dimension(:)            :: pco2atm, windspeed, hice
    real(rk), pointer, dimension(:,:)          :: surf_flux, bott_flux, bott_source, Izt, pressure
    real(rk), pointer, dimension(:,:,:)        :: t, s, u_x
    real(rk), pointer, dimension(:,:,:)        :: vv, dvv, cc, cc_out, dcc, dcc_R, wbio

    !Surface and bottom forcings, used within brom-transport only
    real(rk), allocatable, dimension(:)        :: Eair
    real(rk), allocatable, dimension(:,:,:)    :: cc_top, cc_bottom

    !Horizontal mixing forcings, used within brom-transport only
    real(rk), allocatable, dimension(:,:,:)    :: hmix_rate
    real(rk), allocatable, dimension(:,:,:,:)  :: cc_hmix

    !Grid parameters and forcings for water column only
    real(rk), allocatable, dimension(:)        :: z_w, dz_w, hz_w
    real(rk), allocatable, dimension(:,:,:)    :: t_w, s_w, kz_w, hmix_rate_w, u_x_w
    real(rk), allocatable, dimension(:,:,:,:)  :: cc_hmix_w

    !Grid parameters and forcings for full column including water and sediments
    real(rk), allocatable, dimension(:)        :: z, dz, hz, x, dx, hx
    real(rk), allocatable, dimension(:,:)      :: bc_top, bc_bottom, kz_bio, alpha, phi, phi1, phi_inv, tortuosity, kztCFL
    real(rk), allocatable, dimension(:,:)      :: wCFL, w_b, u_b, wat_content
    real(rk), allocatable, dimension(:,:,:)    :: kz, kzti, fick, fick_per_day, sink, sink_per_day, bcpar_top, bcpar_bottom
    real(rk), allocatable, dimension(:,:,:)    :: kz_mol, pF1, pF2, wti, pWC
    integer, allocatable, dimension(:)         :: is_solid, k_wat, k_sed, k_sed1
    integer, allocatable, dimension(:,:)       :: bctype_top, bctype_bottom, hmixtype
    character(len=attribute_length), allocatable, dimension(:)    :: par_name

    !Parameters for the Urmia Lake box model
    real(rk), pointer, dimension(:,:)            :: volume, area, dvolume
    real(rk), pointer, dimension(:)              :: zeta, dzeta
    real(rk), pointer, dimension(:)              :: rv_O2, rv_NO2, rv_NO3, rv_DON, rv_PO4, rv_Si, rv_Na, rv_K, rv_Mg, rv_Ca
    real(rk), pointer, dimension(:)              :: rv_Cl, rv_SO4, rv_DIC, rv_NH4, rv_Phy, rv_Br, rv_Alk
    real(rk), pointer, dimension(:)              :: rv_DOP, rv_CO3, rv_Chla, rv_TOC, rv_TOM
    !Constant forcings that can be read as parameters from brom.yaml
    real(rk) :: wind_speed, pco2_atm, mu0_musw, dphidz_SWI , dx_adv, rho_def
    real(rk), allocatable, dimension(:)        :: rho

    !Counters
    integer   :: i, k, ip, ip_sol, ip_par, i_dummy


    contains

!=======================================================================================================================
    subroutine init_brom_transport()

    !Initialises the offline vertical transport model BROM-transport

    use ids         !Provides access to variable indices id_O2 etc

    implicit none


    !Reading brom.yaml
    call init_common()

    !Get grid and numerical solution parameters from from brom.yaml
    dt = get_brom_par("dt")
    freq_turb = get_brom_par("freq_turb")
    freq_sed  = get_brom_par("freq_sed ")
    last_day = get_brom_par("last_day")
    water_layer_thickness = get_brom_par("water_layer_thickness")
    k_min = get_brom_par("k_min")
    k_wat_bbl = get_brom_par("k_wat_bbl")
    k_points_below_water = get_brom_par("k_points_below_water")
    i_min = get_brom_par("i_min") 
    i_water = get_brom_par("i_water") !Number of i for water column
    i_max = get_brom_par("i_max")     !number of columns
    max_hz_inc = get_brom_par("max_hz_inc")   ! URMIA: maximum increase of a layer before spliting 
    !max_porosity = get_brom_par("max_porosity")   ! URMIA: maximum porosity to stop sinking
    year = get_brom_par("year")
    days_in_yr = get_brom_par("days_in_yr")
    diff_method = get_brom_par("diff_method")
    bioturb_across_SWI = get_brom_par("bioturb_across_SWI")
    input_type = get_brom_par("input_type")
    use_Eair = get_brom_par("use_Eair")
    use_hice = get_brom_par("use_hice")
    port_initial_state = get_brom_par("port_initial_state")
    icfile_name = get_brom_name("icfile_name")
    outfile_name = get_brom_name("outfile_name")
    ncoutfile_name = get_brom_name("ncoutfile_name")
    ncoutfile_type = get_brom_par("ncoutfile_type")    
    K_O2s = get_brom_par("K_O2s")
    h_adv =  get_brom_par("h_adv")
    h_relax =  get_brom_par("h_relax")
    h_turb =  get_brom_par("h_turb")
    dx_adv = get_brom_par("dx_adv")
    !Initialize FABM model from fabm.yaml
    call fabm_create_model_from_yaml_file(model)
    par_max = size(model%state_variables)

    !Allocate biological variables now that par_max is known
    allocate(surf_flux(i_max,par_max))     !surface flux (tracer unit * m/s, positive for tracer entering column)
    allocate(bott_flux(i_max,par_max))     !bottom flux (tracer unit * m/s, positive for tracer entering column)
    allocate(bott_source(i_max,k_max))     !surface flux (tracer unit * m/s, positive for tracer entering column)
    allocate(bc_top(i_max,par_max))
    allocate(bc_bottom(i_max,par_max))
    allocate(bctype_top(i_max,par_max))
    allocate(bctype_bottom(i_max,par_max))
    allocate(bcpar_top(i_max,par_max,3))
    allocate(bcpar_bottom(i_max,par_max,3))
    allocate(par_name(par_max))
    allocate(pco2atm(i_max))
    allocate(windspeed(i_max))
    allocate(cc_top(i_max,par_max,days_in_yr))
    allocate(cc_bottom(i_max,par_max,days_in_yr))
    allocate(is_solid(par_max))
    allocate(hmixtype(i_max,par_max))
    allocate(rho(par_max))


    !Retrieve the parameter names from the model structure
    do ip=1,par_max
        par_name(ip) = model%state_variables(ip)%name
    end do
    !Make the named parameter indices id_O2 etc.
    call get_ids(par_name)


    !Get boudary condition parameters from brom.yaml:
    !bctype = 0, 1, 2, 3 for no flux (default), Dirichlet constant, Dirichlet sinusoid, and Dirichlet netcdf input respectively
    do ip=1,par_max
        bctype_top(i_water,ip) = get_brom_par('bctype_top_' // trim(par_name(ip)),0.0_rk)
        if (bctype_top(i_water,ip).eq.1) then
            bc_top(i_water,ip) = get_brom_par('bc_top_' // trim(par_name(ip)))
            write(*,*) "Constant Dirichlet upper boundary condition for " // trim(par_name(ip))
            write(*,'(a, es10.3)') " = ", bc_top(i_water,ip)
        else if (bctype_top(i_water,ip).eq.2) then     !Model: bc_top = a1top + a2top*sin(omega*(julianday-a3top))
            bcpar_top(i_water,ip,1) = get_brom_par('a1top_' // trim(par_name(ip)))
            bcpar_top(i_water,ip,2) = get_brom_par('a2top_' // trim(par_name(ip)))
            bcpar_top(i_water,ip,3) = get_brom_par('a3top_' // trim(par_name(ip)))
            write(*,*) "Sinusoidal Dirichlet upper boundary condition for " // trim(par_name(ip))
            write(*,'(a, es10.3, a, es10.3, a, es10.3, a)') " = ", bcpar_top(i_water,ip,1), " + ", &
                  bcpar_top(i_water,ip,2), "*sin(omega*(julianday -", bcpar_top(i_water,ip,3), "))"
        else if (bctype_top(i_water,ip).eq.3) then     !Read from netcdf
            write(*,*) "NetCDF specified Dirichlet upper boundary condition for " // trim(par_name(ip))
        else if (bctype_top(i_water,ip).eq.4) then     !Read from ascii or calc from calintiiy
            write(*,*) "Upper boundary condition from NODC " // trim(par_name(ip))
        end if        
        bctype_bottom(i_water,ip) = get_brom_par('bctype_bottom_' // trim(par_name(ip)),0.0_rk)
        if (bctype_bottom(i_water,ip).eq.1) then
            bc_bottom(i_water,ip) = get_brom_par('bc_bottom_' // trim(par_name(ip)))
            write(*,*) "Constant Dirichlet lower boundary condition for " // trim(par_name(ip))
            write(*,'(a, es10.3)') " = ", bc_bottom(i_water,ip)
        else if (bctype_bottom(i_water,ip).eq.2) then  !Model: bc_bottom = a1bottom + a2bottom*sin(omega*julianday-a3bottom))
            bcpar_bottom(i_water,ip,1) = get_brom_par('a1bottom_' // trim(par_name(ip)))
            bcpar_bottom(i_water,ip,2) = get_brom_par('a2bottom_' // trim(par_name(ip)))
            bcpar_bottom(i_water,ip,3) = get_brom_par('a3bottom_' // trim(par_name(ip)))
            write(*,*) "Sinusoidal Dirichlet lower boundary condition for " // trim(par_name(ip))
            write(*,'(a, es10.3, a, es10.3, a, es10.3, a)') " = ", bcpar_bottom(i_water,ip,1), " + ", &
            bcpar_bottom(i_water,ip,2), "*sin(omega*(julianday -", bcpar_bottom(i_water,ip,3), "))"
        else if (bctype_bottom(i_water,ip).eq.3) then  !Read from netcdf
            write(*,*) "NetCDF specified Dirichlet lower boundary condition for " // trim(par_name(ip))
        end if
        bctype_top(:,ip) = bctype_top(i_water,ip)
        bc_top(:,ip) = bc_top(i_water,ip)
        bctype_bottom(:,ip) = bctype_bottom(i_water,ip)
        bc_bottom(:,ip) = bc_bottom(i_water,ip)
        do k=1,3
          bcpar_top(:,ip,k) = bcpar_top(i_water,ip,k)
          bcpar_bottom(:,ip,k) = bcpar_bottom(i_water,ip,k)
        enddo
    end do
    
    write(*,*) "All other boundary conditions use surface and bottom fluxes from FABM"

    !Input forcing data
    if (input_type.eq.0) then !Input sinusoidal seasonal changes (hypothetical)
        call input_primitive_physics(z_w, dz_w, hz_w, k_wat_bbl, water_layer_thickness, t_w, s_w, kz_w, i_water, i_max, days_in_yr)
        allocate(Eair(days_in_yr))
        allocate(hice(days_in_yr))
        Eair = 0.0_rk
        hice = 0.0_rk
        write(*,*) "Done sinusoidal input"
    end if
    if (input_type.eq.1) then !Input physics from ascii
        call input_ascii_physics(z_w, dz_w, hz_w, k_wat_bbl, water_layer_thickness, t_w, s_w, kz_w, i_water, i_max, days_in_yr)
        allocate(Eair(days_in_yr))
        allocate(hice(days_in_yr))
        allocate(cc_hmix_w(i_max,par_max,k_wat_bbl,days_in_yr))
        allocate(hmix_rate_w(i_max,k_wat_bbl,days_in_yr))
        cc_hmix_w = 0.0_rk
        hmix_rate_w = 0.0_rk
        Eair = 0.0_rk
        hice = 0.0_rk
        write(*,*) "Done ascii input"
    end if
    if (input_type.eq.2) then !Input water column physics from netcdf
        call input_netcdf(z_w, dz_w, hz_w, t_w, s_w, kz_w, hmix_rate_w, Eair, use_Eair, hice, use_hice, year, i_water, i_max, &
        days_in_yr, k_wat_bbl, par_name, par_max, bctype_top, bctype_bottom, cc_top, cc_bottom, hmixtype, cc_hmix_w, h_adv, u_x_w)        
        write(*,*) "Done netcdf input"
        !Note: This uses the netCDF file to set z_w = layer midpoints, dz_w = increments between layer midpoints, hz_w = layer thicknesses
    end if

        !Determine total number of vertical grid points (layers) now that k_wat_bbl is determined
    k_max = k_wat_bbl + k_points_below_water
      if (k_points_below_water==0) then  !This is to "unlock" BBL and sediments for a "classical" water column model
!        k_max=k_wat_bbl-1
      endif
    !Allocate full grid variables now that k_max is knownk
    allocate(z(k_max))
    allocate(dz(k_max))
    allocate(hz(k_max))
    allocate(x(i_max))
    allocate(dx(i_max))
    allocate(hx(i_max))
    allocate(t(i_max,k_max,days_in_yr))
    allocate(s(i_max,k_max,days_in_yr))
    allocate(kz(i_max,k_max+1,days_in_yr))
    allocate(hmix_rate(i_max,k_max,days_in_yr))
    allocate(u_x(i_max,k_max,days_in_yr))
    allocate(cc_hmix(i_max,par_max,k_max,days_in_yr))
    allocate(kz_mol(i_max,k_max+1,par_max))
    allocate(kz_bio(i_max,k_max+1))
    allocate(pF1(i_max,k_max,par_max))
    allocate(pF2(i_max,k_max+1,par_max))
    allocate(pWC(i_max,k_max+1,par_max))
    allocate(alpha(i_max,k_max))
    allocate(phi(i_max,k_max))
    allocate(wat_content(i_max,k_max))
    allocate(phi1(i_max,k_max+1))
    allocate(phi_inv(i_max,k_max))
    allocate(tortuosity(i_max,k_max+1))
    allocate(w_b(i_max,k_max+1))
    allocate(u_b(i_max,k_max+1))
    allocate(wti(i_max,k_max+1,par_max))
    allocate(cc(i_max,k_max,par_max))
    allocate(cc_out(i_max,k_max,par_max))
    allocate(dcc(i_max,k_max,par_max))
    allocate(dcc_R(i_max,k_max,par_max))
    allocate(fick(i_max,k_max+1,par_max))
    allocate(fick_per_day(i_max,k_max+1,par_max))
    allocate(wbio(i_max,k_max,par_max))    !vertical velocity (m/s, negative for sinking)
    allocate(sink(i_max,k_max+1,par_max))  !sinking flux (mmol/m2/s, positive downward)
    allocate(sink_per_day(i_max,k_max+1,par_max))
    allocate(vv(i_max,k_max,3))
    allocate(dvv(i_max,k_max,3))
    allocate(volume(i_max+2,k_max+1))
    allocate(dvolume(i_max+2,k_max+1))
    allocate(area(i_max+2,k_max+1))
    allocate(Izt(i_max,k_max))
    allocate(pressure(i_max,k_max))
    allocate(kzti(i_max,k_max+1,par_max))
    allocate(kztCFL(k_max-1,par_max))
    allocate(wCFL(k_max-1,par_max))
! Urmia
    allocate(rv_O2(10000))
    allocate(rv_NO2(10000))
    allocate(rv_NO3(10000))
    allocate(rv_DON(10000))
    allocate(rv_PO4(10000))
    allocate(rv_Si(10000))
    allocate(rv_Na(10000))
    allocate(rv_K(10000))
    allocate(rv_Mg(10000))
    allocate(rv_Ca(10000))
    allocate(rv_Cl(10000))
    allocate(rv_SO4(10000))
    allocate(rv_DIC(10000))
    allocate(rv_NH4(10000))
    allocate(rv_Phy(10000))
    allocate(rv_Br(10000))
    allocate(rv_Alk(10000))
    allocate(rv_DOP(10000))
    allocate(rv_CO3(10000))
    allocate(rv_Chla(10000))
    allocate(rv_TOC(10000))
    allocate(rv_TOM(10000))
    
    allocate(zeta(10000))
    allocate(dzeta(10000))

    if (k_points_below_water==0) then  !This is to "unlock" BBL and sediments for a "classical" water column model
!        k_max=k_wat_bbl-1
        z=z_w
        dz=dz_w
        hz=hz_w
        k_bbl_sed=k_wat_bbl !needed for Irradiance calculations
    else
        !Construct the full vertical grid
        call make_vert_grid(z, dz, hz, z_w, dz_w, hz_w, k_wat_bbl, k_max, k_bbl_sed)
        write(*,*) "Made vertical grid"
        allocate(k_wat(k_bbl_sed))
        allocate(k_sed(k_max-k_bbl_sed))
        allocate(k_sed1(k_max+1-k_bbl_sed))
        k_wat = (/(k,k=1,k_bbl_sed)/)       !Index vector for all points in the water column
        k_sed = (/(k,k=k_bbl_sed+1,k_max)/) !Index vector for all points in the sediments
        k_sed1 = (/(k,k=k_bbl_sed+1,k_max+1)/) !Indices of layer interfaces in the sediments (including the SWI)
    endif
   
    !Specify horizontal transport
    dx(:)= dx_adv !horizontal resolution in m
    if (h_adv.eq.1.or.h_turb.eq.1) then
        x(1)=0.0_rk
        !u_x_w(:,:,:)=0.02  !horizontal velocity in m/s if not read from forcing file...
        do i=2,i_max            !works when we have more than 2 columns
            x(i)=x(i-1)+dx(i)
        enddo
    else
        allocate(u_x_w(i_max,k_wat_bbl,days_in_yr))
        x=0.0_rk
        u_x_w=0.0_rk
    endif
    !Initialize tridiagonal matrix if necessary
    if (diff_method.gt.0) then
        call init_tridiagonal(k_max)
        write(*,*) "Initialized tridiagonal matrix"
    end if

    !Set model domain
    call fabm_set_domain(model,i_max,k_max)

    !Specify vertical index of surface and bottom
    call model%set_surface_index(k_min)
    call model%set_bottom_index(k_max)

    !Make theý (full) pressure variable to pass to FABM
    !This is used by brom_eqconst.F90 to compute equilibrium constants for pH calculations
    !and by brom_carb.F90 to compute equilibrium constants for saturation states (subroutine CARFIN)
    do i=1,i_max
        pressure(i,:) = z(:) + 10.0_rk
    end do

    !Point FABM to array slices with biogeochemical state.
    do ip=1,par_max
        call fabm_link_bulk_state_data(model, ip, cc(:,:,ip))
    end do


    !Link temperature and salinity data to FABM (needs to be redone every time julianday is updated below)
    call fabm_link_bulk_data(model, standard_variables%temperature, t(:,:,1))
    call fabm_link_bulk_data(model, standard_variables%practical_salinity, s(:,:,1))


    !Link other data needed by FABM
    call fabm_link_bulk_data(model, standard_variables%downwelling_photosynthetic_radiative_flux, Izt)  !W m-2
    call fabm_link_bulk_data(model, standard_variables%pressure, pressure)                              !dbar
    call fabm_link_horizontal_data(model, standard_variables%wind_speed, windspeed)                     !m s-1
    call fabm_link_horizontal_data(model, standard_variables%mole_fraction_of_carbon_dioxide_in_air, pco2atm)  !ppm
    if (use_hice.eq.1) call fabm_link_horizontal_data(model, ice_thickness, hice)
    call fabm_link_bulk_data(model, volume_of_cell, vv(:,:,3))


    !Check FABM is ready
    call fabm_check_ready(model)


    !Allow FABM models to use their default initialization (this sets cc)
    do k=1,k_max
        call fabm_initialize_state(model, 1, i_max, k)
    end do


    !Read initial values from ascii file if req'd
    if (port_initial_state.eq.1) call porting_initial_state_variables(trim(icfile_name), year, i_day, i_max, k_max, par_max, par_name, cc, vv)


    !Initialize output
    call init_netcdf(trim(ncoutfile_name), i_max, k_max, model, use_Eair, use_hice, year)


    !Establish which variables will be treated as solid phase in the sediments, based on the biological velocity (sinking/floating) from FABM.
    wbio = 0.0_rk
    call fabm_get_vertical_movement(model, i_water, i_water, k_wat_bbl, wbio(i_water:i_water,k_wat_bbl,:))
    do i=i_min,i_water
        wbio(i,:,:)=wbio(i_water,:,:)
    enddo
    wbio = -1.0_rk * wbio !FABM returns NEGATIVE wbio for sinking; sign change here means that wbio is POSITIVE for sinking
    is_solid = 0
    ip_sol = 0
    ip_par = 0
    write(*,*) "The following variables are assumed to join the solid phase in the sediments"
    do ip=1,par_max
      do i=i_min,i_water
        if (wbio(i,k_wat_bbl,ip).gt.0.0_rk) then
            is_solid(ip) = 1 !Any variables that SINKS (wbio>0) in the bottom cell of the water column will become "solid" in the sediments
            if(i.eq.i_min) write(*,*) trim(par_name(ip))
            if (ip_par.eq.0) ip_par = ip !Set ip_par = index of first particulate variable
        else
            if (ip_sol.eq.0) ip_sol = ip !Set ip_par = index of first solute variable
        end if
      enddo
    end do


    !Complete hydrophysical forcings
  !  cc_hmix=0.0_rk
  !if (k_points_below_water.gt.0) then
  !  call make_physics_bbl_sed(t, s, kz, hmix_rate, cc_hmix, u_x, u_x_w, t_w, s_w, kz_w, hmix_rate_w, cc_hmix_w, kz_mol, kz_bio, &
  !      z, dz, hz, i_min, i_max, i_water, k_wat_bbl, k_bbl_sed, k_max, par_max, days_in_yr, alpha, is_solid, phi, phi1, phi_inv, &
  !      pF1, pF2, mu0_musw, tortuosity, w_b, u_b, rho, dt, freq_turb, par_name, diff_method, bioturb_across_SWI, &
  !      ip_sol, ip_par, dphidz_SWI, wat_content, pWC)
  !else
        !t(:,:,:) = t_w(i_water,:,:)
        !s(:,:,:) = s_w(i_water,:,:)
        !kz(:,:,:)= kz_w(i_water,:,:)
    do i = i_min, i_water    
       do i_day=1,days_in_yr
        !Salinity (s)
        s(i,1:k_max,i_day)       = s_w(i_water,1:k_max,i_day)
        !Temperature (t)
        t(i,1:k_max,i_day)       = t_w(i_water,1:k_max,i_day)
        !Vertical diffusivity in water column (kz)
        kz(i,1:k_max,i_day)      = kz_w(i_water,1:k_max,i_day)        
       enddo
    enddo
          kz_mol = 0.0_rk
          kz_bio = 0.0_rk
    !Density of particles
    rho_def = get_brom_par("rho_def")
    rho = 0.0_rk
    do ip=1,par_max
        if (is_solid(ip).eq.1) then
            rho(ip) = get_brom_par('rho_' // trim(par_name(ip)),rho_def)
            write(*,*) "Assumed density of ", trim(par_name(ip)), " = ", rho(ip)
        end if
    end do
!  endif

    write(*,*) "Made physics of BBL, sediments"

    if(h_relax.gt.0) then   
    !Get horizontal relaxation parameters from brom.yaml:
    !hmixtype = 0, 1 or 2  for no horizontal relaxation (default), box model mixing respectively
        do ip=1,par_max
            hmixtype(i_water,ip) = get_brom_par('hmix_' // trim(par_name(ip)),0.0_rk)
            hmixtype(:,ip)=hmixtype(i_water,ip)   
                        
            if (hmixtype(i_water,ip).eq.1) then
                write(*,*) "Horizontal relaxation assumed for " // trim(par_name(ip))            
            end if
            if (hmixtype(i_water,ip).eq.2) then
    !            hmix_file = get_brom_name("hmix_filename_" // trim(par_name(ip)(13:))) !niva_oxydep_NUT
                write(*,*) "Horizontal relaxation (ASCII) assumed for " // trim(par_name(ip))
            !else if (hmixtype(i_water,ip).eq.2) then  ! read relaxation files in two cases, also for top bondary condition  #bctype_top(i_water,ip).eq.4.or.  
                open(20, file= get_brom_name("hmix_filename_" // trim(par_name(ip))))!'' // hmix_file
                write(*,*) ("hmix_filename_" // trim(par_name(ip)))
                do k=1,k_wat_bbl
                    do i_day=1,days_in_yr
                        read(20, *) i_dummy,i_dummy,cc_hmix(i_water,ip,k,i_day) ! NODC data (i_max,par_max,k_max,days_in_yr))
                    end do
                end do 
                close(20)
                do i=i_min,i_water-1
                    cc_hmix(i,:,:,:)=cc_hmix(i_water,:,:,:)
                enddo
               
            else if (bctype_top(i_water,ip).eq.4.) then    
                open(40, file = get_brom_name("bc_top_filename_" // trim(par_name(ip)))) !to boundary condition file  
                do i_day=1,days_in_yr
                    read(40, *) cc_top(i,ip,i_day)
                end do   
                close(40)
                do i=i_min,i_water-1
                    cc_top(i,:,:)=cc_top(i_water,:,:)
                enddo
            end if
        end do
    endif
    
    open(8,FILE = 'temp.dat')
    !Set constant forcings
    wind_speed = get_brom_par("wind_speed")    ! 10m wind speed [m s-1]
    pco2_atm   = get_brom_par("pco2_atm")      ! CO2 partical pressure [ppm]
    pco2atm(:) = pco2_atm                      ! Same pCO2 values for all the surface gridpoints
    windspeed(:) = wind_speed                  ! Same wind speed values to all the surface gridpoints

    end subroutine init_brom_transport
!=======================================================================================================================






!=======================================================================================================================
    subroutine do_brom_transport()

    !Executes the offline vertical transport model BROM-transport

    use calculate, only: calculate_light, calculate_phys, calculate_sed, calculate_sed_eya

    implicit none

    integer      :: id, idt, idf                     !time related
    integer      :: surf_flux_with_diff              !1 to include surface fluxes in diffusion update, 0 to include in bgc update
    integer      :: bott_flux_with_diff              !1 to include bottom fluxes in diffusion update, 0 to include in bgc update
    integer      :: model_w_sed                      !1 to assume porosity effects for solutes and solids, 0 - sumplified approach
    integer      :: constant_w_sed                   !1 to assume constant burial (advection) velocities in the sediments
    integer      :: dynamic_w_sed                    !1 to assume dynamic burial (advection) velocities in the sediments depending on dVV(k_bbl_sed)
    integer      :: show_maxmin, show_kztCFL, show_wCFL, show_nan, show_nan_kztCFL, show_nan_wCFL     !options for runtime output to screen
    integer      :: bc_units_convert, sediments_units_convert !options for conversion of concentrations units in the sediment 
    integer      :: julianday, model_year
    integer      :: k_inj,i_inj,inj_swith,inj_num,start_inj,stop_inj    !#number of layer and column to inject into, start day, stop day number 
    integer      :: i_solid, i_air                   !#number of columns for solid and air volumes 
    integer      :: i_flag                           !Urmia related
    real(rk)     :: cnpar                            !"Implicitness" parameter for GOTM vertical diffusion (set in brom.yaml)
    real(rk)     :: cc0                              !Resilient concentration (same for all variables)
    real(rk)     :: omega                            !angular frequency of sinusoidal forcing = 2*pi/365 rads/day
    real(rk)     :: O2stat                           !oxygen status of sediments (factor modulating the bioirrigation rate)
    real(rk)     :: a1_bioirr                        !to detect whether or not bioirrigation is activated
    real(rk)     :: hmix_rate_uniform                !uniform horizontal relaxation if prescribed
    real(rk)     :: fresh_PM_poros                   ! porosity of fresh precipitated PM (i.e. dVV) 
    real(rk)     :: max_poros                   ! threshold porosuty to stop sinkin into the botom layer
    real(rk)     :: w_binf                   ! 
    real(rk)     :: bu_co                   ! "Burial coeficient" for setting velosity exactly to the SWI proportional to the 
                                           !   settling velocity in the water column (0<bu_co<1), 0 - for no setting velosity, (nd)
    real(rk)     :: injection_rate_1, injection_rate_2, injection_rate_3 ! injection rate
    real(rk)     :: vol_PM ! volume of PM in 1 l
    real(rk), parameter :: pi=3.141592653589793_rk
    !character(len=attribute_length), allocatable, dimension(:)    :: inj_var_name
    character(len=attribute_length):: inj_var_name_1
    character(len=attribute_length):: inj_var_name_2
    character(len=attribute_length):: inj_var_name_3
    character(300)                 :: comments

    omega   = 2.0_rk*pi/365.0_rk
    i_solid = 2
    i_air   = 3

    !Get parameters for the time-stepping and vertical diffusion / sedimentation
    cnpar = get_brom_par("cnpar")
    model_w_sed = get_brom_par("model_w_sed")
    dynamic_w_sed = get_brom_par("dynamic_w_sed")
    constant_w_sed = get_brom_par("constant_w_sed")
    fresh_PM_poros = get_brom_par("fresh_PM_poros")
    max_poros = get_brom_par("max_poros")
    w_binf = get_brom_par("w_binf")
    bu_co = get_brom_par("bu_co")
    cc0 = get_brom_par("cc0")
    a1_bioirr = get_brom_par("a1_bioirr")
    surf_flux_with_diff = get_brom_par("surf_flux_with_diff")
    bott_flux_with_diff = get_brom_par("bott_flux_with_diff")
    show_maxmin = get_brom_par("show_maxmin")
    show_kztCFL = get_brom_par("show_kztCFL")
    show_wCFL = get_brom_par("show_wCFL")
    show_nan = get_brom_par("show_nan")
    show_nan_kztCFL = get_brom_par("show_nan_kztCFL")
    show_nan_wCFL = get_brom_par("show_nan_wCFL")
    bc_units_convert = get_brom_par("bc_units_convert")
    sediments_units_convert = get_brom_par("sediments_units_convert")
    hmix_rate_uniform = get_brom_par("hmix_rate_uniform")
    
    inj_swith = get_brom_par("injection_swith")
    inj_var_name_1 = get_brom_name("inj_var_name_1")
    injection_rate_1 = get_brom_par("injection_rate_1")
    inj_var_name_2 = get_brom_name("inj_var_name_2")
    injection_rate_2 = get_brom_par("injection_rate_2")
    inj_var_name_3 = get_brom_name("inj_var_name_3")
    injection_rate_3 = get_brom_par("injection_rate_3")
    k_inj = get_brom_par("k_injection")
    i_inj = get_brom_par("i_injection") 
    start_inj = get_brom_par("start_inj")
    stop_inj = get_brom_par("stop_inj")    

    idt = int(1._rk/dt)                                      !number of cycles per day
    model_year = 0
    kzti = 0.0_rk

    dVV(:,:,:) = 0.0_rk
    volume(:,:) = 0.0_rk
    area(:,:) = 0.0_rk

! READ  Initial volumes and areas of layers, they have "i_max+1" dimension for 1-sediment, 2-water, 3-air:
      open(44, file = "Evolume.dat") 
        do k=k_max,1,-1
          read(44, *) volume(i_water,k) !initially all the layer'volumes is filled with water
!           volume(i_water,k)= 500000000.0_rk
          area(i_water,k)=volume(i_water,k)/dz(k)
        end do   
        !repeat again to read the shallower 50 depths 
        !do k=k_max,1,-1
        !  read(44, *) volume(i_water,k) !initially all the layer'volumes is filled with water
        !  area(:,k)=volume(i_water,k)/dz(k)
        !end do   
        volume(i_water,k_max+1)=0.0_rk
        area(:,k_max+1)=0.0_rk
        hz=dz ! at initial conditions thickness of the box is the same as the distance between their centers!
      close(44)
!          Smoothing is needed....
       do i=1,5
           do k = 2,k_max 
               volume(i_water,k)=(volume(i_water,k-1)+volume(i_water,k)+volume(i_water,k+1))/3.0_rk
               area(2,k)=(area(2,k-1)+area(2,k)+area(2,k+1))/3.0_rk
           end do   
       enddo

! LOAD RIVER CHEMICAL DISCHARGES DATA
      open(48, file = "chem_discharhes_new.dat") 
          read(48, *) comments
!          write (*,*) comments
      do i=1, 6570
          read(48, *) rv_PO4(i),rv_NO3(i),rv_NO2(i),rv_O2(i),rv_Si(i),rv_DON(i),rv_Na(i),rv_K(i),rv_Mg(i),rv_Ca(i),rv_SO4(i),rv_Cl(i),rv_CO3(i),rv_DIC(i),rv_NH4(i),rv_Phy(i),rv_Br(i),rv_TOM(i),rv_Alk(i)
!          write (*,*) rv_PO4(i),rv_NO3(i),rv_NO2(i),rv_O2(i),rv_Si(i),rv_DON(i),rv_Na(i),rv_K(i),rv_Mg(i),rv_Ca(i),rv_SO4(i),rv_Cl(i),rv_CO3(i),rv_DIC(i),rv_NH4(i),rv_Phy(i),rv_Br(i),rv_TOM(i),rv_Alk(i)
      end do   

!        po4	no3	no2	o2	sio2	tn	na	k	mg	ca	so4	cl	co3	hco3	nh4	chl_a	br	toc	alk

        do i=1, 6935
            rv_DON(i)=max(0.0_rk,rv_TOM(i)*0.25_rk*16.0_rk/106.0_rk-rv_NO3(i)-rv_NO2(i)-rv_NH4(i)) !calculate from TOC
            rv_NH4(i)=0.0_rk !#rv_NH4(i)/2._rk  !don't believe their technique
            rv_DIC(i)=rv_DIC(i)+rv_CO3(i)   ! from mg Chl/l into uM N = 1.*0.26
            rv_Phy(i)=rv_Phy(i)*14.0_rk*0.26_rk ! transfer back to Chl-a and than to uM N  
            rv_SO4(i)=rv_SO4(i) *3.0_rk   ! coreect for units
        end do   
      close(48)

! LOAD LAKE SURFACE CHANGES DATA
      open(50, file = "Ezeta_casp3.dat") 
        do i=1, 9855 ! < 19 yrs(6935 days) real + >8 yrs forcast, i.e. 27 years
!          read(48, *) comments
          read(50, *) zeta(i)
!          write (*,*) i, zeta(i)
        end do

       do k=1,5
!           do i=2, 6934
           do i = 2,9854 
    !          Smoothing is needed....
               dzeta(i)=(zeta(i-1)+zeta(i)+zeta(i+1))/3.0_rk
           end do   
           zeta(:)=dzeta(:)
       enddo
        
        dzeta(1)= 0.0_rk
        do i=2, 9855
          dzeta(i)= -1.0_rk*(zeta(i+100)-zeta(i+99))
!          dzeta(i)= -1.0_rk*(zeta(i)-zeta(i-1))          
        end do   
      close(50)

!      extrapolate river's chemistry discharge for the forecast period
      do i=6571, 9855 !<<<<<<<<<<<<---------------
          rv_PO4(i) = rv_PO4(i-365)
          rv_NO3(i) = rv_NO3(i-365)
          rv_NO2(i) = rv_NO2(i-365)
          rv_O2(i) = rv_O2(i-365)
          rv_Si(i) = rv_Si(i-365)
          rv_DON(i) = rv_DON(i-365)
          rv_Na(i) = rv_Na(i-365)
          rv_K(i) = rv_K(i-365)
          rv_Mg(i) = rv_Mg(i-365)
          rv_Ca(i) = rv_Ca(i-365)
          rv_SO4(i) = rv_SO4(i-365)
          rv_Cl(i) = rv_Cl(i-365)
          rv_DIC(i) = rv_DIC(i-365)
          rv_NH4(i) = rv_NH4(i-365)
          rv_Phy(i) = rv_Phy(i-365)
          rv_Br(i) = rv_Br(i-365)
          rv_Alk(i) = rv_Alk(i-365)
      end do   

        !convert bottom boundary values from 'mass/pore water ml' for dissolved and 'mass/mass' for solids into 'mass/total volume'
    if (bc_units_convert.eq.1) then
        do i=i_min,i_water
            do ip=1,par_max
                if (bctype_bottom(i,ip).eq.1) bc_bottom(i,ip) = bc_bottom(i,ip)/pF1(i,k_max,ip)
            enddo
          enddo    
    end if
    
    !Uniform horizontal relaxation if prescribed
    if(hmix_rate_uniform>0.0_rk)    hmix_rate=hmix_rate_uniform

    !Master time step loop over days
    write(*,*) "Starting time stepping"
    k_air_wat = 15
        do k=1,k_air_wat-1
            volume(i_air,k)=volume(i_water,k)
            volume(i_water,k)=0.0_rk
            cc(:,k,:)=0.0_rk
        enddo
    i_flag =0

 !_______BIG Cycle ("i_day"=0,...,last_day-1)________!

  do i_day=0,(last_day-1)

    julianday = i_day - int(i_day/days_in_yr)*days_in_yr + 1    !"julianday" (1,2,...,days_in_yr)
    if (julianday==1) model_year = model_year + 1
    write (*,'(a, i4, a, i4)') " model year:", model_year, "; julianday:", julianday
 !Calculate Izt = <PAR(z)>_24hr for this day
    do i=i_min, i_water
      call calculate_light(julianday, i, k_bbl_sed, k_max, par_max, hz, Eair, use_Eair, hice, cc, is_solid, rho, Izt)
    enddo
 ! Reload daily cheanges in t and s
    call fabm_link_bulk_data(model, standard_variables%temperature, t(:,:,julianday))
    call fabm_link_bulk_data(model, standard_variables%practical_salinity, s(:,:,julianday))

 ! calculate salinity as a function of main ions (convert from umol/l into g/l -> CC*(m.v.)/1000/1000)
    do i=i_min, i_water
      do k=1,k_max
        s(i,k,julianday) = min(900.0_rk,cc(i,k,id_Cl)*35.453e-6 + cc(i,k,id_SO4)*96.e-6 &
          +cc(i,k,id_DIC)*61.e-6 + cc(i,k,id_Na)*22.990e-6 + cc(i,k,id_Mg)*24.305e-6 &
          +cc(i,k,id_K)*39.098e-6 + cc(i,k,id_Ca_u)*40.078e-6)
      enddo
    enddo
        
 ! Set time-varying Dirichlet boundary conditions for current julianday
    do ip=1,par_max
      do i=i_min, i_water
 ! Sinusoidal variations
        if (bctype_top(i,ip).eq.2) bc_top(i,ip) = bcpar_top(i,ip,1) &
          +bcpar_top(i,ip,2)*sin(omega*(julianday-bcpar_top(i,ip,3)))
        if (bctype_bottom(i,ip).eq.2) bc_bottom(i,ip) = bcpar_bottom(i,ip,1) &
          +bcpar_bottom(i,ip,2)*sin(omega*(julianday-bcpar_bottom(i,ip,3)))
 ! Variations read from netcdf
        if (bctype_top(i,ip).eq.3) bc_top(i,ip) = cc_top(i,ip,julianday)
        if (bctype_bottom(i,ip).eq.3) bc_bottom(i,ip) = cc_bottom(i,ip,julianday)
 ! Variations read from ascii file and/or calculated as a function of something
        if (bctype_top(i,ip).eq.4) bc_top(i,ip) = cc_hmix(i,ip,1,julianday)
 ! Calculations of Alk, SO4 from salinity
                !SO4 in mmol/m3, SO4/Salt from Morris, A.W. and Riley, J.P.(1966) quoted in Dickson et al.(2007)              
!        if (bctype_top(i,ip).eq.5) bc_top(i,ip)=(0.1400_rk/96.062_rk)*(s(1,1,julianday)/1.80655_rk)*1.e6_rk !.and.ip.eq.id_SO4
!        if (bctype_top(i,ip).eq.4.and.ip.eq.id_Alk) bc_top(i,ip)=0.068*s(1,1,julianday)  !Alk in mmol/m3, Alk/Salt from Murray, 2014
      enddo
    enddo

!####################################################
 !Subloop over timesteps in the course of one day
 !Note: The numerical approach here is Operator Splitting with tracer transport processes assumed to be
 !numerically more demanding than the biogeochemistry (hence freq_turb, freq_sed >= 1) (Butenschon et al., 2012)

    do id=1,idt
      k_min=k_air_wat+5
     call model%set_surface_index(k_min)
     call model%set_bottom_index(k_max)
      dVV(:,:,1)=0.0
      dvolume=0.0_rk
      dcc=0.0_rk

!===================================================
!_______calculate Kz and Wsed________! 

 !Compute vertical sinking velocity in water column (sinking/floating) using the FABM.
      wbio = 0.0_rk
      do k=1,k_max
        call fabm_get_vertical_movement(model, i_water, i_water, k, wbio(i_water:i_water,k,:))  !Note: wbio is on layer midpoints
      end do
      wbio = -1.0_rk * wbio !FABM returns NEGATIVE wbio for sinking; sign change here means that wbio is POSITIVE for sinking
      wti = 0.0_rk
 ! change values of Kz and Wd in case of formation of sediment formation.
      do k=1,70 !k_air_wat+20
        wti(i_water,k,:)        = wbio(i_water,k,:)
        kz(i_water,k,julianday) = kz_w(i_water,k,julianday)
!        kz(i_solid,k,julianday) = kz_w(i_water,k,julianday)
      enddo           
      do k=71, k_max !k_air_wat+21,k_max
!            if(volume(i_solid,k)/volume(i_water,k).lt.(1.0_rk-max_poros)) then
        if(volume(i_solid,k)/(volume(i_water,k)+volume(i_solid,k)).lt.(1.0_rk-max_poros)) then !!!! if water
          wti(i_water,k,:)        = wbio(i_water,k,:)
          kz(i_water,k,julianday) = kz_w(i_water,k,julianday)
!          kz(i_solid,k,julianday) = kz_w(i_water,k,julianday)
        else                                                        !!!! if sediment
          wti(i_water,k,:)        = 0.0_rk !
          wti(i_water,k-1,:)      = 0.0_rk !
          kz(i_water,k,julianday) = 1.0e-11 !
!          kz(i_solid,k,julianday) = 0.0_rk !
!          if(volume(i_solid,k-1)/(volume(i_water,k-1)+volume(i_solid,k-1)).lt.(1.0_rk-max_poros)) then !!!! if water
!            wti(i_water,k-1,:)      =  wbio(i_water,k,:)
!            kz(i_water,k,julianday) = kz_w(i_water,k,julianday)
!            kz(i_solid,k,julianday) = kz_w(i_water,k,julianday)/2.0_rk
!          endif
        endif
      enddo
!      wti(i_water,k_max,:)   = 0.0_rk 
      wti(i_water,k_max+1,:) = 0.0_rk !1.0E-10  ! for Urmia !we need to have burial velosity to remove nutrients and salts from the model domain
!      kz(i_solid,k_max,julianday) = 0.0_rk 
!      kz(i_solid,k_max+1,julianday) = 0.0_rk 
      kz(i_water,k_max,julianday) = 0.0_rk 
      kz(i_water,k_max+1,julianday) = 0.0_rk 
!===================================================
!_______air/water interface and volumes________!

      dvolume(i_air,:) = 0.0_rk ! change of "air" volume we used for air/water interface, function of dzeta
 !-------------------
      if (dzeta(i_day).ge.0.0_rk) then   ! dzeta is in [m/day], POSITIVE when level goes DOWN vvvv
 !-------------------
        if (hz(k_air_wat).ge.(hz(k_air_wat+1)/max_hz_inc)) then !.or.volume(i_water,k_air_wat).lt.0.0_rk) then ! if the top layer very thin, we merge it with the following deeper one
          dvolume(i_air,k_air_wat) = dt*dzeta(i_day)*area(i_water,k_air_wat) ! change of volume
          hz(k_air_wat) = hz(k_air_wat)-dt*dzeta(i_day)           !change of thickness
 ! change of concentration
          cc(i_water,k_air_wat,:)=cc(i_water,k_air_wat,:)*volume(i_water,k_air_wat)/(volume(i_water,k_air_wat)-dvolume(i_air,k_air_wat))
          volume(i_water,k_air_wat) = volume(i_water,k_air_wat) - dvolume(i_air,k_air_wat) ! water
          if (hz(k_air_wat).le.hz(k_air_wat+2)) then
            volume(i_air-1,  k_air_wat) = volume(i_air-1,  k_air_wat) + dvolume(i_air,k_air_wat) ! air
          else
            volume(i_air,  k_air_wat) = volume(i_air,  k_air_wat) - dvolume(i_air,k_air_wat) ! air
          endif
        else ! we MERGE the thin top layer the following deeper one
          dvolume(i_air,k_air_wat) = dt*dzeta(i_day)*area(i_water,k_air_wat) ! change of volume
          hz(k_air_wat) = hz(k_air_wat)+dt*dzeta(i_day)           !change of thickness
          cc(i_water,k_air_wat,:)=cc(i_water,k_air_wat,:)*volume(i_water,k_air_wat)/(volume(i_water,k_air_wat)-dvolume(i_air,k_air_wat))
          cc(i_water,k_air_wat+1,:)=(cc(i_water,k_air_wat+1,:)*volume(i_water,k_air_wat+1)+cc(i_water,k_air_wat,:)*volume(i_water,k_air_wat))/ (volume(i_water,k_air_wat+1)+volume(i_water,k_air_wat))  !here we merge 
          volume(i_water,k_air_wat+1)=volume(i_water,k_air_wat+1)+hz(k_air_wat)*area(i_water,k_air_wat) !  increase the WATER volume() of the next(deeper) layer
          volume(i_air,k_air_wat)= min((hz(k_air_wat+2)-hz(k_air_wat))*area(i_water,k_air_wat),hz(k_air_wat)*area(i_water,k_air_wat))  ! we increase the AIR volume of the previous(upper) layer
          volume(i_water,k_air_wat) = 0.0_rk  !-hz(k_air_wat)*area(i_water,k_air_wat)  ! and set to =0. volume of previous(upper) layer 
          cc(i_water,k_air_wat,:) = 0.0_rk  ! and concentration of previous(upper) layer 
          hz(k_air_wat) = hz(k_air_wat+2) !0.0_rk  ! and thickness of previous(upper) layer sjoud be standard
          hz(k_air_wat+1)=volume(i_water,k_air_wat+1)/area(i_water,k_air_wat+1)  ! new next(deeper) layer thickness
          k_air_wat=k_air_wat+1 ! finally, change index
!          i_flag=1
        endif
!-------------------
      else                             !  dzeta() is NEGATIVE, the level goes UP  ^^^^^
!-------------------
        if (hz(k_air_wat).lt.(hz(k_air_wat+2)*(1.0_rk+1.0_rk/real(max_hz_inc)))) then !.or.volume(i_water,k_air_wat).lt.0.0_rk) then ! if the top layer very thin, we merge it with the following deeper one
          dvolume(i_air,k_air_wat) = dt*dzeta(i_day)*area(i_water,k_air_wat)  ! change of volume, negative
          hz(k_air_wat) = hz(k_air_wat)-dt*dzeta(i_day)  !change of thickness   !decreases
          cc(i_water,k_air_wat,:)=cc(i_water,k_air_wat,:)*volume(i_water,k_air_wat)/(volume(i_water,k_air_wat)-dvolume(i_air,k_air_wat)) !change of concentration
          volume(i_water,k_air_wat) = volume(i_water,k_air_wat) - dvolume(i_air,k_air_wat) ! water
          if (hz(k_air_wat).le.hz(k_air_wat+2)) then
            volume(i_air,  k_air_wat) = volume(i_air,  k_air_wat) + dvolume(i_air,k_air_wat) ! air ?????
          else
            volume(i_air,  k_air_wat-1) = volume(i_air,  k_air_wat-1) + dvolume(i_air,k_air_wat) ! air ?????
          endif
        else ! dzeta() is NEGATIVE, MAKE NEW LEVEL; dzeta() is negative
          dvolume(i_air,k_air_wat) = dt*dzeta(i_day)*area(i_water,k_air_wat)
          hz(k_air_wat) = hz(k_air_wat)-dt*dzeta(i_day)               !change of thickness
          cc(i_water,k_air_wat,:)=cc(i_water,k_air_wat,:)*volume(i_water,k_air_wat)/(volume(i_water,k_air_wat)-dvolume(i_air,k_air_wat)) !change of concentration
          cc(i_water,k_air_wat-1,:)=cc(i_water,k_air_wat,:) ! we set the same concentration
          volume(i_water,k_air_wat) = volume(i_water,k_air_wat) - dvolume(i_air,k_air_wat) ! water
 ! now fill in the new shallower layer
          volume(i_water,k_air_wat-1) = volume(i_water,k_air_wat) - area(i_water,k_air_wat)*hz(k_air_wat+2)
          volume(i_water,k_air_wat) = area(i_water,k_air_wat)*hz(k_air_wat+2)
          hz(k_air_wat)=hz(k_air_wat+2)
          hz(k_air_wat-1)=volume(i_water,k_air_wat-1)/area(i_water,k_air_wat-1)
          volume(i_air,  k_air_wat) = 0.0_rk ! air
          volume(i_air,  k_air_wat-1) = area(i_water,k_air_wat-1)*hz(k_air_wat+1)- volume(i_water,k_air_wat-1)! air
          k_air_wat=k_air_wat-1 !new index
        endif
!-------------------
      endif

!===================================================
!_______vertical mixing________!

      fick = 0.0_rk !NEW
 ! Air-sea interface
      fick(i_water,k_air_wat,:) = 0.0_rk !surf_flux(i_water,:)
 ! calculate Fick fluxes between the boxes
      do k=k_air_wat+1,k_max
        do ip=1,par_max !Sum over contributions from each particulate variable
 ! calculate new volume  occupied by particulate matter
!          if (is_solid(ip).eq.1) then
!            fick(i_water,k,:) = dt*86400.0_rk*kz(i_solid,k,julianday)*area(i_water,k)*(cc(i_water,k-1,:)-cc(i_water,k,:))/(0.5_rk*(hz(k-1)+hz(k)))
!          else
            fick(i_water,k,:) = dt*86400.0_rk*kz(i_water,k,julianday)*area(i_water,k)*(cc(i_water,k-1,:)-cc(i_water,k,:))/(0.5_rk*(hz(k-1)+hz(k)))
!          end if
        end do
        fick(i_water,k,:) = dt*86400.0_rk*kz(i_water,k,julianday)*area(i_water,k)*(cc(i_water,k-1,:)-cc(i_water,k,:))/(0.5_rk*(hz(k-1)+hz(k)))
      end do
        fick(i_water,k_max+1,:) = 0.0_rk
 ! calculate changes in time and integrate
      do k=k_air_wat,k_max
        dcc(i_water,k,:) = (fick(i_water,k,:)-fick(i_water,k+1,:))/volume(i_water,k)
      end do
 ! cc of boxes
      do k=k_air_wat,k_max
        cc(i_water,k,:) = cc(i_water,k,:) + dcc(i_water,k,:)
      end do
 ! calculate changes of volume of water due to transport of particles during mixing
!      do k=1,k_max
!        do ip=1,par_max !Sum over contributions from each particulate variable
!          if (is_solid(ip).eq.1) then
!            dvolume(i_water,k)= dvolume(i_water,k)+(fick(i_water,k,ip)-fick(i_water,k,ip))/rho(ip) ! dcc(i,k-1,ip)/rho(ip) should be 
!          endif
!        enddo
!      enddo

!===================================================
!_____water_biogeochemistry_______!

     dcc = 0.0_rk
     do i=i_water, i_water
       do k=k_air_wat,k_max
         call fabm_do(model, i, i, k, dcc(i:i,k,:))   
 ! Note: We MUST pass the range "i:i" to fabm_do -- a single value "i" will produce compiler error
       end do
 ! Add surface and bottom fluxes if treated here
       if (surf_flux_with_diff.eq.0) then
         surf_flux = 0.0_rk
         call fabm_do_surface(model, i, i, surf_flux(i:i,:))
         fick(i,k_min,:) = surf_flux(i,:)
         do ip=1,par_max
           dcc(i,k_air_wat+1,ip) = dcc(i,k_air_wat+1,ip) + surf_flux(i,ip) / hz(k_air_wat+1)
         end do
       end if
       if (bott_flux_with_diff.eq.0) then
         bott_flux = 0.0_rk
         bott_source = 0.0_rk                
         call fabm_do_bottom(model, i, i, bott_flux(i:i,:),bott_source(i:i,:))
         sink(i,k_max+1,:) = bott_flux(i,:)
         do ip=1,par_max
           dcc(i,k_max,ip) = dcc(i,k_max,ip) + bott_flux(i,ip) / hz(k_max)
         end do
       end if
       dcc_R(i,:,:) = dcc(i,:,:) !for Urmia if (dynamic_w_sed.eq.1) dcc_R(i,:,:) = dcc(i,:,:) !Record biological reaction terms for use in calculate_sed
 ! Euler time step due to FABM biogeochemistry
       do k=1,k_max
         do ip=1,par_max
           cc(i,k,ip) = cc(i,k,ip) + 86400.0_rk*dt*dcc(i,k,ip)
         end do
       end do
     enddo
     i=i_water
 ! Reassert Dirichlet BCs
     do ip=1,par_max
       if (bctype_top(i,ip).gt.0) then
         cc(i,k_min,ip) = bc_top(i,ip)
       end if
       if (bctype_bottom(i,ip).gt.0) then
         cc(i,k_max,ip) = bc_bottom(i,ip)
       end if
     end do
       cc(i,:,:) = max(cc0, cc(i,:,:)) !Impose resilient concentration

!===================================================
 !_______Changes of volumes of "water"/"solid" boxes_________!

      cc(:,:,id_V_air)=volume(i_air,k)
      cc(:,:,id_V_wat)=volume(i_water,k)
      cc(:,:,id_V_sed)=volume(i_solid,k)
      cc(:,:,id_dV_ch)= 0.0_rk
      cc(:,:,id_dV_sink)= 0.0_rk

!!!!      do k=k_air_wat,k_max
!!!!        do ip=1,par_max !Sum over contributions from each particulate variable
!!!!!          if (is_solid(ip).eq.1) then
!!!!!              cc(i_water,k,id_dV_ch)=cc(i_water,k,id_dV_ch) +( &
!!!!!                  + 86400.0_rk*dt*dcc_R(i_water,k,ip)/rho(ip)*volume(i_water,k) &
!!!!!                   ) 
!!!!!              cc(i_water,k,id_dV_sink)=cc(i_water,k,id_dV_sink) +( &
!!!!!                  + sink(i_water,k-1,ip)*area(i_water,k)/area(i_water,k-1)-sink(i_water,k,ip) &
!!!!!                   ) /rho(ip)
!!!!! ! calculate changes of volume of water due to mixing (done before) + biology + sinking (from the previous time step)
!!!!! ! dvolume(i_water,k)= dvolume(i_water,k)+(dcc_R(i_water,k,ip)*volume(i_water,k) &
!!!!! !        + sink(i_water,k-1,ip)*area(i_water,k)/area(i_water,k-1)-sink(i_water,k,ip))/rho(ip)
!!!!!!            dvolume(i_water,k)= dvolume(i_water,k) +( &
!!!!!!                  + dcc_R(i_water,k,ip)*volume(i_water,k) &
!!!!!!!                  + (fick(i_water,k,ip)-fick(i_water,k-1,ip)) &
!!!!!!                  + sink(i_water,k-1,ip)*area(i_water,k)/area(i_water,k-1)-sink(i_water,k,ip) &
!!!!!!                   ) /rho(ip)
!!!!!          end if
!!!!          if (ip.eq.id_NaCl) then
!!!!              cc(i_solid,k,id_dV_ch)=cc(i_solid,k,id_dV_ch) +( &
!!!!                  + 86400.0_rk*dt*dcc_R(i_water,k,ip)/rho(ip)*(volume(i_water,k)+volume(i_solid,k)) &
!!!!                   )
!!!!              cc(i_solid,k,id_dV_sink)=cc(i_solid,k,id_dV_sink) +( &
!!!!                  + sink(i_water,k-1,ip)*area(i_water,k)/area(i_water,k-1)-sink(i_water,k,ip) &
!!!!                   ) /rho(ip)
!!!!          end if
!!!!        end do
!!!!!        dvolume(i_water,k)=cc(i_water,k,id_dV_ch)+cc(i_water,k,id_dV_sink)
!!!!        dvolume(i_water,k)=cc(i_solid,k,id_dV_ch)+cc(i_solid,k,id_dV_sink)
!!!!      enddo
!!!!      do k=k_air_wat,k_max
!!!!        volume(i_solid,k)=max(0.0_rk,volume(i_solid,k)+dvolume(i_water,k)) ! particulate
!!!!        volume(i_water,k)=max(0.0_rk,volume(i_water,k)-dvolume(i_water,k)) ! water
!!!!      enddo
!!!      do k=k_air_wat,k_max
!!! ! (from the previous time step)
!!!        volume(i_water,k)=volume(i_water,k)+volume(i_solid,k)
!!!        volume(i_solid,k)=0.0_rk
!!! !              vol_PM=0.0_rk
!!!        do ip=1,par_max !Sum over contributions from each particulate variable
!!! ! calculate new volume  occupied by particulate matter
!!!          if (is_solid(ip).eq.1) then
!!!!            volume(i_solid,k)= volume(i_solid,k)+cc(i_water,k,ip)*volume(i_water,k)/rho(ip)
!!!            volume(i_solid,k)= volume(i_solid,k)+(cc(i_water,k,ip)+cc(i_solid,k,ip))*volume(i_water,k)/rho(ip)
!!!!                     vol_PM= vol_PM+(cc(i_water,k,ip)+cc(i_water,k,ip))/rho(ip)
!!!          end if
!!!        end do
!!! ! re-calculate new volume of water
!!!        volume(i_water,k)=volume(i_water,k)-volume(i_solid,k) ! water
!!! !            volume(i_solid,k)=vol_PM*volume(i_water,k)
!!!!             volume(i_water,k)=volume(i_water,k)-volume(i_solid,k) ! water
!!!      enddo
!       if (id.eq.1) write (8,'(a, i4, a, i4, a, e)') " model year:", model_year, "; julianday:", julianday, "; dVV(k_bbl_sed):", dVV(1,k_bbl_sed,1)
           !Air (total-others)/Water/Sediment volumes of layers_______!
     do k=k_air_wat,k_max
 ! (from the previous time step)
        volume(i_water,k)=volume(i_water,k)+volume(i_solid,k)
        volume(i_solid,k)=0.0_rk
 !              vol_PM=0.0_rk
        do ip=1,par_max !Sum over contributions from each particulate variable
 ! calculate new volume  occupied by particulate matter
          if (is_solid(ip).eq.1) then
!            volume(i_solid,k)= volume(i_solid,k)+cc(i_water,k,ip)*volume(i_water,k)/rho(ip)
            volume(i_solid,k)= volume(i_solid,k)+(cc(i_water,k,ip))*volume(i_water,k)/rho(ip)
!                     vol_PM= vol_PM+(cc(i_water,k,ip)+cc(i_water,k,ip))/rho(ip)
          end if
        end do
 ! re-calculate new volume of water
        volume(i_water,k)=volume(i_water,k)-volume(i_solid,k) ! water
 !            volume(i_solid,k)=vol_PM*volume(i_water,k)
!             volume(i_water,k)=volume(i_water,k)-volume(i_solid,k) ! water
     enddo

     do k=1,k_max
 ! in m3
        !cc(2,k,id_V_air)=volume(i_air,k)
        !cc(2,k,id_V_wat)=volume(i_water,k)
        !cc(2,k,id_V_sed)=volume(i_solid,k)
        !cc(2,k,id_hz_m)=hz(k)
        !cc(2,k,id_z_m)=z(k)
        !cc(2,k,id_area)=area(i_water,k)
 ! in %
        cc(1,k,id_V_air)=min(100.0_rk,max(0.0_rk,volume(i_air,  k)*100.0_rk/(volume(i_air,k)+volume(i_water,k)+volume(i_solid,k))))
        cc(1,k,id_V_wat)=max(0.0_rk,min(100.0_rk,volume(i_water,k)*100.0_rk/(volume(i_air,k)+volume(i_water,k)+volume(i_solid,k))))
        cc(1,k,id_V_sed)=max(0.0_rk,min(100.0_rk,volume(i_solid,k)*100.0_rk/(volume(i_air,k)+volume(i_water,k)+volume(i_solid,k))))
        cc(1,k,id_hz_m)=100.0_rk*hz(k)/hz(k_max)
        cc(1,k,id_z_m)=z(k)
        cc(1,k,id_area)=max(0.0_rk,100.0_rk -100.0_rk*volume(i_solid,k)/max(0.000001,volume(i_water,k))) !porosity○
      enddo

!===================================================
 !_______Particles sinking_________!

              sink = 0.0_rk
              dcc = 0.0_rk
           !Water column and sediment layer interfaces
           do ip=1,par_max
             if (is_solid(ip).eq.1) then
            !!Calculate fluxes at the boundaries between the boxes, sink(in mmol/(boundary-between-layers)/sec )
            !! Sunk = Cc() * (Wc *Area / Vol) *dt
               do k=k_air_wat+1,k_max
                  sink(i_water,k,ip) = cc(i_water,k,ip)*86400.0_rk*dt*wti(i_water,k,ip)*area(i_water,k) !/volume(i_water,k)
!                  sink(i_water,k,ip) = volume(i_water,k)*cc(i_water,k,ip)*86400.0_rk*dt*wti(i_water,k,ip)/(0.5_rk*(hz(k-1)+hz(k)))
               enddo
                  sink(i_water,k_air_wat,ip)=0.0_rk
            !!Calculate tendencies dcc = dcc/dt = -dF/dz in the "water" boxes (#2) and "solid" boxes (#1)
               do k=k_air_wat+1,k_max-1
!                    dcc(i_water,k,ip) = (sink(i_water,k-1,ip)* area(i_water,k)/area(i_water,k-1) - sink(i_water,k,ip))/volume(i_water,k)
                    dcc(i_water,k,ip) = (sink(i_water,k-1,ip)- sink(i_water,k,ip))/volume(i_water,k)
                   ! dcc(i_solid,k,ip) = 0.0_rk ! (sink(i_water,k-1,ip)*(area(i_water,k-1)-area(i_water,k))/area(i_water,k-1))/ (max(0.01_rk, volume(i_water,k)))
               end do
               ! from the top layer just sink down
!                    dcc(i_water,1,ip)= - sink(i_water,1,ip)/volume(i_water,1)
               ! from the deepest layer no sink down 
                    dcc(i_water,k_max,ip) = sink(i_water,k_max-1,ip) /volume(i_water,k_max)
            !!Time integration
                do i=1,i_water
                    do k=1,k_max
                        cc(i,k,ip) = cc(i,k,ip) + dcc(i,k,ip)
                    end do
                end do
              endif
           enddo

!===================================================
!________Injection____________________!

 !            !Source of "acetate" 1292 mmol/sec, should be devided to the volume of the grid cell, i.e. dz(k)*dx(i)*dx(i)
                i_inj=1
                k_inj=k_air_wat+2
            cc(i_inj,k_inj,id_PO4) = cc(i_inj,k_inj,id_PO4) &
                            +86400.0_rk*dt*rv_PO4(i_day)/volume(i_inj,k_inj)
            cc(i_inj,k_inj,id_NO3) = cc(i_inj,k_inj,id_NO3) &
                            +86400.0_rk*dt*rv_NO3(i_day)/volume(i_inj,k_inj)
            cc(i_inj,k_inj,id_NO2) = cc(i_inj,k_inj,id_NO2) &
                            +86400.0_rk*dt*rv_NO2(i_day)/volume(i_inj,k_inj)
            cc(i_inj,k_inj,id_Si)  = cc(i_inj,k_inj,id_Si) &
                            +86400.0_rk*dt*rv_Si(i_day)/volume(i_inj,k_inj)
            cc(i_inj,k_inj,id_DON) = cc(i_inj,k_inj,id_DON) &
                            +86400.0_rk*dt*rv_DON(i_day)/volume(i_inj,k_inj)
            cc(i_inj,k_inj,id_Na)  = cc(i_inj,k_inj,id_Na) &
                            +86400.0_rk*dt*rv_Na(i_day)/volume(i_inj,k_inj)
            cc(i_inj,k_inj,id_Ca_u)  = cc(i_inj,k_inj,id_Ca_u) &
                            +86400.0_rk*dt*rv_Ca(i_day)/volume(i_inj,k_inj)
            cc(i_inj,k_inj,id_Mg)  = cc(i_inj,k_inj,id_Mg) &
                            +86400.0_rk*dt*rv_Mg(i_day)/volume(i_inj,k_inj)
            cc(i_inj,k_inj,id_K)  = cc(i_inj,k_inj,id_K) &
                            +86400.0_rk*dt*rv_K(i_day)/volume(i_inj,k_inj)
            cc(i_inj,k_inj,id_Cl) = cc(i_inj,k_inj,id_Cl) &
                            +86400.0_rk*dt*rv_Cl(i_day)/volume(i_inj,k_inj)
            cc(i_inj,k_inj,id_SO4) = cc(i_inj,k_inj,id_SO4) &
                            +86400.0_rk*dt*rv_SO4(i_day)/volume(i_inj,k_inj)
            cc(i_inj,k_inj,id_DIC) = cc(i_inj,k_inj,id_DIC) &
                            +86400.0_rk*dt*rv_DIC(i_day)/volume(i_inj,k_inj)
            cc(i_inj,k_inj,id_NH4) = cc(i_inj,k_inj,id_NH4) &
                            +86400.0_rk*dt*rv_NH4(i_day)/volume(i_inj,k_inj)
            cc(i_inj,k_inj,id_Phy) = cc(i_inj,k_inj,id_Phy) &
                            +86400.0_rk*dt*rv_Phy(i_day)/volume(i_inj,k_inj)
            cc(i_inj,k_inj,id_Alk) = cc(i_inj,k_inj,id_Alk) &
                            +86400.0_rk*dt*rv_Alk(i_day)/volume(i_inj,k_inj)
10 continue
!===================================================
    ! end of Subloop over timesteps in the course of one day     
    end do

!####################################################
!   if ((i_day/2).eq.((i_day+1)/2)) then
    !Save output to netcdf every day
    fick_per_day = 86400.0_rk * fick
    sink_per_day = 86400.0_rk * sink

    if (sediments_units_convert.eq.0) then 

        if (ncoutfile_type == 1) then
            call save_netcdf(i_max, k_max, max(1,julianday), cc, t, s, kz, kzti, wti, model, z, hz, Eair, use_Eair, hice, use_hice, &
            fick_per_day, sink_per_day, ip_sol, ip_par, x, u_x, i_day+1)
        else
            call save_netcdf(i_max, k_max, max(1,julianday), cc, t, s, kz, kzti, wti, model, z, hz, Eair, use_Eair, hice, use_hice, &
            fick_per_day, sink_per_day, ip_sol, ip_par, x, u_x, julianday)
        endif

    else
        !convert into observarional units in the sediments for dissolved (mass/pore water) and solids (mass/mass) with an exception for biota
        cc_out(:,:,:)=cc(:,:,:)
        do ip=1,par_max
            if (ip.ne.id_Phy.or.ip.ne.id_Het.or.ip.ne.id_Baae.or.ip.ne.id_Baan.or.ip.ne.id_Bhae.or.ip.ne.id_Bhan) then 
                cc_out(:,k_bbl_sed+1:k_max,:)=cc(:,k_bbl_sed+1:k_max,:)*pF1(:,k_bbl_sed:k_max-1,:)
            endif
        enddo 

        if (ncoutfile_type == 1) then
            call save_netcdf(i_max, k_max, max(1,julianday), cc_out, t, s, kz, kzti, wti, model, z, hz, Eair, use_Eair, hice, use_hice, &
            fick_per_day, sink_per_day, ip_sol, ip_par, x, u_x, i_day+1)
        else
            call save_netcdf(i_max, k_max, max(1,julianday), cc_out, t, s, kz, kzti, wti, model, z, hz, Eair, use_Eair, hice, use_hice, &
            fick_per_day, sink_per_day, ip_sol, ip_par, x, u_x, julianday)
        endif        

    endif

!   endif 
    
    !Save .dat files for plotting with Grapher for Sleipner for days 72 and 240
    !Note: saving to ascii every day causes an appreciable decrease in speed of execution
    if (julianday == 364) then
        vv=1._rk
        call saving_state_variables(trim(outfile_name), model_year, julianday, i_max, k_max, par_max, par_name, z, hz, cc, vv, t, s, kz)
    endif
    
    !Write daily output to screen if required
    do i=i_min,i_water
        if (show_maxmin.eq.1) then
            write(*,*) "maxval(cc(i,:,:),1) = ", maxval(cc(i,:,:),1) !Good to keep an eye on max/min values for each parameter
            write(*,*) "minval(cc(i,:,:),1) = ", minval(cc(i,:,:),1)
        end if
        if (show_kztCFL.gt.0) then
            do ip=1,par_max
                kztCFL(:,ip) = (kzti(i,2:k_max,ip)*86400.0_rk*dt/freq_turb)/(dz(1:k_max-1)**2)
            end do
            if (show_kztCFL.eq.2) then
                write(*,*) "maxval(kzti,2) = ", maxval(kzti,2)
                write(*,*) "maxval(kztCFL,2) = ", maxval(kztCFL,2)
            else
                write(*,*) "maxval(kzti) = ", maxval(kzti)
                write(*,*) "maxval(kztCFL) = ", maxval(kztCFL)
            end if
        end if
        if (show_wCFL.gt.0) then
            do ip=1,par_max
                wCFL(:,ip) = (abs(wti(i,2:k_max,ip))*86400.0_rk*dt/freq_sed)/dz(1:k_max-1)
            end do
            if (show_wCFL.eq.2) then
                write(*,*) "maxval(wti(i,:,:),2) = ", maxval(wti(i,:,:),2)
                write(*,*) "maxval(wCFL,2) = ", maxval(wCFL,2)
            else
                write(*,*) "maxval(wti) = ", maxval(wti)
                write(*,*) "maxval(wCFL) = ", maxval(wCFL)
            end if
        end if
    enddo
    
    ! end of BIG cycle    
    end do

    end subroutine do_brom_transport
!=======================================================================================================================







!=======================================================================================================================
    subroutine clear_brom_transport()

    deallocate(z_w)
    deallocate(dz_w)
    deallocate(hz_w)
    deallocate(t_w)
    deallocate(u_x_w)
    deallocate(s_w)
    deallocate(kz_w)
    deallocate(hmix_rate_w)
    deallocate(cc_hmix_w)
    deallocate(z)
    deallocate(dz)
    deallocate(hz)
    deallocate(t)
    deallocate(s)
    deallocate(u_x)
    deallocate(kz)
    deallocate(hmix_rate)
    deallocate(cc_hmix)
    deallocate(kz_bio)
    deallocate(kz_mol)
    deallocate(alpha)
    deallocate(phi)
    deallocate(phi1)
    deallocate(phi_inv)
    deallocate(tortuosity)
    deallocate(w_b)
    deallocate(u_b)
    deallocate(wti)
    deallocate(pF1)
    deallocate(pF2)
    deallocate(pWC)
    deallocate(cc)
    deallocate(cc_out)
    deallocate(dcc)
    deallocate(vv)
    deallocate(dvv)
    deallocate(fick)
    deallocate(fick_per_day)
    deallocate(sink)
    deallocate(sink_per_day)
    deallocate(wbio)
    deallocate(surf_flux)
    deallocate(bott_flux)
    deallocate(bc_top)
    deallocate(bc_bottom)
    deallocate(bctype_top)
    deallocate(bctype_bottom)
    deallocate(par_name)
    deallocate(Izt)
    deallocate(pressure)
    deallocate(Eair)
    deallocate(hice)
    deallocate(cc_top)
    deallocate(cc_bottom)
    deallocate(kzti)
    deallocate(kztCFL)
    deallocate(wCFL)
    deallocate(k_wat)
    deallocate(k_sed)
    deallocate(k_sed1)
    deallocate(x)
    deallocate(dx)
    deallocate(hx)
    if (diff_method.gt.0) call clean_tridiagonal()

    end subroutine clear_brom_transport
!=======================================================================================================================


    end module brom_transport
