! set date types for WECO based on TECO structure:
!   1. run attributes 
!   2. site based parameters
!   3. constant paramters
!   4. input and output varibles (file names)
!
! editor: Zhou Jian
! --------------------------------------------------------------------
module mod_dataTypes
    implicit none
    character(len=50) :: siteID = 'SPRUCE'             ! which site to run, for define the input and output file names
    character(len=50) :: namelistfile, strParas
    character(len=50) :: climfile = "input/forcing.txt"
    ! public :: run_nml, siteParas_nml, constParas_nml, outFiles_nml
    public :: run_nml, siteParas_nml, constParas_nml, initStates_nml
    ! ------------------------------------------------------------------
    ! contant parameters
    real, parameter :: pi = 3.14159256
    !-------------------------------------------------------------------
    ! parameters for simulation
    integer :: simuType ! which type to simulate: normal; spin-up; MCMC
    integer, public, save :: startYear, endYear, nyear
    real,    public, save :: Ttreat,    CO2treat
    logical :: doSoilphy = .True. 
    logical :: doSnow    = .True. 
    logical :: doMethane = .True.
    !-------------------------------------------------------------------
    ! parameters from parafile
    real lat,longi,wsmax,wsmin
    real LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax
    real SapR,SapS,SLA,GLmax,GRmax,Gsmax,stomN
    real a1,Ds0,Vcmax0,extkU,xfang,alpha
    real tauL,tauW,tauR,tauFL,tauCL
    real tauFS,tauSS,tauPS
    real gddonset,Q10,Rl0,Rs0,Rr0  
    real rMe,Q10pro,kCH4,Omax,ThreCH4,Tveg,TproMe,Toxi  ! for methane DA:parameters added for MEMCMC 
    !--------------------------------------------------------------------
    ! for consts parameteres
    real,dimension(3):: trnsmtL,rhoL,rhoS
    real emleaf,emsoil,Rconst,sigma,cpair,Patm,Trefk
    real H2OLv0,airMa,H2OMw,chi,Dheat,wleaf,gsw0,theta
    real conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,Edjm
    real Entrpy,gam0,gam1,gam2
    ! for soil thermal, snow DA  ..int
    real shcap_snow,condu_snow,condu_b,depth_ex,diff_s,diff_snow
    real albedo_snow,resht,thd_snow_depth,b_bound
    real infilt_rate 
    real fa,fsub,rho_snow,decay_m
    !--------------------------------------------------------------------
    ! for force data
    integer, public, save :: ilines ! = 4*366*24
    ! integer,dimension(nLines0):: yearSeq,doySeq,hourSeq
    ! real forceData(nLines0,rowForce) ! zhou: change to (line, row)
    ! real water_table(nLines0),snow_in(nLines0)
    !--------------------------------------------------------------------
    ! for initial state 
    real initSto, initNSC, initQC(8)
    real fwsoil, topfws, omega
    real CN0(8), thksl(10), FRLEN(10), CH4_V(10), CH4(10)
    real sftmp, Tsnow, Twater,Tice
    real G,snow_dsim,dcount,dcount_soil, ice_tw, zwt
    real Tsoill(10), ice(10), liq_water(10) 
    real N_deposit, N_fert 
    real infilt, accumulation, SNvcmax ! infilt_rate ???
    real alphaN, NSN, QNminer, N_deficit! initial values of Nitrogen pools and C/N ratio
    ! no namelist 
    real wcl(10), Esoil, water_tw, TauC(8)
    real WILTPT, FILDCP
    real stor_use, LAI, bmleaf, bmstem, bmroot,bmplant 
    real CN(8), QN(8), QNplant, QC(8), OutC(8), Rh_pools(5), OutN(8)
    real add,difference,gpp, npp, Rauto, Rgrowth, Rmain, Rnitrogen, store
    real RaLeaf, RaStem, RaRoot, RgLeaf, RgRoot, rgstem, rmleaf, rmroot, rmstem
    real pheno
    real storage
    real Hsoil
    real L_fall,alpha_L,alpha_W,alpha_R
    real N_uptake
    real tsoil_layer(11)
    !--------------------------------------------------------------------
    ! forcing variables for simulations
    real    :: Tair          ! air temperature,  degree
    real    :: Tsoil         ! soil temperature, degree
    real    :: RH            ! relative humidity
    real    :: VPD           
    real    :: rain          ! kgH2O m-2 s-1
    real    :: windU         ! wind velocity (m s-1)
    real    :: PAR           ! umol m-2 s-1
    real    :: radiation     ! W/m2
    real    :: soilwater     ! soil moisture, vol/vol
    real    :: snowDepth     !
    real    :: co2ca
    ! why?
    real TairK
    real eairP, Dair, radsol         ! Added for SPRUCE, due to lack of VPD data
    !-------------------------------------------------------------------
    ! maybe output?
    real ET, transp, evap, runoff
    real rain_yr, transp_yr, evap_yr, runoff_yr     ! yearly output
    ! -------------------------------------------------------------------
    ! simulation variable for driver process
    real fnsc, NSC
    real vcmx0, eJmx0 ! what difference with vcmax0? vcmax0 has default value. 
    real sps ! plantgrowth
    real snow_depth ! for soil_thermal
    real melt
    real wsc(10), phi ! from soil water to methane
    real SNRauto ! from growth to respiration, zhou: why not initial?
    ! -------------------------------------------------------------------
    type, public :: siteParaTypes
        real lat,longi,wsmax,wsmin
        real LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax
        real SapR,SapS,SLA,GLmax,GRmax,Gsmax,stomN
        real a1,Ds0,Vcmax0,extkU,xfang,alpha
        real tauL,tauW,tauR,tauFL,tauCL
        real tauFS,tauSS,tauPS
        real gddonset,Q10,Rl0,Rs0,Rr0  
        real rMe,Q10pro,kCH4,Omax,ThreCH4,Tveg,TproMe,Toxi  ! for methane DA:parameters added for MEMCMC 
    end type siteParaTypes

    type :: constParaTypes
        real,dimension(3):: trnsmtL,rhoL,rhoS
        real emleaf,emsoil,Rconst,sigma,cpair,Patm,Trefk
        real H2OLv0,airMa,H2OMw,chi,Dheat,wleaf,gsw0,theta
        real conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,Edjm
        real Entrpy,gam0,gam1,gam2
        real shcap_snow,condu_snow,condu_b,depth_ex,diff_s,diff_snow
        real albedo_snow,resht,thd_snow_depth,b_bound
        real infilt_rate
        real fa,fsub,rho_snow,decay_m
    end type constParaTypes

    ! type :: initStateTypes
    !     real initSto, initNSC, initQC,         &
    !     real fwsoil, topfws, omega,                                  &
    !     real CN0, thksl, FRLEN, CH4_V, CH4,                          &
    !     real sftmp, Tsnow, Twater,Tice,                              &
    !     real G,snow_dsim,dcount,dcount_soil, ice_tw, zwt,            &
    !     real Tsoill, ice, liq_water,                                 &
    !     real N_deposit, N_fert,                                      & 
    !     real infilt, accumulation, SNvcmax,                          & 
    !     real alphaN, NSN, QNminer, N_deficit 
    ! end type initStateTypes
    !--------------------------------------------------------------------
    ! force data
    type :: forceDataTypes !
        integer :: year          ! Year
        integer :: doy           ! day of the year
        real    :: hod           ! hour of the day
        real    :: Tair          ! air temperature,  degree
        real    :: Tsoil         ! soil temperature, degree
        real    :: RH            ! relative humidity
        real    :: VPD           
        real    :: rain          ! kgH2O m-2 s-1
        real    :: windU         ! wind velocity (m s-1)
        real    :: PAR           ! umol m-2 s-1
        real    :: radiation     ! W/m2
        real    :: soilwater     ! soil moisture, vol/vol
        real    :: snowDepth     !
        ! real    :: P_air         ! pa
        ! real    :: CO2           ! ppm
    end type forceDataTypes

    type :: fluxPoolTypes
        real vegCPool(3)  ! leaf; stem; root
        real soilCpool(5) ! fine litter; coarse litter; fast, slow, passive pools
    end type fluxPoolTypes

    ! type :: vegParasTypes
    !     real StemSap
    ! end type vegParasTypes

    ! type :: vegStateTypes
    !     real :: GPP
    !     real :: NPP
    !     real :: NEE
    !     real :: Re
    !     real :: Ra
    !     real :: Rh
    !     ! daily
    !     real :: dailyGPP
    !     real :: dailyNPP
    !     real :: dailyNEE
    !     real :: dailyRe
    !     real :: dailyRa
    !     real :: dailyRh
    !     ! yearly
    !     real :: annualGPP
    !     real :: annualNPP
    !     real :: annualNEE
    !     real :: annualRe
    !     real :: annualRa
    !     real :: annualRh
    ! end type vegStateTypes
    !====================================================================

    namelist /run_nml/ simuType, startYear, endYear, Ttreat, CO2treat

    namelist /siteParas_nml/                                    &
        lat,longi,wsmax,wsmin,                                  &   
        LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax,                   &
        SapR,SapS,SLA,GLmax,GRmax,Gsmax,stomN,                  &
        a1,Ds0,Vcmax0,extkU,xfang,alpha,                        &
        tauL,tauW,tauR,tauFL,tauCL,tauFS,tauSS,tauPS,           &
        gddonset,Q10,Rl0,Rs0,Rr0,                               &
        rMe,Q10pro,kCH4,Omax,ThreCH4,Tveg,TproMe,Toxi

    namelist /constParas_nml/                                   &
        trnsmtL,rhoL,rhoS,                                      &
        emleaf,emsoil,Rconst,sigma,cpair,Patm,Trefk,         &
        H2OLv0,airMa,H2OMw,chi,Dheat,wleaf,gsw0,theta,    &
        conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,Edjm,         &
        Entrpy,gam0,gam1,gam2,                                  &
        shcap_snow,condu_snow,condu_b,depth_ex,diff_s,diff_snow,&
        albedo_snow,      &
        resht,thd_snow_depth,b_bound,                           &
        infilt_rate,fa,fsub,rho_snow,decay_m
    
    namelist /initStates_nml/ initSto, initNSC, initQC,         &
        fwsoil, topfws, omega,                                  &
        CN0, thksl, FRLEN, CH4_V, CH4,                          &
        sftmp, Tsnow, Twater,Tice,                              &
        G,snow_dsim,dcount,dcount_soil, ice_tw, zwt,            &
        Tsoill, ice, liq_water,                                 &
        N_deposit, N_fert,                                      & 
        infilt, accumulation, SNvcmax,                          & 
        alphaN, NSN, QNminer, N_deficit                        
           
    ! type(siteParaTypes),  save :: sitePnDemo
    type(forceDataTypes), pointer, save :: forceData(:)
    type(siteParaTypes),  pointer, save :: sitePn
    type(constParaTypes), pointer, save :: constPn
    ! type(initStateTypes), pointer, save :: initStatePn
    type(fluxPoolTypes),  pointer, save :: fpVars
    ! type runParaTypes

    ! end type runParaTypes

    ! namelist /test_nml/ &
    !     testA

    contains    ! run initial processes or other update processes/diagnostics
    ! =============================================================================
    ! just one nml file
    subroutine init_parasNml(nmlName, strParas)
        implicit none
        character(len=50),intent(in) :: nmlName
        character(len=50),intent(in) :: strParas
        ! ---- local vars -----
        integer :: io           ! i/o status for the namelist
        integer :: nml_unit
        nml_unit = 999
        open(nml_unit, file=nmlName, form='formatted', action='read', status='old')
        select case(strParas)
            case('runVars')
                read(nml_unit, nml=run_nml,        iostat=io, end=10)
            case('siteVars')
                read(nml_unit, nml=siteParas_nml,  iostat=io, end=10)
            case('constVars')
                read(nml_unit, nml=constParas_nml, iostat=io, end=10)
            case('initStateVars')
                read(nml_unit, nml=initStates_nml, iostat=io, end=10)
            case default
                print *, "invalid nml name. you must choose: runVars, siteVars, constVars OR initStateVars!"
                stop strParas//"Please check strParas..."
        end select
    10  close (nml_unit)
    end subroutine init_parasNml
    ! ------------------------------------------------------------------------------------------
    ! initialize the site based parameters
    subroutine init_sitePn()!siteParas)
        implicit none
        ! type(siteParaTypes),pointer,intent(inout) :: siteParas
        type(siteParaTypes), pointer :: siteData
        allocate(siteData)  
        siteData= siteParaTypes(& 
            lat,longi,wsmax,wsmin,                                  &   
            LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax,                   &
            SapR,SapS,SLA,GLmax,GRmax,Gsmax,stomN,                  &
            a1,Ds0,Vcmax0,extkU,xfang,alpha,                        &
            tauL,tauW,tauR,tauFL,tauCL,tauFS,tauSS,tauPS,           &
            gddonset,Q10,Rl0,Rs0,Rr0,                               &
            rMe,Q10pro,kCH4,Omax,ThreCH4,Tveg,TproMe,Toxi) 
        ! siteParas%lat = lat
        sitePn => siteData
    end subroutine init_sitePn
    !--------------------------------------------------------------------------------------------
    ! initialize the constant parameters
    subroutine init_constPn() !constParas)
        implicit none
        ! type(constParaTypes),pointer,intent(inout) :: constParas
        type(constParaTypes), pointer :: constData
        allocate(constData)  
        constData= constParaTypes(& 
            trnsmtL,rhoL,rhoS,                                      &
            emleaf,emsoil,Rconst,sigma,cpair,Patm,Trefk,         &
            H2OLv0,AirMa,H2OMw,chi,Dheat,wleaf,gsw0,theta,    &
            conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,Edjm,         &
            Entrpy,gam0,gam1,gam2,                                  &
            shcap_snow,condu_snow,condu_b,depth_ex,diff_s,diff_snow, albedo_snow,      &
            resht,thd_snow_depth,b_bound,                           &
            infilt_rate,fa,fsub,rho_snow,decay_m) 
        ! siteParas%lat = lat
        constPn => constData
    end subroutine init_constPn
    !--------------------------------------------------------------------------------------------
    ! ! initialize the states, include the C pools
    ! subroutine init_states()
    !     implicit none
    !     type(constParaTypes), pointer :: constData
    !     allocate(constData)  
    !     constData= constParaTypes(& 
    !         trnsmtL,rhoL,rhoS,                                      &
    !         emleaf,emsoil,Rconst,sigma,cpair,Patm,Trefk,         &
    !         H2OLv0,AirMa,H2OMw,chi,Dheat,wleaf,gsw0,theta,    &
    !         conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,Edjm,         &
    !         Entrpy,gam0,gam1,gam2,                                  &
    !         shcap_snow,condu_snow,condu_b,depth_ex,diff_s,diff_snow, albedo_snow,      &
    !         resht,thd_snow_depth,b_bound,                           &
    !         infilt_rate,fa,fsub,rho_snow,decay_m) 
    !     ! siteParas%lat = lat
    !     constPn => constData
    ! end subroutine init_states
    !--------------------------------------------------------------------------------------------
    ! to read forcing data
    subroutine readForceData(forcingData,datalines,days_data,yr_data,timestep)
        implicit none ! revised from Biomass Model by Weng laoshi
        type(forceDataTypes),pointer,intent(inout) :: forcingData(:)
        integer,intent(inout) :: datalines,days_data,yr_data
        real, intent(inout)   :: timestep
        !------------local var -------------------
        type(forceDataTypes), pointer :: climateData(:)
        character(len=80)  commts
        integer, parameter :: niterms = 9       ! MDK data for Oak Ridge input
        ! integer nyear
        ! integer, parameter :: ilines = 4*366*24 ! the maxmum records of forcing data
        integer,dimension(ilines) :: year_data
        real,   dimension(ilines) :: doy_data,hour_data
        real :: input_data(niterms,ilines)
        ! real inputstep
        integer :: istat1,istat2,istat3
        integer :: doy,idays
        integer :: i,j,k
        integer :: m,n
        write(*,*) "zhoujian: ilines", ilines
        ! climfile=trim(filepath_in)//trim(climfile)
        ! open forcing data
        open(11,file=climfile,status='old',ACTION='read',IOSTAT=istat2)
        write(*,*)istat2
        ! skip 2 lines of input met data file
        read(11,'(a160)') commts
        ! read(11,'(a160)') commts ! MDK data only has one line comments
        m       = 0  ! to record the lines in a file
        idays   = 1  ! the total days in a data file
        yr_data = 0  ! to record years of a dataset
        do    ! read forcing files
            m=m+1
            read(11,*,IOSTAT=istat3)year_data(m),doy_data(m),hour_data(m),   &
                                    (input_data(n,m),n=1,niterms)
            if(istat3<0)exit
            if(m == 1) then
                doy = doy_data(m)
            else
                doy = doy_data(m-1)
            endif
            if(doy /= doy_data(m)) idays = idays + 1
            ! write(*,*)year_data(m),doy_data(m),hour_data(m)!, input_data(:,m)
        enddo ! end of reading the forcing file

        timestep = hour_data(2) - hour_data(1)
        ! write(*,*)"forcing",datalines,yr_data,timestep !,dt_fast_yr
        if (timestep==1.0)then
            write(*,*)"the data freqency is hourly"
        elseif(timestep==0.5)then
            write(*,*)"the data freqency is half hourly"
        else
            write(*,*)"Please check time step!"
            stop
        endif
        close(11)    ! close forcing file
        ! Put the data into forcing 
        datalines = m - 1
        days_data = idays
        yr_data   = year_data(datalines-1) - year_data(1) + 1

        allocate(climateData(datalines))
        do i=1,datalines
            climateData(i)%year      = year_data(i)          ! Year
            climateData(i)%doy       = doy_data(i)           ! day of the year
            climateData(i)%hod       = hour_data(i)          ! hour of the day
            climateData(i)%Tair      = input_data(1,i) !+ 273.16  ! air temperature, K
            climateData(i)%Tsoil     = input_data(2,i) !+ 273.16  ! soil temperature, K
            climateData(i)%RH        = input_data(3,i) !* 0.01    ! relative humidity (0.xx)
            climateData(i)%VPD       = input_data(4,i) 
            climateData(i)%rain      = input_data(5,i)  !/(timestep * 3600)! ! kgH2O m-2 s-1
            climateData(i)%windU     = input_data(6,i)        ! wind velocity (m s-1)
            climateData(i)%PAR       = input_data(7,i)       ! umol/m2/s
            climateData(i)%radiation = input_data(7,i)       ! W/m2
            climateData(i)%soilwater = input_data(8,i)!0.8    ! soil moisture, vol/vol
            climateData(i)%snowDepth = input_data(9,i)
            ! climateData(i)%P_air     = input_data(8,i)        ! pa
            ! climateData(i)%CO2       = input_data(9,i) * 1.0e-6       ! mol/mol
        enddo
        forcingData => climateData
        write(*,*)"forcing", datalines,days_data,yr_data
        return
    end subroutine readForceData
    !----------------------------------------------------------------------------------
    ! subroutine writeDailyOutputs()
    !     implicit none

    ! end 
    ! ----------------------------------------------------------------------------------
end module mod_dataTypes