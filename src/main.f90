! This is the main module of WECO, include:
!   1. read preset parameters
!   2. read input (force, parameters) and output (file name)
!   3. initial data 
!   4. choose which process to run: normal; Spin-UP; MCMC
!  
! Editor: Zhou Jian
! --------------------------------------------------------------

program WECO
    use mod_dataTypes
    use mod_driver
    ! ----------------------------------
    integer datalines,days_data,yr_data
    real    timestep
    type(fluxPoolTypes), pointer :: statePn 
    ! ---------------------------------------------------------------------
    ! namelistFile = 'parameters_default.nml'
    namelistfile = 'parameters_'//trim(siteID)//'.nml'                     ! set parameters file
    strParas     = "runVars";   call init_parasNml(namelistfile, strParas) ! site based parameters
    strParas     = "siteVars";  call init_parasNml(namelistfile, strParas) ! const parameters
    strParas     = "constVars"; call init_parasNml(namelistfile, strParas) ! run setup parameters
    ! force data
    nyear  = endyear - startYear + 1
    ilines = nyear * 366 * 24
    call readForceData(forceData,datalines,days_data,yr_data,timestep)
    write(*,*)"zhou: force:", forceData(1)%Tair
    call init_sitePn()!sitePn)
    call init_constPn()
    write(*,*)"zhou: const:",diff_s
    ! initial state values
    strParas = "initStateVars";   call init_parasNml(namelistfile, strParas) ! site based parameters
    ! write(*,*) "zhou:", initSto
    ! write(*,*)"test-Ttreat:", Ttreat
    ! Ttreat = 10000.
    ! strParas = "initStateVars";   call init_parasNml(namelistfile, strParas) ! site based parameters
    ! write(*,*)"test-Ttreat:", Ttreat
    do i=1,10
        wcl(i) = wsmax/100.
    enddo 
    water_tw = zwt*0.001
    tauC     = (/tauL,tauW,tauR,tauFL,tauCL,tauFS,tauSS,tauPS/)*8760. ! the unit of residence time is transformed from yearly to hourly
    SLA      = SLA/10000.         ! Convert unit from cm2/g to m2/g
    ! growth rates of plant
    GLmax    = GLmax/8760.
    GRmax    = GRmax/8760.
    Gsmax    = Gsmax/8760.
    ! FOR what ?
    WILTPT   = wsmin/100.0
    FILDCP   = wsmax/100.0
    ! 
    stor_use = initSto/720.       ! times_storage_use:  720 hours, 30 days
    NSC      = initNSC
    LAI      = LAIMIN
    bmleaf   = initQC(1)/0.48
    bmstem   = initQC(2)/0.48
    bmroot   = initQC(3)/0.48
    bmplant  = bmstem+bmroot+bmleaf
    QC       = initQC
    !---------------------------
    CN       = CN0
    QN       = initQC/CN0
    QNplant  = QN(1) + QN(2) + QN(3)
    ! -----------------------------
    allocate(statePn)
    statePn%vegCPool(:)  = initQC(1:3)
    statePn%soilCpool(:) = initQC(4:8)
    ! ======================================================
    select case(simuType)
        case (1)
            write(*,*)"normal simulation ...."
            call driver(forceData, sitePn, constPn, statePn)
        case (2)
            write(*,*)"spin-up"
        case (3)
            write(*,*)"MCMC"
    end select
    ! contains
    ! subroutine 
    deallocate(statePn)
end program WECO