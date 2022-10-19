module mod_driver
    use mod_dataTypes
    use mod_vegetation
    use mod_soil
    use mod_transferCpools
    use mod_methane
    implicit none

    contains
    !-----------------------------------------------------
    subroutine driver(forceData, sitePn, constPn, statePn)
        ! character(len=50) :: strParas
        type(forceDataTypes),pointer,intent(inout) :: forceData(:)
        type(siteParaTypes), pointer,intent(inout) :: sitePn
        type(constParaTypes),pointer,intent(inout) :: constPn
        type(fluxPoolTypes), intent(inout) :: statePn
        ! ------- local var --------- !
        integer iyr, iday, lenDoy, ihour, iforce, dlayer
        integer, parameter :: lenHod = 24 ! hourly simulation
        ! ------- 
        real GDD5, accumTa      !Accumulated temperature => Ta
        integer onset, phenoset, day_mod
        ! ------- maybe local ----------
        real StemSap, RootSap, NSCmin, NSCmax, rain_d ! not sure rain_d be used in other place
        real snow_depth_e
        real gpp_d
        ! write(*,*)"This is driver module ...", sitePn%lat, constPn%trnsmtL
        ! write(*,*)"test-Ttreat2:", QNplant
        ! write(*,*) "zhou:", statePn%soilCpool
        ! write(*,*)"test-NSN:", nsn
        ! stop
        ! ===============================================================
        accumTa = 0.  
        iforce  = 0
        do iyr = 1, nyear   !nyear
            if(MOD(startYear+iyr-1,4).eq.0)then  !! leap year
                lenDoy = 366
            else
                lenDoy = 365
            endif
            GDD5       = 0.0    ! initial phenoset variables
            onset      = 0
            phenoset   = 0
            ! lenDoy     = 101 ! zhou: for test
            !-----------------------------------------
            do iday = 1, lenDoy
                ! zhou :  no Nitrogen fertilization 
                ! if(days==75 .OR. days==105)then  ! Nitrogen fertilization since 1999 in Duke
                !     QNminer = QNminer + N_fert   !(5.6 gN/yr/m2,N fertiliztion in March and Apr)
                ! endif
                StemSap = AMIN1(Stemmax, SapS*bmStem)   ! Stemmax and SapS were input from parameter file, what are they? Unit? Maximum stem biomass? -JJJJJJJJJJJJJJJJJJJJJJ 
                RootSap = AMIN1(Rootmax, SapR*bmRoot)
                NSCmin  = 5. 
                NSCmax  = 0.05*(StemSap + RootSap + statePn%vegCPool(1))
                if(accumTa.gt.5.0) GDD5 = GDD5 + accumTa
                gpp_d = 0.
                ! ********* for daily initials in soil thermal module
                ! soilt_d_simu = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./) 
                ! ice_d_simu   = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0./) 
                ! soilt_d_obs  = (/0.,0.,0.,0.,0.,0.,0./) 
                ! zwt_d        = 0.0
                ! obs_counter  = (/0,0,0,0,0,0,0/)
                if (doSnow) then 
                    if (iyr .eq. 1. .and. iday .eq. 1.) then ! zhou: ta -> accumTa, why change Ta to -12.85?
                        accumTa = -12.85    !ta     = -12.85    ! since changed the ta criteria (0. to 1.e-10)) in calculating melt
                        rain_d = 0.        !dbmemo
                    endif
                    call snow_daily(rain_d,lat,iday,accumTa,snow_dsim,fa,fsub,rho_snow,melt,dcount,decay_m)                            
                    snow_depth_e = snow_dsim  !  snow_dsim, fa, fsub, rho_snow are const parameter, melt is output
                endif 
                do ihour = 1, lenHod
                    ! forcing data
                    iforce      = iforce + 1
                    Tair        = forceData(iforce)%Tair          ! air temperature,  degree
                    Tsoil       = forceData(iforce)%Tsoil         ! soil temperature, degree
                    RH          = forceData(iforce)%RH            ! relative humidity
                    VPD         = forceData(iforce)%VPD  
                    rain        = forceData(iforce)%rain          ! kgH2O m-2 s-1
                    windU       = forceData(iforce)%windU         ! wind velocity (m s-1)
                    PAR         = forceData(iforce)%PAR           ! umol m-2 s-1
                    radiation   = forceData(iforce)%radiation     ! W/m2 = >radsol
                    soilwater   = forceData(iforce)%soilwater     ! soil moisture, vol/vol
                    snowDepth   = forceData(iforce)%snowDepth     !
                    co2ca       = 380.0*1.0E-6
                    Tairk       = Tair + 273.2
                    ! --------------------------------------------------------------------
                    ! *** int added for soil thermal/ soil water
                    day_mod=mod(iforce,24) ! m -> iforce
                    if (doSnow) then
                        snow_depth = snow_depth_e
                    else
                        snow_depth = 0 !snow_in(m) ! zhou: just for test 20220211
                    endif
                    if (snow_depth .lt. 0.0) snow_depth = 0.0   
                    snow_depth = snow_depth*100.   ! change from m to cm    
                    ! Ajust some unreasonable values
                    RH     = AMAX1(0.01,AMIN1(99.99,RH))
                    eairP  = esat(Tair)*RH/100.             ! Added for SPRUCE, due to lack of VPD data
                    Dair   = esat(Tair)-eairP               
                    radsol = AMAX1(radiation,0.01)
                    ! *** int if do soil thermal G is not given a value here
                    !             else G will be given a value below
                    if (doSoilphy) then 
                        GOTO 160
                    endif
                    if(radsol.gt.10.0) then
                        G=-25.0
                    else
                        G=20.5											
                    endif
                    Esoil = 0.05*radsol               ! energy of soil?
                    if(radsol.LE.10.0) Esoil=0.5*G
                160 continue
                    accumTa = accumTa + Tair/24.0             ! sum of a day, for calculating daily mean temperature
                    ! calculating scaling factor of NSC
                    if(NSC.le.NSCmin) fnsc = 0.0
                    if(NSC.ge.NSCmax) fnsc = 1.0
                    if((NSC.lt.NSCmax).and.(NSC.gt.NSCmin))then 
                        fnsc = (NSC-NSCmin)/(NSCmax-NSCmin)
                    endif
                    Vcmx0 = Vcmax0*SNvcmax*1.0e-6   ! update vcmx0 and eJmx0 according to C/N of leaves
                    eJmx0 = 1.67*Vcmx0 ! ! eJmx0 = 2.7*Vcmx0  ! original; Weng 02/21/2011 Medlyn et al. 2002
                    call veg_canopy(iday, ihour, gpp)
                    ! write(*,*)"test-gpp:",iday, ihour,gpp, Rmain, LAI
                    ! write(*,*)"test-Rmain", Rmain, Rmleaf,Rmstem, Rmroot
                    ! write(*,*)"test-LAI", bmleaf, SLA
                    ! write(*,*)"test-QC:", L_fall
                    if(gpp >0)then
                        stop
                    endif
                    ! write(*,*)"test:", evap, transp, rain, evap, runoff
                    call soil_water(transp, evap, liq_water, ice, accumTa, &
                                    & wsc, phi)
                    ET        = evap      + transp
                    rain_yr   = rain_yr   + rain
                    transp_yr = transp_yr + transp
                    evap_yr   = evap_yr   + evap
                    runoff_yr = runoff_yr + runoff
                    call respiration(GPP, StemSap, RootSap, SNRauto, bmleaf, &
                                    & Rmleaf, rmstem, rmroot, Rmain) ! SNRauto update in growth module
                    call veg_growth(GDD5,NSCmin,NSCmax, onset, &
                        & StemSap, RootSap, Rgrowth,Rgroot,Rgleaf,Rgstem, L_fall,alpha_L,alpha_W,alpha_R, SNRauto)
                    call TCS_CN()
                    if(doMethane) call methane()
                    ! update NSC
                    Rauto      = Rmain+Rgrowth+Rnitrogen
                    NSC        = NSC+GPP-Rauto-(NPP-add)-store
                    Difference = GPP-Rauto-NPP
                    if(NSC<0)then
                        bmstem = bmstem+NSC/0.48
                        NPP    = NPP+NSC
                        NSN    = NSN-NSC/CN(2)
                        NSC    = 0.
                    endif
                    ! write(*,*)"test-NSC:", NSC, gpp, Rauto, NPP, add, store
                    ! write(*,*)"test-Rauto:", Rmain, Rgrowth, Rnitrogen
                    ! update daily states
                    ! GL_d    = GL_d+NPP*alpha_L
                    ! GW_d    = GW_d+NPP*alpha_W
                    ! GR_d    = GR_d+NPP*alpha_R
                    ! LFALL_d = LFALL_d+L_fall
                    ! update some variables
                    RaLeaf  = RgLeaf + RmLeaf
                    RaStem  = RgStem + RmStem
                    RaRoot  = RgRoot + RmRoot + Rnitrogen
                    ! WFALL_d = WFALL_d+OutC(2) !_wood
                    ! RFALL_d = RFALL_d+OutC(3) !_root
                    ! N_LG_d  = N_LG_d+N_leaf
                    ! N_WG_d  = N_WG_d+N_wood
                    ! N_RG_d  = N_RG_d+N_root
                    ! N_LF_d  = N_LF_d+N_LF
                    ! N_WF_d  = N_WF_d+N_WF
                    ! N_RF_d  = N_RF_d+N_RF

                    ! N_up_d    = N_up_d+N_uptake
                    ! N_fix_d   = N_fix_d+N_fixation
                    ! N_dep_d   = N_dep_d+N_deposit
                    ! N_leach_d = N_leach_d+N_leach
                    ! N_vol_d   = N_vol_d+N_vol

                    ! N_up_yr    = N_up_yr+N_uptake
                    ! N_fix_yr   = N_fix_yr+N_fixation
                    ! N_dep_yr   = N_dep_yr+N_deposit
                    ! N_leach_yr = N_leach_yr+N_leach
                    ! N_vol_yr   = N_vol_yr+N_vol

                    ! R_Ntr_yr   = R_Ntr_yr + Rnitrogen
                    ! initial some variables
                    ! do dlayer=1,10
                    !     ice_d_simu(dlayer)=ice_d_simu(dlayer)+ice(dlayer) 
                    ! enddo       
                    ! do dlayer=1,11
                    !     soilt_d_simu(dlayer)=soilt_d_simu(dlayer)+testout(dlayer)  
                    !     ! first = surface soil temperature 2:11=1:10 layer soil temperatures 
                    ! enddo                    
                    ! do dlayer=1,10
                    !     CH4V_d(dlayer)=CH4V_d(dlayer)+CH4_V(dlayer) 
                    ! enddo                  
                    ! zwt_d=zwt_d+zwt    ! ..int I doubt it... mean for zwt?     check later  Shuang 
                    ! Rhetero=Rh_f + Rh_c + Rh_Micr + Rh_Slow + Rh_Pass
                    ! Rhetero = Rh_pools(1)+Rh_pools(2)+Rh_pools(3) &
                    !     &     +Rh_pools(4)+Rh_pools(5)
                    ! Rsoil   = Rhetero+RmRoot+RgRoot+Rnitrogen
                    ! NEE     = Rauto+Rhetero - GPP
                    ! Q_soil  = QC(6) + QC(7) + QC(8)

                    bmleaf  = QC(1)/0.48
                    bmstem  = QC(2)/0.48
                    bmroot  = QC(3)/0.48
                    bmplant = bmleaf+bmroot+bmstem
                    LAI     = bmleaf*SLA
                    ! LAI     = 5.8 ! FOR TEST
                    ! write(*,*)"test-QC:", QC
                    ! NMIN_d  = NMIN_d+N_miner
                    ! ! output hourly
                    ! Recoh   = Rhetero+Rauto
                    ! ETh     = ET     !*1000.
                    ! Th      = transp !*1000.
                    ! Eh      = evap   !*1000.
                    ! INTh    = -9999
                    ! VPDh    = Dair/1000.
                    ! ROh     = runoff !*1000.
                    ! DRAINh  = -9999
                    ! LEh     = ETh*((2.501-0.00236*Tair)*1000.0)/3600.
                    ! SHh     = -9999
                    ! LWh     = -9999
                    ! NEP     = -NEE

                    ! ! sums of a day
                    ! diff_d   = diff_d+difference
                    gpp_d    = gpp_d + GPP
                    ! gpp_ra   = gpp_ra+Rauto
                    ! NPP_d    = NPP_d+NPP
                    ! NEP_d    = NEP_d+NEP
                    ! NEE_d    = NEE_d+NEE
                    ! RECO_d   = RECO_d+Recoh
                    ! Rh_d     = Rh_d + Rhetero
                    ! Ra_d     = Reco_d-Rh_d
                    ! RLEAV_d  = RLEAV_d+RmLeaf+RgLeaf
                    ! RWOOD_d  = RWOOD_d+RmStem+RgStem
                    ! RROOT_d  = RROOT_d+RmRoot+RgRoot+Rnitrogen
                    ! Rsoil_d  = Rh_d+RROOT_d
                    ! NUP_d    = NUP_d+N_uptake
                    ! NVOL_d   = NVOL_d+N_vol
                    ! NLEACH_d = NLEACH_d+N_leach
                    ! transp_d = transp_d + transp*(24./dtimes)
                    ! evap_d   = evap_d + evap*(24./dtimes)
                    ! ET_d     = transp_d + evap_d
                    ! LE_d     = LE_d+LEh/24.
                    ! Hcanop_d = Hcanop_d+Hcanop/(24./dtimes)
                    ! runoff_d = runoff_d+runoff
                    ! !   *** .int
                    ! ! added for MEMCMC also for generation of daily methane emission                  
                    ! simuCH4_d = simuCH4_d+simuCH4
                    ! Pro_sum_d = Pro_sum_d+Pro_sum
                    ! Oxi_sum_d = Oxi_sum_d+Oxi_sum
                    ! Fdifu1_d  = Fdifu1_d+Fdifu(1)
                    ! Ebu_sum_d = Ebu_sum_d+Ebu_sum
                    ! Pla_sum_d = Pla_sum_d+Pla_sum                                              
                    ! !   ***                  
                    ! ! sum of the whole year
                    ! diff_yr = diff_yr+difference
                    ! gpp_yr  = gpp_yr+gpp
                    ! NPP_yr  = NPP_yr+NPP
                    ! Rh_yr   = Rh_yr +Rhetero
                    ! Ra_yr   = Ra_yr+Rauto
                    ! Rh4_yr  = Rh4_yr+Rh_pools(1)
                    ! Rh5_yr  = Rh5_yr+Rh_pools(2)
                    ! Rh6_yr  = Rh6_yr+Rh_pools(3)
                    ! Rh7_yr  = Rh7_yr+Rh_pools(4)
                    ! Rh8_yr  = Rh8_yr+Rh_pools(5)
                    ! Pool1   = Pool1+QC(1)/8760.
                    ! Pool2   = Pool2+QC(2)/8760.
                    ! Pool3   = Pool3+QC(3)/8760.
                    ! Pool4   = Pool4+QC(4)/8760.
                    ! Pool5   = Pool5+QC(5)/8760.
                    ! Pool6   = Pool6+QC(6)/8760.
                    ! Pool7   = Pool7+QC(7)/8760.
                    ! Pool8   = Pool8+QC(8)/8760.
                    ! out1_yr = out1_yr+OutC(1)
                    ! out2_yr = out2_yr+OutC(2)
                    ! out3_yr = out3_yr+OutC(3)
                    ! out4_yr = out4_yr+OutC(4)
                    ! out5_yr = out5_yr+OutC(5)
                    ! out6_yr = out6_yr+OutC(6)
                    ! out7_yr = out7_yr+OutC(7)
                    ! out8_yr = out8_yr+OutC(8)
                    ! NEE_yr  = NEE_yr+NEE
                    ! GL_yr   = GL_yr+NPP*alpha_L
                    ! GW_yr   = GW_yr+NPP*alpha_W
                    ! GR_yr   = GR_yr+NPP*alpha_R
                    ! ! numbering         
                    ! n=n+1
                    ! !   *** .int
                    ! ! added for soil thermal      unknown function check later   Shuang
                    ! ! if((yr+first_year-1).eq.obs_soilwater(1,k1) .and.    &
                    ! !     &   days .eq. obs_soilwater(2,k1) .and.      &
                    ! !     &   (i-1).eq. obs_soilwater(3,k1))then
                    ! !     Simu_soilwater(1:10,k1)=wcl(1:10)
                    ! !     Simu_soiltemp(1:11,k1) =testout
                    ! !     Simu_watertable(1,k1)  =zwt
                    ! !     k1=k1+1
                    ! ! endif

                    
                    ! !   ***                  
                    if(isnan(gpp))then
                        write(*,*)'gpp is nan'
                        return
                    endif
                enddo ! end of all hour per day
                
                if((GDD5.gt.gddonset) .and. phenoset.eq.0) then
                    pheno    = iday
                    phenoset = 1
                endif
                ! write(*,*)"test-LAI:", iday, LAI
                write(*,*)"test:", iday, gpp_d !,npp
            enddo ! end of all day per year
            ! write(*,*)"test-gpp:", gpp, evap, Esoil
            ! stop "this is here ..."
            initSto      = accumulation
            stor_use     = storage/720. ! times_storage_use = 720
            accumulation = 0.0
            onset        = 0
        enddo  ! end of all year
    end subroutine driver
    ! -------------------------------------------------------------------------------
    ! real function esat(T)
    !     real T
    !     ! returns saturation vapour pressure in Pa
    !     esat = 610.78*exp(17.27*T/(T+237.3))
    ! return
    ! end

end module mod_driver