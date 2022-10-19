!
!
module mod_vegetation
    ! USE mod_driver
    use mod_dataTypes
    USE mod_soil
    implicit none
    ! logical :: doSoilphy = .False.
    !-----------------------------------------------------------
    contains
    subroutine veg_canopy(idoy, ihod, gpp)
    ! This subroutine to calculate the canopy variables, including GPP, evop.
    ! inputs:       idoy, ihod          --> idoy: day of year to calculate radiation.
    ! ouputs:       gpp, evap, transp   --> carbon, water and energy variables
    ! mid-ouputs:   Rsoilabs            --> to soil thermal module: total radiation absorbed by soil
    ! ------------------------------------------------------------------------------------------------------
        implicit none
        ! ---- input vars  ---- !
        integer idoy, ihod
        ! ---- output vars ---- !
        real gpp, transp, evap
        ! ---- local vars  ---- !
        ! -----------------------------------------------
        real fbeam                          ! output from yrday module: beam fraction in incoming solar radiation
        real coszen                         ! from cal_sinbet
        real FLAIT, Radabv(2)
        real Acan1, Acan2, Ecan1, Ecan2     ! from xlayers
        real Rsoilab1, Rsoilab2, QLleaf, QLair, raero ! from xlayers to soil_thermal
        real Acanop, Ecanop
        ! -------------------------------------------------------------------
        call  yrday(idoy, ihod, fbeam)     ! calculate beam fraction in incoming solar radiation
        coszen = cal_sinbet(idoy,ihod)     ! cos zenith angle of sun
        if(windU.lt.0.01) windU = 0.01      ! set windspeed to the minimum speed to avoid zero Gb
        if(topfws.gt.0.5) then              ! calculate soil albedo for NIR as a function of soil water (Garratt pp292)
            rhoS(2) = 0.18                  ! zhou: rhoS is a const parameter, and it is updated?
        else
            rhoS(2) = 0.52-0.68*topfws
        endif
        
        ! assign plant biomass and leaf area index at time t, assume leaf biomass = root biomass
        FLAIT     = LAI                        ! initial LAI = LAIMIN
        radabv(1) = 0.5*radsol                 !(1) - solar radn
        radabv(2) = 0.5*radsol                 !(2) - NIR
        ! call multilayer model of Leuning - uses Gaussian integration but radiation scheme is that of Goudriaan
        call xlayers(coszen, FLAIT, fbeam, Radabv, &
                    &   Rsoilab1, Rsoilab2, QLleaf, QLair, raero, & ! outputs for soil_thermal
                    &   Acan1, Acan2, Ecan1, Ecan2, Esoil, Hsoil) ! outputs 
        ! write(*,*)"before Esoil:", Esoil
        if (doSoilphy) then 
            call soil_thermal(Rsoilab1, Rsoilab2, QLleaf, QLair, FLAIT, raero, &
                              & Esoil, Hsoil)
        endif
        ! write(*,*)"after Esoil:", Esoil
        Acanop = Acan1+Acan2
        Ecanop = Ecan1+Ecan2
        gpp    = Acanop*3600.0*12.0                                   ! every hour, g C m-2 h-1
        transp = AMAX1(Ecanop*3600.0/(1.0e6*(2.501-0.00236*Tair)),0.) ! mm H2O /hour
        evap   = AMAX1(Esoil*3600.0/(1.0e6*(2.501-0.00236*Tair)),0.)
        ! write(*,*)"test-gpp:", gpp, evap, Esoil
        ! write(*,*)"test-Acanop", Acanop, Acan1, Acan2
    end subroutine veg_canopy

    ! -----------------------------------------------------------
    subroutine respiration(GPP, StemSap, RootSap, SNRauto,bmleaf, &
                            & Rmleaf, rmstem, rmroot, Rmain) ! outputs
        !   autotrophic respiration
        !   calculate plant and soil respiration by the following equation:
        !   RD=BM*Rd*Q10**((T-25)/10) (Sun et al. 2005. Acta Ecologica Sinica)
        implicit none
        real GPP
        real StemSap,RootSap
        real SNRauto                ! update in transfer  
        real bmleaf !, bmstem,bmroot  ! no used bmstem and bmroot?
        real RmLeaf,RmStem,RmRoot,Rmain
        real conv                  ! converter from "umol C /m2/s" to "gC/m2/hour"
        conv = 3600.*12./1000000.    ! umol C /m2/s--> gC/m2/hour
        if(LAI.gt.LAIMIN) then
            RmLeaf = Rl0*SNRauto*bmleaf*0.48*SLA*0.1     &               
                        &   *Q10**((Tair-10.)/10.)*fnsc*conv
            RmStem = Rs0*SNRauto*StemSap*0.001*Q10**((Tair-25.)/10.)*fnsc*conv
            RmRoot = Rr0*SNRauto*RootSap*0.001*Q10**((Tair-25.)/10.)*fnsc*conv
        else
        RmLeaf = 0.3*GPP
        RmStem = 0.3*GPP
        RmRoot = 0.4*GPP
        endif
        Rmain  = Rmleaf+Rmstem+Rmroot
        if(Rmain > 0.0015*NSC)then             ! If Total autotropic respiration greater than 0.15% of Nonstructure Carbon, rescale. 
            Rmleaf = Rmleaf/Rmain*0.0015*NSC
            Rmstem = Rmstem/Rmain*0.0015*NSC
            Rmroot = Rmstem/Rmain*0.0015*NSC
            Rmain  = Rmleaf+Rmstem+Rmroot
        endif
        ! write(*,*)"test-RmStem:", Rs0, SNRauto, StemSap, Q10, ((Tair-25.)/10.), fnsc, conv
        ! write(*,*)"test-RmRoot", Rr0, SNRauto, RootSap, Q10, ((Tair-25.)/10.), fnsc, conv
        ! write(*,*)"test-RmLeaf:", Rl0, SNRauto,bmleaf,SLA, Q10, ((Tair-10.)/10.), fnsc, conv
        return
    end subroutine respiration

    ! ----------------------------------------------------------------------------------
    subroutine veg_growth(GDD5,NSCmin,NSCmax, onset, &
                        & StemSap, RootSap, Rgrowth,Rgroot,Rgleaf,Rgstem, L_fall,alpha_L,alpha_W,alpha_R, SNRauto)
        ! plant growth model
        implicit none
        integer onset
        real GDD5,NSCmin,NSCmax  ! from driver, put it to dataTypes?
        real SNRauto
        ! ----- local vars ----- !
        real bmL,bmR,bmP,bmS !,StemSap,RootSap
        ! real store
        real GrowthP,GrowthL,GrowthR,GrowthS
        real ht,hmax,hl0,CNP0
        REAL LAIMAX0,la0,GPmax,acP,c1,c2 ! this must be define in site based parameters
        real,save :: addaccu=0,GrowthLaccu=0,GrowthSaccu=0,GrowthRaccu=0
        real St,Sw,Ss,SL_rs,SR_rs,Slai
        real Sps,phiN ! Sps is not assigned previous, something is wrong. -JJJJJJJJJJJJJJJJJJJJJ
        real RS, beta_T,Tcold !,Twarm,Topt
        real gamma_W,gamma_Wmax,gamma_T,gamma_Tmax,gamma_N
        real bW,bT,W, RRs0
        ! ----- ouputs   ------- !
        real StemSap, RootSap ! update in this part and maybe use in respiration
        real Rgrowth,Rgroot,Rgleaf,Rgstem
        real L_fall ! the carbon leaving the pools
        real alpha_L,alpha_W,alpha_R
        integer i

        ! Twarm = 35.0
        Tcold = 5.0     !   Tcold=0.0       ! For SPRUCE
        ! Topt  = 30.
        phiN  = 0.33
        bmL   = bmleaf*0.48   ! Carbon
        bmR   = bmRoot*0.48
        bmS   = bmStem*0.48
        if(bmL.lt.NSC/0.333) bmL = NSC/0.333
        if(bmR.lt.NSC/0.333) bmR = NSC/0.333
        if(bmS.lt.NSC/0.334) bmS = NSC/0.334
        StemSap = SapS*bmS  ! Weng 12/05/2008
        RootSap = SapR*bmR
        if(StemSap.lt.0.001)StemSap = 0.001
        if(RootSap.lt.0.001)RootSap = 0.001
        bmP     = bmL+bmR+bmS					! Plant C biomass 
        acP     = bmL+StemSap+bmS					! Plant available sapwood C  
        CNp0    = bmP/(bmL/CN0(1)+bmR/CN0(3)+bmS/CN0(2))		! Plant CN ratio
        hmax    = 24.19   ! m
        hl0     = 0.00019  ! m2/kg C
        LAIMAX0 = 6.
        la0     = 0.2
        ht      = hmax*(1.-exp(-hl0*bmP))				! Scaling plant C biomass to height
        LAIMAX  = AMAX1(LAIMAX0*(1.-exp(-la0*ht)),LAIMIN+0.1)  ! Scaling plant height to maximum LAI
        !   Phenology
        if((GDD5.gt.gddonset).and.onset.eq.0.and.storage.gt.stor_use) then
            onset = 1
        endif
        if((onset.eq.1).and.(storage.gt.stor_use))then
            if(LAI.lt.LAIMAX) add = stor_use
            storage = storage-add
        else
            add     = 0.0
            onset   = 0
        endif
        if(accumulation.lt.(NSCmax+0.005*RootSap))then
            store   = AMAX1(0.,0.005*NSC)			! 0.5% of nonstructure carbon is stored
        else
            store   = 0.0
        endif
        accumulation = accumulation+store
        ! Scalars for plant growth
        Sps          = Sps*(1.-exp(-phiN*NSN))							! Sps is not assigned previous, something is wrong. -JJJJJJJJJJJJJJJJJJJJJ
        Ss           = AMIN1(1.0,2.*fnsc)
        RRs0          = 1.0                                  ! this seem not the parameter of Rs0, thus change it to RRS0
        RS           = bmR/bmL
        SL_rs        = RS/(RS+RRs0*(2.-W))
        SR_rs        = (RRs0*(2.-W))/(RS+RRs0*(2.-W))
        Slai         = amin1(1.0,2.333*(LAIMAX-LAI)/(LAIMAX-LAIMIN))
        St           = AMAX1(0.0, 1.0-exp(-(Tair-gddonset/10.)/5.0))  !0.5 !
        ! Sw=AMAX1(0.333, 0.333+omega)
        Sw           = AMIN1(0.5, AMAX1(0.333, 0.333+omega))
        W            = AMIN1(1.0,3.333*omega)
        ! Plant growth and allocation, based on LM3V
        GPmax        = (GLmax*bmL+GSmax*StemSap+GRmax*bmR) !/acP					
        GrowthP      = AMIN1(GPmax*fnsc*St*(1.-exp(-NSN)), 0.004*NSC, 0.004*NSN*CNp0)
        GrowthL      = MAX(0.0,GrowthP*0.5)      ! updated when QC leaf and wood changed due to the change of plot area for tree biomass
        GrowthR      = MIN(GrowthP*0.4,MAX(0.0,0.75/Sw*bmL-bmR))  ! *c1/(1.+c1+c2)
        GrowthS      = MAX(0.0,GrowthP - (GrowthL+GrowthR) )         ! *c2/(1.+c1+c2)

        NPP         = GrowthL + GrowthR + GrowthS + add       ! Modified by Jiang Jiang 2015/10/13
        addaccu     = addaccu+add
        GrowthLaccu = GrowthLaccu+GrowthL
        GrowthRaccu = GrowthRaccu+GrowthR
        GrowthSaccu = GrowthSaccu+GrowthS
        if(NPP.eq.0.0)then
            alpha_L = 0.333
            alpha_W = 0.333
            alpha_R = 0.333
        else
            alpha_L = (GrowthL+add)/NPP     
            alpha_W = GrowthS/NPP
            alpha_R = GrowthR/NPP
        endif
        ! Carbon cost for growth
        ! Rgrowth,Rgroot,Rgleaf,Rgstem, 0.5 is from IBIS and Amthor, 1984
        Rgleaf  = 0.5*GrowthL
        Rgstem  = 0.5*GrowthS
        Rgroot  = 0.5*GrowthR
        Rgrowth = Rgleaf+Rgstem+Rgroot
        ! Leaf litter 
        gamma_Wmax = 0.12/24. ! maxmum leaf fall rate per hour
        gamma_Tmax = 0.12/24.
        bW         = 4.0
        bT         = 2.0
        if(Tair.gt.(Tcold+10.)) then
            beta_T = 1.
        else 
            if(Tair.gt.Tcold) beta_T = (Tair-Tcold)/10.
            if(Tair.LE.Tcold) beta_T = 0.0
        endif

        if (tauC(1) < 8760.)then
                gamma_W = (1. - W)     **bW * gamma_Wmax
                gamma_T = (1. - beta_T)**bT * gamma_Tmax
        else
                gamma_W = 0.
                gamma_T = 0.
        endif
        gamma_N = 1.0/tauC(1)*Sw      ! Modify by Jiang Jiang 2015/10/20
        if(LAI < LAIMIN) then
            gamma_W = 0.
            gamma_T = 0.
            gamma_N = 0.
        endif
        ! L_fall=bmleaf*0.48*AMIN1((gamma_T+gamma_N),0.99)
        L_fall = bmleaf*0.48*gamma_N
        ! write(*,*)"test-Lfall:", L_fall, bmleaf, gamma_N, tauC(1), Sw
        ! write(*,*)"test-sw:", omega
        ! write(*,*)"test-RootSap:",SapR,bmR, bmRoot, NSC
    return
    end subroutine veg_growth
    
    ! ==========================================================
    subroutine  yrday(doy, hour, fbeam)
        implicit none
        integer doy, hour
        real pidiv, slatx, sindec, cosdec
        real a, b, sinbet, solext, tmprat, tmpR, tmpK
        real fdiff
        real fbeam  ! output
        ! ---------------------------------------------
        pidiv  = pi/180.0
        slatx  = lat*pidiv
        sindec = -sin(23.4*pidiv)*cos(2.0*pi*(doy+10.0)/365.0)
        cosdec = sqrt(1.-sindec*sindec)
        a      = sin(slatx)*sindec
        b      = cos(slatx)*cosdec
        sinbet = a+b*cos(2*pi*(hour-12.)/24.)
        solext = 1370.0*(1.0+0.033*cos(2.0*pi*(doy-10.)/365.0))*sinbet 
        tmprat = radsol/solext
        tmpR   = 0.847-1.61*sinbet+1.04*sinbet*sinbet
        tmpK   = (1.47-tmpR)/1.66
        if(tmprat.le.0.22) fdiff=1.0
        if(tmprat.gt.0.22.and.tmprat.le.0.35) then
            fdiff = 1.0-6.4*(tmprat-0.22)*(tmprat-0.22)
        endif
        if(tmprat.gt.0.35.and.tmprat.le.tmpK) then
            fdiff = 1.47-1.66*tmprat
        endif
        if(tmprat.ge.tmpK) then
            fdiff = tmpR
        endif
        fbeam     = 1.0-fdiff
        if(fbeam.lt.0.0) fbeam=0.0
    return
    end subroutine yrday

    ! -----------------------------------------------------------
    subroutine xlayers(coszen, FLAIT, fbeam, Radabv, &
                    &   Rsoilab1, Rsoilab2, QLleaf, QLair, raero, & ! outputs for soil_thermal
                    &   Acan1, Acan2, Ecan1, Ecan2, Esoil, Hsoil)   ! outputs
        implicit none
        real coszen, FLAIT, fbeam, Radabv(2) ! from canopy
        !-----local vars -------
        real    Gaussx(5),Gaussw(5), vcmxx, eJmxx
        real    xphi1, xphi2, funG, extKb
        real    pi180, cozen15, cozen45, cozen75, xK15, xK45, xK75
        real    transd, extkd, extkn
        integer nw, ng
        real    rhoc(3,2),reff(3,2),kpr(3,2),scatt(2)       !Goudriaan
        real    rhoch, rhoc15, rhoc45, rhoc75
        real    flai
        real    emair, Rnstar(2), grdn ! from Radiso
        real    Qd0, Qb0, Qabs(3,2) ! origial goundriaan module
        real    windux, scalex
        ! ---- NOT local vars ----
        real raero
        ! -----------------------------------------------------------------------
        real fslt, fshd
        !------------------------------------------------------------------------
        ! maybe summary?
        real Tleaf(2), Tleaf1, Tleaf2, Tlk1, Tlk2 ! Tleaf from agsean_day or agsean_ngt
        real FLAIT1, Rsoilabs, Rsoilab1, Rsoilab2, Rsoilab3
        real QLair, QLleaf, QLsoil
        ! outputs
        real Aleaf(2), Eleaf(2) ! from agsean 
        real Acan1, Acan2, Ecan1, Ecan2, Esoil, Hsoil ! Esoil is inout variable
        ! ! from agsean_day to photosyn
        real rhocp, H2OLv, slope, psyc ,Cmolar
        ! real CO2Cs, Qapar, Tlk, weighJ, gsc0, weighR, Gbc, Dleaf
        real fw1, WILTPT, FILDCP, Rsoil, rLAI
        ! --------------------------------------------------------------------------
        ! Normalised Gaussian points and weights (Goudriaan & van Laar, 1993, P98)
        !* 5-point
        data Gaussx/0.0469101,0.2307534,0.5,0.7692465,0.9530899/
        data Gaussw/0.1184635,0.2393144,0.2844444,0.2393144,0.1184635/
        ! soil water conditions
        WILTPT = wsmin/100.
        FILDCP = wsmax/100.
        Tleaf1 = 0.
        Tleaf2 = 0.     
        Acan1  = 0.0        !CO2
        Acan2  = 0.0
        Ecan1  = 0.0        !Evap
        Ecan2  = 0.0                                  
        raero  = 50./windU  ! aerodynamic resistance
        xphi1  = 0.5 - 0.633*xfang -0.33*xfang*xfang            ! Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
        xphi2  = 0.877 * (1.0 - 2.0*xphi1)
        funG   = xphi1 + xphi2*coszen                           ! G-function: Projection of unit leaf area in direction of beam
         if(coszen.gt.0) then                                   ! check if day or night
            extKb = funG/coszen                                 ! beam extinction coeff - black leaves
        else
            extKb = 100.
        end if
        ! Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
        ! Effective extinction coefficient for diffuse radiation Goudriaan & van Laar Eq 6.6)
        pi180   = 3.1416/180.
        cozen15 = cos(pi180*15)
        cozen45 = cos(pi180*45)
        cozen75 = cos(pi180*75)
        xK15    = xphi1/cozen15+xphi2
        xK45    = xphi1/cozen45+xphi2
        xK75    = xphi1/cozen75+xphi2
        transd  = 0.308*exp(-xK15*FLAIT)+0.514*exp(-xK45*FLAIT)+     &
                  &       0.178*exp(-xK75*FLAIT)
        extkd   = (-1./FLAIT)*alog(transd)
        extkn   = extkd                        !N distribution coeff 
        !canopy reflection coefficients (Array indices: first;  1=VIS,  2=NIR
        !                                               second; 1=beam, 2=diffuse
        do nw = 1,2                                                         ! nw:1=VIS, 2=NIR
            scatt(nw)  = trnsmtL(nw) + rhoL(nw)                             ! scattering coeff
            if((1.-scatt(nw))<0.0) scatt(nw) = 0.9999                       ! Weng 10/31/2008
            kpr(nw,1)  = extKb*sqrt(1.-scatt(nw))                           ! modified k beam scattered (6.20)
            kpr(nw,2)  = extkd*sqrt(1.-scatt(nw))                           ! modified k diffuse (6.20)
            rhoch      = (1.-sqrt(1.-scatt(nw)))/(1.+sqrt(1.-scatt(nw)))    ! canopy reflection black horizontal leaves (6.19)
            rhoc15     = 2.*xK15*rhoch/(xK15+extkd)                         ! canopy reflection (6.21) diffuse
            rhoc45     = 2.*xK45*rhoch/(xK45+extkd)
            rhoc75     = 2.*xK75*rhoch/(xK75+extkd)
            rhoc(nw,2) = 0.308*rhoc15+0.514*rhoc45+0.178*rhoc75
            rhoc(nw,1) = 2.*extKb/(extKb+extkd)*rhoch                       ! canopy reflection (6.21) beam 
            reff(nw,1) = rhoc(nw,1)+(rhoS(nw)-rhoc(nw,1))   &               ! effective canopy-soil reflection coeff - beam (6.27)
                         &  *exp(-2.*kpr(nw,1)*FLAIT) 
            reff(nw,2) = rhoc(nw,2)+(rhoS(nw)-rhoc(nw,2))   &               ! effective canopy-soil reflection coeff - diffuse (6.27)
                         &  *exp(-2.*kpr(nw,2)*FLAIT)    
        enddo
        ! isothermal net radiation & radiation conductance at canopy top - needed to calc emair
        call Radiso(flai, FLAIT, extkd, fbeam, Qabs, emair, Rnstar, grdn) ! emair, Rnstar, grdn  are outputs
        do ng = 1, 5
            flai = gaussx(ng)*FLAIT
            ! -------------------------------------------------------------------------------
            ! zhou: put goudriaan subroutine here ....No need the subroutine of goudriaan
            ! radiation absorption for visible and near infra-red
            do nw = 1,2
                Qd0        = (1.-fbeam)*radabv(nw)                                           !diffuse incident radiation
                Qb0        = fbeam*radabv(nw)                                                !beam incident radiation
                Qabs(nw,2) = Qd0*(kpr(nw,2)*(1.-reff(nw,2))*exp(-kpr(nw,2)*FLAI))+         & !absorbed radiation - shaded leaves, diffuse
                                &   Qb0*(kpr(nw,1)*(1.-reff(nw,1))*exp(-kpr(nw,1)*FLAI)-   & !beam scattered
                                &   extKb*(1.-scatt(nw))*exp(-extKb*FLAI))
                Qabs(nw,1) = Qabs(nw,2)+extKb*Qb0*(1.-scatt(nw))                             !absorbed radiation - sunlit leaves 
            end do
            !----------------------------------------------------------------------------------
            ! isothermal net radiation & radiation conductance at canopy top
            call Radiso(flai, FLAIT, extkd, fbeam, Qabs, emair, Rnstar, grdn)
            windUx = windU*exp(-extkU*flai)             !windspeed at depth xi
            scalex = exp(-extkn*flai)                   !scale Vcmx0 & Jmax0
            Vcmxx  = Vcmx0*scalex
            eJmxx  = eJmx0*scalex
            if(radabv(1).ge.10.0) then                          !check solar Radiation > 10 W/m2
                ! leaf stomata-photosynthesis-transpiration model - daytime
                call agsean_day(windUx, Qabs, grdn, raero, Rnstar,Vcmxx, eJmxx, &
                          & Aleaf, Tleaf, Eleaf)       
            else
                call agsean_ngt(windUx, Qabs, grdn, raero, Rnstar, vcmxx, eJmxx, &
                          & Aleaf, Tleaf, Eleaf)
            endif
            fslt      = exp(-extKb*flai)                        !fraction of sunlit leaves
            fshd      = 1.0-fslt                                !fraction of shaded leaves
            ! Rnst1     = Rnst1+fslt*Rnstar(1)*Gaussw(ng)*FLAIT  !Isothermal net rad`
            ! Rnst2     = Rnst2+fshd*Rnstar(2)*Gaussw(ng)*FLAIT
            ! RnstL(ng) = Rnst1+Rnst2
            !
            ! Qcan1 = Qcan1+fslt*Qabs(1,1)*Gaussw(ng)*FLAIT  !visible
            ! Qcan2 = Qcan2+fshd*Qabs(1,2)*Gaussw(ng)*FLAIT
            ! QcanL(ng)=Qcan1+Qcan2
            ! !
            ! Rcan1=Rcan1+fslt*Qabs(2,1)*Gaussw(ng)*FLAIT  !NIR
            ! Rcan2=Rcan2+fshd*Qabs(2,2)*Gaussw(ng)*FLAIT
            ! RcanL(ng)=Rcan1+Rcan2
            !
            ! write(*,*)"test-Aleaf:", Aleaf
            if(Aleaf(1).lt.0.0) Aleaf(1)=0.0      !Weng 2/16/2006
            if(Aleaf(2).lt.0.0) Aleaf(2)=0.0      !Weng 2/16/2006

            Acan1 = Acan1+fslt*Aleaf(1)*Gaussw(ng)*FLAIT*stomN    !amphi/hypostomatous
            Acan2 = Acan2+fshd*Aleaf(2)*Gaussw(ng)*FLAIT*stomN
            ! AcanL(ng)=Acan1+Acan2

            ! layer1(ng)=Aleaf(1)
            ! layer2(ng)=Aleaf(2)

            Ecan1     = Ecan1+fslt*Eleaf(1)*Gaussw(ng)*FLAIT
            Ecan2     = Ecan2+fshd*Eleaf(2)*Gaussw(ng)*FLAIT
            ! EcanL(ng) = Ecan1+Ecan2
            ! !
            ! Hcan1=Hcan1+fslt*Hleaf(1)*Gaussw(ng)*FLAIT
            ! Hcan2=Hcan2+fshd*Hleaf(2)*Gaussw(ng)*FLAIT
            ! HcanL(ng)=Hcan1+Hcan2
            ! !
            ! Gbwc1=Gbwc1+fslt*gbleaf(1)*Gaussw(ng)*FLAIT*stom_n
            ! Gbwc2=Gbwc2+fshd*gbleaf(2)*Gaussw(ng)*FLAIT*stom_n
            ! !
            ! Gswc1=Gswc1+fslt*gsleaf(1)*Gaussw(ng)*FLAIT*stom_n
            ! Gswc2=Gswc2+fshd*gsleaf(2)*Gaussw(ng)*FLAIT*stom_n
            ! !
            Tleaf1=Tleaf1+fslt*Tleaf(1)*Gaussw(ng)*FLAIT
            Tleaf2=Tleaf2+fshd*Tleaf(2)*Gaussw(ng)*FLAIT    ! Tleaf IS INF
            ! write(*,*)"test-Tleaf:", Tleaf1, Tleaf2, fslt,fshd, Tleaf(1), Tleaf(2), Gaussw(ng) 
        enddo
        ! write(*,*)"test - Qabs:", Qabs, FLAIT
        ! write(*,*)"test-Acan", Acan1, Acan2, fslt, Aleaf, Flait, stomN
        FLAIT1 = (1.0-exp(-extKb*FLAIT))/extkb
        Tleaf1 = Tleaf1/FLAIT1
        Tleaf2 = Tleaf2/(FLAIT-FLAIT1)
        ! write(*,*)"test-Tleaf", FLAIT1, FLAIT ! NO error in FLAIT1 AND FLAIT
        ! Soil surface energy and water fluxes
        ! Radiation absorbed by soil
        Rsoilab1 = fbeam*(1.-reff(1,1))*exp(-kpr(1,1)*FLAIT)        &
                    &  +(1.-fbeam)*(1.-reff(1,2))*exp(-kpr(1,2)*FLAIT)          !visible
        Rsoilab2 = fbeam*(1.-reff(2,1))*exp(-kpr(2,1)*FLAIT)        &
                    &  +(1.-fbeam)*(1.-reff(2,2))*exp(-kpr(2,2)*FLAIT)          !NIR
        Rsoilab1 = Rsoilab1*Radabv(1)
        Rsoilab2 = Rsoilab2*Radabv(2)
      
        Tlk1     = Tleaf1+273.2
        Tlk2     = Tleaf2+273.2
        ! write(*,*)"test-Tlk", Tleaf1, Tleaf1
        ! temp1=-extkd*FLAIT
        QLair    = emair*sigma*(TairK**4)*exp(-extkd*FLAIT)
        QLleaf   = emleaf*sigma*(Tlk1**4)*exp(-extkb*FLAIT)           &
                    &   +emleaf*sigma*(Tlk2**4)*(1.0-exp(-extkb*FLAIT))
        ! write(*,*)"test-QLleaf1", QLleaf, emleaf, sigma, Tlk1, extkb, FLAIT, Tlk2 ! Tlk1 and Tlk2 IS NAN
        QLleaf   = QLleaf*(1.0-exp(-extkd*FLAIT)) 
        ! write(*,*)"test-QLleaf1", extkd, FLAIT
        QLsoil   = emsoil*sigma*(TairK**4)
        Rsoilab3 = (QLair+QLleaf)*(1.0-rhoS(3))-QLsoil
        ! write(*,*)"test-Rsoilabs", QLair, QLleaf, rhoS(3), QLsoil  ! Qleaf IS INF
        ! Total radiation absorbed by soil    
        Rsoilabs = Rsoilab1+Rsoilab2+Rsoilab3  ! Rsoilab3 is inf
        ! write(*,*)"test-Rsoilabs", Rsoilab1, Rsoilab2, Rsoilab3

        rhocp    = cpair*Patm*AirMa/(Rconst*TairK)
        H2OLv    = H2oLv0-2.365e3*Tair
        slope    = (esat(Tair+0.1)-esat(Tair))/0.1
        psyc     = Patm*cpair*AirMa/(H2OLv*H2OMw)
        Cmolar   = Patm/(Rconst*TairK)

        fw1      = AMIN1(AMAX1((FILDCP-wcl(1))/(FILDCP-WILTPT),0.05),1.0)
        Rsoil    = 30.*exp(0.2/fw1)
        rLAI     = exp(FLAIT)
        ! latent heat flux into air from soil
        !       Eleaf(ileaf)=1.0*
        ! &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    !2* Weng 0215
        ! &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
        Esoil    = (slope*(Rsoilabs-G)+rhocp*Dair/(raero+rLAI))/       &
                    &       (slope+psyc*(rsoil/(raero+rLAI)+1.))
        ! sensible heat flux into air from soil
        Hsoil    = Rsoilabs-Esoil-G
        ! ESoil and Rsoilabs IS inf
        ! write(*,*)"test - Esoil:", Esoil, slope, Rsoilabs, G, rhocp, Dair, raero, rLAI, rsoil, psyc
    return
    end subroutine xlayers

    ! -----------------------------------------------------------
    subroutine Radiso(flai, FLAIT, extkd, fbeam, Qabs, emair, Rnstar, grdn)
        implicit none
        real flai, FLAIT, extkd, fbeam ! zhou: no 
        ! ---- local vars ------
        real emcloud
        real rhocp, emsky, ep8z, tau8
        real Bn0, Bnxi
        ! ---- output? ----
        real emair
        real Rnstar(2)
        real Qabs(3,2)
        real grdn
        rhocp   = cpair*Patm*airMa/(Rconst*TairK)   !volumetric heat capacity (J/m3/K)
        emsky   = 0.642*(eairP/Tairk)**(1./7)       !! apparent atmospheric emissivity for clear skies (Brutsaert, 1975) note eair in Pa
        ep8z    = 0.24+2.98e-12*eairP*eairP*exp(3000/TairK)     ! apparent emissivity from clouds (Kimball et al 1982)
        tau8    = amin1(1.0,1.0-ep8z*(1.4-0.4*ep8z))            !ensure tau8<1
        emcloud = 0.36*tau8*(1.-fbeam)*(1-10./TairK)**4      !10 from Tcloud = Tair-10
        ! apparent emissivity from sky plus clouds      
        !      emair=emsky+emcloud
        ! 20/06/96
        emair   = emsky
        if(emair.gt.1.0) emair = 1.0
        ! net isothermal outgoing longwave radiation per unit leaf area at canopy
        ! top & thin layer at flai (Note Rn* = Sn + Bn is used rather than Rn* = Sn - Bn in Leuning et al 1985)
        Bn0  = sigma*(TairK**4.)
        Bnxi = Bn0*extkd*(exp(-extkd*flai)*(emair-emleaf)       &
                &    + exp(-extkd*(flait-flai))*(emsoil-emleaf))
        ! isothermal net radiation per unit leaf area for thin layer of sunlit and
        ! shaded leaves
        Rnstar(1) = Qabs(1,1)+Qabs(2,1)+Bnxi
        Rnstar(2) = Qabs(1,2)+Qabs(2,2)+Bnxi
        ! radiation conductance (m/s) @ flai
        grdn = 4.*sigma*(TairK**3.)*extkd*emleaf*               &       ! corrected by Jiang Jiang 2015/9/29
                &    (exp(-extkd*flai)+exp(-extkd*(flait-flai)))       &
                &    /rhocp
    return
    end subroutine Radiso

    ! ----------------------------------------------------------------------------

    ! subroutine goudriaan(FLAI,coszen,radabv,fbeam,reff,kpr,   &
    !                         & scatt,Qabs)
    !     ! for spheric leaf angle distribution only
    !     ! compute within canopy radiation (PAR and near infra-red bands)
    !     ! using two-stream approximation (Goudriaan & vanLaar 1994)
    !     ! tauL: leaf transmittance
    !     ! rhoL: leaf reflectance
    !     ! rhoS: soil reflectance
    !     ! sfang XiL function of Ross (1975) - allows for departure from spherical LAD
    !     !     (-1 vertical, +1 horizontal leaves, 0 spherical)
    !     ! FLAI: canopy leaf area index
    !     ! funG: Ross' G function
    !     ! scatB: upscatter parameter for direct beam
    !     ! scatD: upscatter parameter for diffuse
    !     ! albedo: single scattering albedo
    !     ! output:
    !     ! Qabs(nwave,type), nwave=1 for visible; =2 for NIR,
    !     !                     type=1 for sunlit;   =2 for shaded (W/m2)
    !     real FLAI, coszen, fbeam
    !     real radabv(2)
    !     real Qabs(3,2),reff(3,2),kpr(3,2),scatt(2)
    !     ! -------------------------------------------------
    !     real xu, xphi1, xphi2, funG, extKb
    !     integer nw
    !     real Qb0, Qd0
    !     xu=coszen                                         !cos zenith angle

    !     !     Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
    !     xphi1 = 0.5 - 0.633*xfang -0.33*xfang*xfang
    !     xphi2 = 0.877 * (1.0 - 2.0*xphi1)
    !     funG  = xphi1 + xphi2*xu                             !G-function: Projection of unit leaf area in direction of beam

    !     if(coszen.gt.0) then                                  !check if day or night
    !         extKb=funG/coszen                                   !beam extinction coeff - black leaves
    !     else
    !         extKb=100.
    !     end if
                        
    !     ! Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
    !     do nw=1,2
    !         Qd0=(1.-fbeam)*radabv(nw)                                          !diffuse incident radiation
    !         Qb0=fbeam*radabv(nw)                                               !beam incident radiation
    !         Qabs(nw,2)=Qd0*(kpr(nw,2)*(1.-reff(nw,2))*exp(-kpr(nw,2)*FLAI))+  & !absorbed radiation - shaded leaves, diffuse
    !         &            Qb0*(kpr(nw,1)*(1.-reff(nw,1))*exp(-kpr(nw,1)*FLAI)-   & !beam scattered
    !         &            extKb*(1.-scatt(nw))*exp(-extKb*FLAI))
    !         Qabs(nw,1)=Qabs(nw,2)+extKb*Qb0*(1.-scatt(nw))                     !absorbed radiation - sunlit leaves 
    !     end do
    !     return
    ! end subroutine goudriaan

    ! -----------------------------------------------------------
    subroutine agsean_day(windUx, Qabs, grdn, raero, Rnstar, Vcmxx, eJmxx,  &
                          & Aleaf, Tleaf, Eleaf) 
        implicit none
        real    windUx, Qabs(3,2), grdn, raero, Rnstar(2) ! from xlayers 
        real    Vcmxx, eJmxx ! from xlayer to photosyn
        real    Aleafx, gscx ! from photosyn
        ! ---- local vars -----
        real    rhocp, H2OLv, slope, psyc
        real    Cmolar, weighJ, gbHu
        integer ileaf, kr1
        real    Tlk, Dleaf, co2cs, Qapar
        real    Gras, gbHf, gbH, rbH, rbw, rbH_L, rrdn
        real    Y, gbc, gsc0, varQc, weighR
        real    Aleaf(2), co2Ci(2), Eleaf(2), Hleaf(2) !
        real    gsw, gswv, rswv, Tlk1
        real    gbleaf(2), gsleaf(2)
        ! ---- use and output
        real    Tleaf(2)
        !-------------------------------------------
        rhocp  = cpair*Patm*AirMa/(Rconst*TairK)
        H2OLv  = H2oLv0-2.365e3*Tair
        slope  = (esat(Tair+0.1)-esat(Tair))/0.1
        psyc   = Patm*cpair*AirMa/(H2OLv*H2OMw)
        Cmolar = Patm/(Rconst*TairK)
        weighJ = 1.0
        ! boundary layer conductance for heat - single sided, forced convection
        ! (Monteith 1973, P106 & notes dated 23/12/94)
        if(windUx/wleaf >= 0.0)then
            gbHu = 0.003*sqrt(windUx/wleaf)    !m/s
        else
            gbHu = 0.003 !*sqrt(-windUx/wleaf)
        endif         ! Weng 10/31/2008
        ! write(*,*)"test_parameters:", Tair, Dair, co2ca, Qabs, Dheat, gsw0
        ! ---------------------------------------------------------------------
        do ileaf = 1,2                             ! loop over sunlit and shaded leaves
            Tleaf(ileaf) = Tair                    ! first estimate of leaf temperature - assume air temp
            Tlk          = Tleaf(ileaf)+273.2      ! Tleaf to deg K
            Dleaf        = Dair                    ! first estimate of deficit at leaf surface - assume Da   !Pa           
            co2cs        = co2ca                   ! first estimate for co2cs !mol/mol
            Qapar        = (4.6e-6)*Qabs(1,ileaf)
            ! -----------------------------------
            kr1          = 0                     !iteration counter for LE
            ! return point for evaporation iteration
            do               !iteration for leaf temperature
                ! single-sided boundary layer conductance - free convection (see notes 23/12/94)
                Gras   = 1.595e8*ABS(Tleaf(ileaf)-Tair)*(wleaf**3.)     !Grashof
                gbHf   = 0.5*Dheat*(Gras**0.25)/wleaf
                gbH    = gbHu+gbHf                         !m/s
                rbH    = 1./gbH                            !b/l resistance to heat transfer
                rbw    = 0.93*rbH                          !b/l resistance to water vapour
                ! Y factor for leaf: stom_n = 1.0 for hypostomatous leaf;  stom_n = 2.0 for amphistomatous leaf
                rbH_L  = rbH*stomN/2.                   !final b/l resistance for heat  
                rrdn   = 1./grdn
                Y      = 1./(1.+ (rbH_L+raero)/rrdn)
                ! boundary layer conductance for CO2 - single side only (mol/m2/s)
                gbc    = Cmolar*gbH/1.32            !mol/m2/s
                gsc0   = gsw0/1.57                 !convert conductance for H2O to that for CO2
                varQc  = 0.0
                weighR = 1.0
                call photosyn(co2cs, Qapar, Tlk, weighJ, gsc0, weighR, Gbc, Dleaf, Vcmxx, eJmxx, &
                        & Aleafx, gscx)  !outputs
                ! choose smaller of Ac, Aq
                Aleaf(ileaf) = Aleafx      !0.7 Weng 3/22/2006          !mol CO2/m2/s
                ! calculate new values for gsc, cs (Lohammer model)
                co2cs        = co2ca-Aleaf(ileaf)/gbc
                co2Ci(ileaf) = co2cs-Aleaf(ileaf)/gscx
                ! scale variables
                gsw  = gscx*1.56       !gsw in mol/m2/s, oreginal:gsw=gscx*1.56,Weng20090226
                gswv = gsw/Cmolar                           !gsw in m/s
                rswv = 1./gswv
                ! calculate evap'n using combination equation with current estimate of gsw
                Eleaf(ileaf) = 1.0*                                         &
                    &  (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/  &   !2* Weng 0215
                    &  (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))        
                Hleaf(ileaf) = Y*(Rnstar(ileaf)-Eleaf(ileaf))     ! calculate sensible heat flux
                Tlk1         = 273.2+Tair+Hleaf(ileaf)*(rbH/2.+raero)/rhocp ! calculate new leaf temperature (K)
                ! write(*,*)"test:", Hleaf, rbH
                ! write(*,*)"test:", Gras, Tleaf(ileaf), Tair, wleaf 
                ! calculate Dleaf use LE=(rhocp/psyc)*gsw*Ds
                Dleaf         = psyc*Eleaf(ileaf)/(rhocp*gswv)
                gbleaf(ileaf) = gbc*1.32*1.075
                gsleaf(ileaf) = gsw
                ! compare current and previous leaf temperatures
                ! write(*,*)"test-Tkl:", Tlk1, Tlk, Tlk1-Tlk
                if(abs(Tlk1-Tlk).le.0.1) exit ! original is 0.05 C Weng 10/31/2008
                ! update leaf temperature  ! leaf temperature calculation has many problems! Weng 10/31/2008
                Tlk          = Tlk1
                Tleaf(ileaf) = Tlk1-273.2
                kr1          = kr1+1
                if(kr1 > 500)then
                    Tlk = TairK
                    exit
                endif
                if(Tlk < 200.)then
                    Tlk=TairK
                    exit 
                endif                     ! Weng 10/31/2008
                ! goto 100                          !solution not found yet
            enddo
            ! write(*,*)"test-day-Tleaf:", Tleaf, Tlk1
    ! 10  continue
        enddo
    end subroutine agsean_day

    !-----------------------------------------------------------
    subroutine agsean_ngt(windUx, Qabs, grdn, raero, Rnstar, vcmxx, eJmxx, &
                          & Aleaf, Tleaf, Eleaf)
        implicit none
        real windUx, Qabs(3,2), grdn, raero, Rnstar(2), vcmxx, eJmxx  ! from xlayers 
        real Aleafx, gscx                               ! from photosyn
        ! ---- local vars -----
        real    rhocp, H2OLv, slope, psyc, Cmolar, weighJ, gbHu
        integer ileaf, kr1
        real    Tlk, Dleaf, co2cs, Qapar
        real    Gras, gbHf, gbH, rbH, rbw, rbH_L, rrdn
        real    Y, gbc, gsc0, varQc, weighR
        real    co2Ci(2), Eleaf(2), Hleaf(2), Aleaf(2)
        real    gsc, gsw, gswv, rswv, Tlk1
        real    gbleaf(2), gsleaf(2)
        ! ---- use and output
        real Tleaf(2)
        ! ----------------------------------------------------
        rhocp  = cpair*Patm*AirMa/(Rconst*TairK)
        H2OLv  = H2oLv0-2.365e3*Tair
        slope  = (esat(Tair+0.1)-esat(Tair))/0.1
        psyc   = Patm*cpair*AirMa/(H2OLv*H2OMw)
        Cmolar = Patm/(Rconst*TairK)
        weighJ = 1.0
        ! boundary layer conductance for heat - single sided, forced convection
        ! (Monteith 1973, P106 & notes dated 23/12/94)
        gbHu   = 0.003*sqrt(windUx/wleaf)    !m/s
        ! raero=0.0                        !aerodynamic resistance s/m
        do ileaf = 1,2                  ! loop over sunlit and shaded leaves
            Tleaf(ileaf) = Tair                 ! first estimate of leaf temperature - assume air temp
            Tlk          = Tleaf(ileaf)+273.2    ! Tleaf to deg K
            Dleaf        = Dair                  ! first estimate of deficit at leaf surface - assume Da !Pa            
            co2cs        = co2ca               ! first estimate for co2cs !mol/mol
            Qapar        = (4.6e-6)*Qabs(1,ileaf)
            ! -------------------------------------------------
            kr1          = 0                     !iteration counter for LE
            do
                ! 100        continue !    return point for evaporation iteration
                ! single-sided boundary layer conductance - free convection (see notes 23/12/94)
                Gras   = 1.595e8*abs(Tleaf(ileaf)-Tair)*(wleaf**3)     !Grashof
                gbHf   = 0.5*Dheat*(Gras**0.25)/wleaf
                gbH    = gbHu + gbHf                         !m/s
                rbH    = 1./gbH                            !b/l resistance to heat transfer
                rbw    = 0.93*rbH                          !b/l resistance to water vapour
                ! Y factor for leaf: stom_n = 1.0 for hypostomatous leaf;  stom_n = 2.0 for amphistomatous leaf
                rbH_L  = rbH*stomN/2.                   !final b/l resistance for heat  
                rrdn   = 1./grdn
                Y      = 1./(1.+ (rbH_L+raero)/rrdn)
                ! boundary layer conductance for CO2 - single side only (mol/m2/s)
                gbc    = Cmolar*gbH/1.32            !mol/m2/s
                gsc0   = gsw0/1.57                        !convert conductance for H2O to that for CO2
                varQc  = 0.0                  
                weighR = 1.0
                ! respiration      
                Aleafx = -0.0089*Vcmxx*exp(0.069*(Tlk-293.2))
                gsc    = gsc0
                ! choose smaller of Ac, Aq
                Aleaf(ileaf) = Aleafx                     !mol CO2/m2/s
                ! calculate new values for gsc, cs (Lohammer model)
                co2cs = co2ca-Aleaf(ileaf)/gbc
                co2Ci(ileaf) = co2cs-Aleaf(ileaf)/gsc
                ! scale variables
                gsw    = gsc*1.56                              !gsw in mol/m2/s
                gswv   = gsw/Cmolar                           !gsw in m/s
                rswv   = 1./gswv
                ! calculate evap'n using combination equation with current estimate of gsw
                Eleaf(ileaf) = (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/   &
                                &      (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
                ! calculate sensible heat flux
                Hleaf(ileaf) = Y*(Rnstar(ileaf)-Eleaf(ileaf))
                ! calculate new leaf temperature (K)
                Tlk1         = 273.2+Tair+Hleaf(ileaf)*(rbH/2.+raero)/rhocp
                ! write(*,*)"test-TLK1:",Tlk1, Tair, Hleaf(ileaf), rbH, raero, rhocp
                ! calculate Dleaf use LE=(rhocp/psyc)*gsw*Ds
                Dleaf         = psyc*Eleaf(ileaf)/(rhocp*gswv)
                gbleaf(ileaf) = gbc*1.32*1.075
                gsleaf(ileaf) = gsw
                ! compare current and previous leaf temperatures
                if(abs(Tlk1-Tlk).le.0.1) exit
                if(kr1.gt.500) exit
                ! update leaf temperature
                Tlk          = Tlk1 
                Tleaf(ileaf) = Tlk1-273.2
                kr1          = kr1+1
            enddo                          !solution not found yet
    ! 10    continue
    ! write(*,*)"test-night-Tleaf:", Tleaf, Tlk1
        enddo
        
    return
    end subroutine agsean_ngt
    ! ----------------------------------------------------------
    subroutine photosyn(co2cs, Qapar, Tlk, weighJ, gsc0, weighR, Gbc, Dleaf, Vcmxx, eJmxx, &
                        & Aleafx, gscx) ! outputs
        ! calculate Vcmax, Jmax at leaf temp (Eq 9, Harley et al 1992)
        ! turned on by Weng, 2012-03-13
        ! VcmxT = Vjmax(Tlkx,Trefk,Vcmx1,Eavm,Edvm,Rconst,Entrpy)
        ! eJmxT = Vjmax(Tlkx,Trefk,eJmx1,Eajm,Edjm,Rconst,Entrpy)
        implicit none
        real CO2Cs, Qapar, Tlk, weighJ, gsc0, weighR, Gbc, Dleaf ! from agsean_day
        real Vcmxx, eJmxx ! from xlayers
        real Acx,Aqx ! from ciandA
        real Aleafx, gscx ! output
        ! ------- local vars -------
        real TminV, TmaxV, ToptV, TminJ, TmaxJ, ToptJ
        real Tlf, VcmxT, eJmxT, eJ, conKcT, conKoT
        real Rd, Tdiff, gammas, gamma, a1, X, Gma, Bta
        ! real sps ! zhou: not sure which parameter is?
        ! ------------------------------------------------------------------------------
        CO2Cs = AMAX1(CO2Cs,0.6*CO2Ca)
        ! check if it is dark - if so calculate respiration and g0 to assign conductance 
        if(Qapar .le.0.) then                            !night, umol quanta/m2/s
            Aleafx = -0.0089*Vcmxx*exp(0.069*(Tlk-293.2))   ! original: 0.0089 Weng 3/22/2006
            Gscx   = gsc0
        endif
        ! calculate  Vcmax, Jmax at leaf temp using Reed et al (1976) function J appl Ecol 13:925
        TminV = gddonset/10.  ! original -5.        !-Jiang Jiang 2015/10/13
        TmaxV = 50.
        ToptV = 35.
        
        TminJ = TminV
        TmaxJ = TmaxV
        ToptJ = ToptV 
        
        Tlf   = Tlk-273.2
        VcmxT = VJtemp(Tlf,TminV,TmaxV,ToptV,Vcmxx)
        eJmxT = VJtemp(Tlf,TminJ,TmaxJ,ToptJ,eJmxx)      

        eJ     = weighJ*fJQres(eJmxT,alpha,Qapar,theta)    ! calculate J, the asymptote for RuBP regeneration rate at given Q
        conKcT = EnzK(Tlk,Trefk,conKc0,Rconst,Ekc)         ! calculate Kc, Ko, Rd gamma*  & gamma at leaf temp
        conKoT = EnzK(Tlk,Trefk,conKo0,Rconst,Eko)
        ! following de Pury 1994, eq 7, make light respiration a fixed proportion of
        ! Vcmax
        Rd     = 0.0089*VcmxT*weighR                              !de Pury 1994, Eq7
        Tdiff  = Tlk-Trefk
        gammas = gam0*(1.+gam1*Tdiff+gam2*Tdiff*Tdiff)       !gamma*
        ! gamma = (gammas+conKcT*(1.+O2ci/conKoT)*Rd/VcmxT)/(1.-Rd/VcmxT)
        gamma = 0.0
        ! ***********************************************************************
        ! Analytical solution for ci. This is the ci which satisfies supply and demand
        ! functions simultaneously
        ! calculate X using Lohammer model, and scale for soil moisture
        a1 = 1./(1.-0.7)
        X  = a1*fwsoil/((co2cs - gamma)*(1.0 + Dleaf/Ds0))
        ! calculate solution for ci when Rubisco activity limits A
        Gma = VcmxT  
        Bta = conKcT*(1.0+ o2ci/conKoT)
        call ciandA(Gma,Bta,gsc0,X,Rd,co2Cs,gammas,Acx)
        ! calculate +ve root for ci when RuBP regeneration limits A
        Gma = eJ/4.
        Bta = 2.*gammas
        ! calculate coefficients for quadratic equation for ci
        call ciandA(Gma,Bta,gsc0,X,Rd,co2Cs,gammas,Aqx)
        ! choose smaller of Ac, Aq
        sps    = AMAX1(0.001,sps)                  !Weng, 3/30/2006
        Aleafx = (amin1(Acx,Aqx) - Rd) !*sps     ! Weng 4/4/2006
        ! if(Aleafx.lt.0.0) Aleafx=0.0    ! by Weng 3/21/2006
        ! calculate new values for gsc, cs (Lohammer model)
        CO2cs  = co2ca-Aleafx/Gbc
        Gscx   = gsc0 + X*Aleafx  ! revised by Weng
        ! write(*,*)"test-Aqx:", Gma, Bta, eJ, gamma
        ! write(*,*)"test-leafx:", Aleafx, Acx, Aqx, Rd
        return
    end subroutine photosyn
    ! -----------------------------------------------------------
    subroutine ciandA(Gma,Bta,g0,X,Rd,co2Cs,gammas,Aquad)
        implicit none
        real Gma,Bta,g0,X,Rd,co2Cs,gammas,ciquad,Aquad
        ! --------- local vars ------
        real b0, b1, b2, bx
        ! calculate coefficients for quadratic equation for ci
        b2 = g0+X*(Gma-Rd)
        b1 = (1.-co2cs*X)*(Gma-Rd)+g0*(Bta-co2cs)-X*(Gma*gammas+Bta*Rd)
        b0 = -(1.-co2cs*X)*(Gma*gammas+Bta*Rd)-g0*Bta*co2cs
        bx = b1*b1-4.*b2*b0
        if(bx.gt.0.0) then 
            ! calculate larger root of quadratic
            ciquad = (-b1+sqrt(bx))/(2.*b2)
        endif
        IF(ciquad.lt.0.or.bx.lt.0.) THEN
            Aquad = 0.0
            ciquad = 0.7 * co2Cs
        ELSE
            Aquad = Gma*(ciquad-gammas)/(ciquad+Bta)
        ENDIF
        ! write(*,*)"test-ciandA:", Gma,Bta,g0,X,Rd,co2Cs,gammas
    return
    end subroutine ciandA
    ! ----------------------------------------------------------
    ! zhou: calculate some common variables
    subroutine reComVars(rhocp, H2OLv, slope, psyc, Cmolar, weighJ)
        implicit none
        real rhocp, H2OLv, slope, psyc, Cmolar, weighJ 
        rhocp  = cpair*Patm*AirMa/(Rconst*TairK)
        H2OLv  = H2oLv0-2.365e3*Tair
        slope  = (esat(Tair+0.1)-esat(Tair))/0.1
        psyc   = Patm*cpair*AirMa/(H2OLv*H2OMw)
        Cmolar = Patm/(Rconst*TairK)
        weighJ = 1.0
    end subroutine reComVars
    ! -----------------------------------------------------------
    ! functions 
    real function cal_sinbet(doy,timeh)
        integer doy, timeh
        real    rad, sinlat, coslat, sindec, cosdec
        real    A, B
        ! sin(bet), bet = elevation angle of sun
        ! calculations according to Goudriaan & van Laar 1994 P30
        rad = pi/180.
        ! sine and cosine of latitude
        sinlat = sin(rad*lat)
        coslat = cos(rad*lat)
        ! sine of maximum declination
        sindec=-sin(23.45*rad)*cos(2.0*pi*(doy+10.0)/365.0)
        cosdec=sqrt(1.-sindec*sindec)
        ! terms A & B in Eq 3.3
        A = sinlat*sindec
        B = coslat*cosdec
        cal_sinbet = A+B*cos(pi*(timeh-12.)/12.)
        return
    end function cal_sinbet
    ! ! ----------------------------------------------------------
    ! real function esat(T)
    !     real T
    !     ! returns saturation vapour pressure in Pa
    !     esat = 610.78*exp(17.27*T/(T+237.3))
    ! return
    ! end
    ! ----------------------------------------------------------
    ! Reed et al (1976, J appl Ecol 13:925) equation for temperature response
    ! used for Vcmax and Jmax
    real function VJtemp(Tlf,TminVJ,TmaxVJ,ToptVJ,VJmax0)
        real Tlf,TminVJ,TmaxVJ,ToptVJ,VJmax0
        real pwr
        if(Tlf.lt.TminVJ) Tlf=TminVJ   !constrain leaf temperatures between min and max
        if(Tlf.gt.TmaxVJ) Tlf=TmaxVJ
        pwr     = (TmaxVJ-ToptVJ)/(ToptVj-TminVj)
        VJtemp  = VJmax0*((Tlf-TminVJ)/(ToptVJ-TminVJ))*     &
                    &  ((TmaxVJ-Tlf)/(TmaxVJ-ToptVJ))**pwr 
    return
    end
    !------------------------------------------------------------
    real function fJQres(eJmx,alpha,Q,theta)
        real eJmx,alpha,Q,theta
        real AX, BX, CX
        AX = theta                                 !a term in J fn
        BX = alpha*Q+eJmx                          !b term in J fn
        CX = alpha*Q*eJmx                          !c term in J fn
        if((BX*BX-4.*AX*CX)>=0.0)then
            fJQres = (BX-SQRT(BX*BX-4.*AX*CX))/(2*AX)
        else
            fJQres = (BX)/(2*AX)                   !Weng 10/31/2008
        endif
    return
    end
    ! -----------------------------------------------------------
    real function EnzK(Tk,Trefk,EnzK0,Rconst,Eactiv)
        real Tk,Trefk,EnzK0,Rconst,Eactiv
        real temp1
        temp1 = (Eactiv/(Rconst* Trefk))*(1.-Trefk/Tk)
        !      if (temp1<50.)then
        EnzK = EnzK0*EXP((Eactiv/(Rconst* Trefk))*(1.-Trefk/Tk))
        !      else
        !      EnzK = EnzK0*EXP(50.)                                          ! Weng 10/31/2008
        !      endif
    return
    end
end module mod_vegetation