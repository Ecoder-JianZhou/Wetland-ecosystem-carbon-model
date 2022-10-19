module mod_transferCpools
    use mod_dataTypes
    implicit none
    
    contains
    subroutine TCS_CN()
        implicit none
        ! carbon transfer according to Xu et al. 2007 
        ! subroutine TCS_CN(Tair,Tsoil,omega,runoff,&
        ! &               NPP,alpha_L,alpha_W,alpha_R,L_fall,&
        ! &               tauC,QC,OutC,Rh_pools,Rnitrogen,NSC,&
        ! &               CNmin,CNmax,NSNmax,NSNmin,alphaN,    &        ! nitrogen
        ! &               NSN,N_uptake,N_miner,QN,QNminer,&
        ! &               CN,CN0,fnsc,rdepth,N_deficit,&
        ! &               N_leaf,N_wood,N_root,N_LF,N_WF,N_RF,&
        ! &               N_deposit,N_fixation,N_leach,N_vol,&
        ! &               SNvcmax,SNgrowth,SNRauto,SNrs,Q10,&
        ! &               tsoill,testout,doSoilphy) 
        
        ! local vars
        real CNmin, CNmax   !, NSNmax, NSNmin
        real N_miner, N_leaf,N_wood,N_root
        real N_LF,N_WF,N_RF ! maybe use to summary
        real N_fixation, N_leach, N_vol,SNgrowth
        real Q10h(5)
        real NPP_L,NPP_W,NPP_R
        real Q_plant
        real etaL,etaW,etaR                ! the percentage of fine litter of the litters from plant parts; zhou: seem no use?
        real f_F2M,f_C2M,f_C2S,f_M2S,f_M2P,f_S2P,f_S2M,f_P2M ! zhou: put them to site-based parameters?
        real CN_foliage, N_demand,N_immob,N_imm(5),Nfix0
        real N_transfer, N_loss,Qroot0,Cfix0
        real Scalar_N_flow,Scalar_N_T, kappaVcmax
        real SNfine,SNcoarse,SNmicr,SNslow,SNpass
        real costCuptake,costCfix,costCreuse
        real Creuse0,Nup0,N_deN0,LDON0
        !     the variables relative to soil moisture calcualtion
        real S_omega    !  average values of the moisture scaling functions
        real S_t(5)     !  average values of temperature scaling functions
        real S_w_min    !  minimum decomposition rate at zero available water
        real ScNloss
        integer i,j,n
        real frac_soc(10)
        ! ------ outputs? ------- !
        ! real Rnitrogen,SNRauto ! SNRauto to respiration
        real testout(11), ksye !,tsoil_layer(11)
                
        tsoil_layer = testout
        frac_soc    = (/0.75,0.2,0.02,0.015,0.005,0.0,0.0,0.0,0.0,0.0/)
        Q10h        = (/Q10,Q10,Q10,Q10,Q10/) ! temperature sensitivity of Rh_pools
        Qroot0      = 500.
        Nfix0       = 1./60.   ! maximum N fix ratio, N/C
        Nup0        = 0.02     ! nitrogen uptake rate
        Cfix0       = 12.      ! C cost per N for fixation
        ksye        = 0.05      ! C cost per N for uptake
        Creuse0     = 2.     ! C cost per N for resorption
        ScNloss     = 1.
        N_deN0      = 1.E-3*ScNloss   ! 1.E-3, 5.E-3, 10.E-3, 20.E-3
        LDON0       = 1.E-3*ScNloss
        Rnitrogen   = 0.
        ! for N scalars
        CNmin       = 40.0
        CNmax       = 200.0
        ! Max and min NSN pool
        ! NSNmax      = QN(1) + 0.2*QN(2) + QN(3)  ! 15.0
        ! NSNmin      = 0.01
        ! partitioning coefficients
        etaL        = 0.6          ! 60% of foliage litter is fine, didn't use
        etaW        = 0.15        ! 15% of woody litter is fine
        etaR        = 0.85         ! 85% of root litter is fine  , didn't use    
        f_F2M       = 0.55        ! *exp((CN0(4)-CN_fine)*0.1)
        f_C2M       = 0.275       ! *exp((CN0(5)-CN_coarse)*0.1)
        f_C2S       = 0.275       ! *exp((CN0(5)-CN_coarse)*0.1)
        f_M2S       = 0.3
        f_M2P       = 0.1
        f_S2P       = 0.2        !0.03 Change by Jiang Jiang 10/10/2015
        f_S2M       = 0.5
        f_P2M       = 0.45
        ! calculating soil scaling factors, S_omega and S_tmperature
        S_w_min     = 0.08 !minimum decomposition rate at zero soil moisture
        S_omega     = S_w_min + (1.-S_w_min) * Amin1(1.0, 0.3*omega)
        ! ------ commented for CWE
        do i=1,5
            S_t(i) = Q10h(i)**((Tsoil-10.)/10.)  ! Duke
        enddo
        if (doSoilphy) then 
            S_t  = (/0.0,0.0,0.0,0.0,0.0/)
            do i = 1,5
                if(i.lt.3) then    ! couarse and fine litter use surface layer soil temperature
                    S_t(i) = Q10h(i)**((tsoil_layer(1)-10.)/10.)  ! Duke
                else 
                    do j = 1,10       ! fast,slow and passive pool use weighed soil temperature in layers according to soc distribution
                        S_t(i) = S_t(i)+frac_soc(j)*Q10h(i)**((tsoil_layer(j+1)-10.)/10.)  ! Duke
                    enddo
                endif
            enddo
        else
            do i=1,5
                S_t(i) = Q10h(i)**((Tsoil-10.)/10.)  ! Duke
            enddo  
        endif 
        ! Calculating NPP allocation and changes of each C pool
        NPP_L       = alpha_L * NPP           ! NPP allocation
        NPP_W       = alpha_W * NPP
        NPP_R       = alpha_R * NPP
        ! N scalars on decomposition
        SNfine      = exp(-(CN0(4)-CN(4))/CN0(4)) 
        SNcoarse    = exp(-(CN0(5)-CN(5))/CN0(5)) 
        SNmicr      = exp(-(CN0(6)-CN(6))/CN0(6)) 
        SNslow      = 1.    !exp(-(CN0(7)-CNC(7))/CN0(7)) 
        SNpass      = exp(-(CN0(8)-CN(8))/CN0(8)) 
        
        ! the carbon leaving the pools
        OutC(1)     = L_fall
        OutC(2)     = QC(2)/tauC(2)*S_omega !*exp(CN(2)/CN0(2)-1.) 
        OutC(3)     = QC(3)/tauC(3)*S_omega
        OutC(4)     = QC(4)/tauC(4)*S_omega* S_T(1)*CN(4)/CN0(4)!*SNfine
        OutC(5)     = QC(5)/tauC(5)*S_omega* S_T(2)*CN(5)/CN0(5)!*SNcoarse
        OutC(6)     = QC(6)/tauC(6)*S_omega* S_T(3)!*SNmicr
        OutC(7)     = QC(7)/tauC(7)*S_omega* S_T(4)!*SNslow
        OutC(8)     = QC(8)/tauC(8)*S_omega* S_T(5)!*SNpass
        ! heterotrophic respiration from each pool
        Rh_pools(1) = OutC(4)* (1. - f_F2M)
        Rh_pools(2) = OutC(5)* (1. - f_C2M - f_C2S)
        Rh_pools(3) = OutC(6)* (1. - f_M2S - f_M2P)
        Rh_pools(4) = OutC(7)* (1. - f_S2P - f_S2M)
        Rh_pools(5) = OutC(8)* (1. - f_P2M)

        !========================================================================
        ! Nitrogen part
        ! nitrogen leaving the pools and resorption
        do i=1,8
            OutN(i) = OutC(i)/CN(i)
        enddo
        ! nitrogen mineralization
        N_miner = OutN(4)* (1. - f_F2M)  &
                    &   +OutN(5)* (1. - f_C2M - f_C2S) &
                    &   +OutN(6)* (1. - f_M2S - f_M2P) &
                    &   +OutN(7)* (1. - f_S2P - f_S2M) &
                    &   +OutN(8)* (1. - f_P2M)
        ! Nitrogen immobilization
        N_imm   = 0.
        N_immob = 0.
        if(QNminer>0)then
            do i=4,8
                if(CN(i)>CN0(i))then
                    N_imm(i-3) = Amin1(QC(i)/CN0(i)-QC(i)/CN(i),0.1*QNminer)
                    N_immob    = N_immob+N_imm(i-3)
                endif
            enddo
        endif
        ! Let plant itself choose the strategy between using C to uptake
        ! or fix N2 by comparing C invest.
        ! N demand
        N_demand    = NPP_L/CN0(1)+NPP_W/CN0(2)+NPP_R/CN0(3) !+N_deficit
        ! Nitrogen input:
        N_transfer  = 0.
        N_uptake    = 0.
        N_fixation  = 0.
        costCuptake = 0.
        costCfix    = 0.
        costCreuse  = 0.
        ! 1. Nitrogen resorption
        N_transfer  = (OutN(1) + OutN(2) +OutN(3))*alphaN
        costCreuse  = Creuse0*N_transfer
        N_demand    = N_demand-N_transfer
        If(N_demand>0.0)then
            ! 2.  N uptake
            if(ksye/QNminer<Cfix0)then
                N_uptake   = AMIN1(N_demand+N_deficit,          &
                             &   QNminer*QC(3)/(QC(3)+Qroot0), &
                             &   Nup0*NSC/(ksye/QNminer))
                costCuptake = N_uptake*(ksye/QNminer)
                N_demand    = N_demand-N_uptake
            elseif(NSN<24.*30.*N_demand)then
                ! 3.  Nitrogen fixation
                N_fixation  = Amin1(N_demand,fnsc*Nfix0*NSC)
                costCfix    = Cfix0*N_fixation
                N_demand    = N_demand-N_fixation
            endif
        endif
        N_deficit = N_deficit+N_demand
        ! update NSN
        NSN       = NSN+N_transfer+N_uptake+N_fixation
        ! write(*,*)"test-NSN:", NSN, N_transfer, N_uptake, N_fixation
        ! Total C cost for nitrogen
        Rnitrogen = costCuptake+costCfix+costCreuse
        ! Oak Ridge N fixation rate: 
        ! asymbiotic: 2 mg N/m2/yr ;  symbiotic: 65 mg/m2/yr, Oak Ridge
        ! N_fixation=0.067/8760. ! Oak Ridge
        ! N_fixation=0.23/8760.  ! Duke
        ! Nitrogen using, non-structural nitrogen pool, NSN
        N_leaf    = AMIN1(NPP*alpha_L/CN(1)+QC(1)/CN0(1)-QC(1)/CN(1),0.2*NSN)
        N_wood    = AMIN1(NPP*alpha_W/CN(2)                         ,0.1*NSN)
        N_root    = AMIN1(NPP*alpha_R/CN(3)+QC(3)/CN0(3)-QC(3)/CN(3),0.2*NSN)
        ! write(*,*)"test-NRoot:", NPP, alpha_R, CN(3), QC(3), CN0(3), NSN
        NSN       = NSN-(N_leaf+N_wood+N_root)
        N_LF      = OutN(1)*(1.-alphaN)
        N_WF      = OutN(2)*(1.-alphaN)
        N_RF      = OutN(3)*(1.-alphaN)
        ! update QNminer
        QNminer   = QNminer+N_miner+N_deposit  &
                    &   -(N_uptake+N_immob)
        ! Loss of mineralized N and dissolved organic N
        Scalar_N_flow = 0.5*runoff/rdepth
        ! Scalar_N_T=0.005*(Tsoil+273.)/(Tsoil+273+333.)
        ! commented line for soil thermal       
        Scalar_N_T    = N_deN0*exp((Tsoil-25.)/10.)
        ! added lines for soil thermal
        if (doSoilphy) then 
            Scalar_N_T = 0.0 
            do j=1,10
                Scalar_N_T = Scalar_N_T + frac_soc(j)*N_deN0*exp((tsoil_layer(j+1)-25.)/10.)  
            enddo
        else
            Scalar_N_T     = N_deN0*exp((Tsoil-25.)/10.)
        endif  
        N_leach = Scalar_N_flow*QNminer+Scalar_N_flow*QN(6)*LDON0
        N_vol   = Scalar_N_T*QNminer
        N_loss  = N_leach + N_vol
        ! update QNminer
        QNminer = QNminer-N_loss
        ! update plant carbon pools, ! daily change of each pool size
        QC(1)   = QC(1) - OutC(1) + NPP_L    
        QC(2)   = QC(2) - OutC(2) + NPP_W
        QC(3)   = QC(3) - OutC(3) + NPP_R
        QC(4)   = QC(4) - OutC(4) + OutC(1)+etaW*OutC(2)+OutC(3)
        QC(5)   = QC(5) - OutC(5) + (1.-etaW)*OutC(2)
        QC(6)   = QC(6) - OutC(6) + f_F2M*OutC(4)+f_C2M*OutC(5)     &         
                    &   + f_S2M*OutC(7)+f_P2M * OutC(8)
        QC(7)   = QC(7) - OutC(7)+f_C2S*OutC(5)+f_M2S*OutC(6)
        QC(8)   = QC(8) - OutC(8)+f_M2P*OutC(6)+f_S2P*OutC(7)
        Q_plant = QC(1) + QC(2) + QC(3)
        ! update nitrogen pools
        QN(1)   = QN(1) - OutN(1) + N_leaf
        QN(2)   = QN(2) - OutN(2) + N_wood
        QN(3)   = QN(3) - OutN(3) + N_root
        QN(4)   = QN(4) - OutN(4)+ N_imm(1)     &
                    &   + (OutN(1) + etaW*OutN(2) + OutN(3))*(1.-alphaN)
        QN(5)   = QN(5) - OutN(5) + N_imm(2)   &
                    &   + (1.-etaW)*OutN(2)*(1.-alphaN)
        QN(6)   = QN(6) - OutN(6) + N_imm(3) - Scalar_N_flow*QN(6)*LDON0  &
                    &   + f_F2M*OutN(4)+f_C2M*OutN(5)  &
                    &   + f_S2M*OutN(7)+f_P2M*OutN(8)
        QN(7)   = QN(7) - OutN(7) + N_imm(4)  &
                    &   + f_C2S*OutN(5) + f_M2S*OutN(6)
        QN(8)   = QN(8) - OutN(8) + N_imm(5) &
                    &   + f_M2P*OutN(6) + f_S2P*OutN(7)
        QNplant = QN(1) + QN(2)+ QN(3)
        ! update C/N ratio
        CN         = QC/QN
        CN_foliage = (QC(1)+QC(3))/(QN(1)+QN(3))
        ! calculate N related scalars for Oak Ridge
        ! SNvcmax =exp(-(CN(1)-CN0(1))) ! /CN0(1) ! Oak Ridge
        ! SNgrowth=exp(-(CN(1)-CN0(1))/CN0(1)) !  AMAX1((CNmax-CN_foliage)/(CNmax-CNmin),0.0)+0.25
        ! SNRauto =AMAX1((CNmax-CN_foliage)/(CNmax-CNmin),0.0)+0.5
        ! SNrs=1.
        ! calculate N related scalars for Duke FACE
        kappaVcmax = CN0(1)/1.
        SNvcmax    = exp(-kappaVcmax*(CN(1)-CN0(1))/CN0(1)) ! /CN0(1) ! Duke
        SNgrowth   = exp(-(CN(1)-CN0(1))/CN0(1)) !  AMAX1((CNmax-CN_foliage)/(CNmax-CNmin),0.0)+0.25
        SNRauto    = exp(-(CN(1)-CN0(1))/CN0(1)) !  AMAX1((CNmax-CN_foliage)/(CNmax-CNmin),0.0)+0.5
        ! SNrs       = 1.
        ! write(*,*)"test-QC2:", L_fall
        ! write(*,*)"test-SNRauto:", CN(1),CN0(1), CN0(1), QC(1), QN(1)
        ! write(*,*)"test-QN-1", QN(1), OutN(1), N_leaf
        ! write(*,*)"test-OUTN(1)",OutC(1), CN(1), NPP, alpha_L, NSN
        ! write(*,*)"test-nsn:",NSN, N_leaf, N_wood, N_root
        ! write(*,*)"test-CN(3:", qc(3), qn(3)
        ! write(*,*)"test-QN3:", OutN(3), N_root
        ! WRITE(*,*)"test-outn3:", OutC(3), CN(3)
        ! write(*,*)"test-outC3:", QC(3), tauC(3), S_omega
        ! write(*,*)"etst-s_omega", S_w_min, S_w_min, 0.3*omega
        return
    end subroutine TCS_CN
    
end module mod_transferCpools