module mod_methane
    use mod_dataTypes
    implicit none

    ! ---------------------------------------------------------------------------
    contains
    subroutine methane()
        ! update single value in a hourly loop when MEMCMC=0
        ! update single value of Rh_pools,Tsoil,zwt,wsc in a hourly loop when MEMCMC=1
        implicit none 
        ! introduce variables and constants used in this subroutine; set soil layers

        integer i
        integer,parameter :: nlayers=10       !use this statement to set the parameter value
        real zwt    
        real consum
        ! set values for MEMCMC
        integer,parameter :: miterms=17
        integer,parameter :: ilines=9000  
        ! CH4 Production    
        real Rh(nlayers),Rh_pools(5),Rh_h,ProCH4(nlayers),Pro_sum
        real r_me         !release ratio of CH4 to CO2
        real Q10pro
        real fSTP(nlayers)         !CH4 production factor of soil temperature
        real vt,xt
        real Tmax_me,Tpro_me
        real fpH          !CH4 production factor of soil pH
        real fEhP         !CH4 production factor of soil redox potential
        real FRLEN(nlayers)        !fraction of root in each layer
        real Tsoil
        ! CH4 Oxidation
        real CH4(nlayers),CH4_V(nlayers)          !both are CH4 concentration: CH4(nlayers)unit gC/m2, CH4_V(nlayers) unit g C/ m3
        real wsc(nlayers)      
        real OxiCH4(nlayers),Oxi_sum       !CH4 oxidation
        real Omax_layers(nlayers),Omax       !maximum oxidation rate
        real kCH4_layers(nlayers),kCH4       !system specific half saturation constant
        real Q10oxi
        real fCH4(nlayers)         !CH4 oxidation factor of CH4 concentration
        real fSTO(nlayers)         !CH4 oxidation factor of soil temperature
        real fEhO         !CH4 oxidation factor of soil redox potential
        real Toxi
        ! CH4 Diffusion
        real Deff(nlayers)     !CH4 effective diffusivity !!used in mineral soil  v1.1 
        real D_CH4_a           !CH4 diffusion coefficient in air  !unit cm2 s-1   diffusivity of CH4 in air
        real D_CH4_w           !CH4 diffusion coefficient in water  !unit cm2 s-1   diffusivity of CH4 in water
        real phi          !soil porosity  also used in water table part
        real fwater(nlayers),fair(nlayers)
        real D_CH4_soil(nlayers),D_CH4_soil_a(nlayers),D_CH4_soil_b(nlayers)      !!used in organic peat soil  v1.2
        real fcoarse      !relative volume of coarse pores depending on soil texture  Zhuang 2004
        real ftort        !tortuosity coefficient with a value of 0.66    Walter and Heimann 2000
        !suggesting that the distance covered by diffusion is about two thirds of the length of the real average path
        real SAND         !relative contents of sand (%) in the soil
        real PVSAND       !relative volume of coarse pores in sandy soils     set to 0.45     value from Walter 2001
        real SILT         !relative contents of silt (%) in the soil
        real PVSILT       !relative volume of coarse pores in silty soils     set to 0.20     value from Walter 2001
        real CLAY         !relative contents of clay (%) in the soil
        real PVCLAY       !relative volume of coarse pores in clayish soils     set to 0.14   value from Walter 2001
        real DEPTH(10)        !depth in soil  will define it inside this subroutine again      resolution 100mm 200mm
        real THKSL(10)        !will define it inside this subroutine again  
        ! real Fdifu(nlayers+1)
        real Fdifu(nlayers)
        real CH4_atm      !concentration of CH4 in atmosphere     seen as 0 cause the value is too low someone use 0.076
        real simuCH4      !simulated CH4 emission
        ! Boundary condition parameters
        real ScCH4                                 !Schmidt numbers for methane Wania
        real pistonv                               !Piston velocity
        real Ceq                                   !equilibrium concentration of gas in the atmosphere
        real kHinv                                 !Henry's coefficient dependent variable on left side of equation, T is the independent variable
        real kH_CH4         !Henry's constant at standard temperature (CH4) Unit L atm mol-1
        real CHinv          !Coefficient in Henry's Law Unit K      
        real Tsta           !standard temperature Unit K
        real Ppartial       !CH4 partial pressure in air Unit atm
        ! Ebullition 
        real CH4_thre,CH4_thre_ly(nlayers),EbuCH4(nlayers),Kebu
        real Ebu_sum_unsat,Ebu_sum_sat,Ebu_sum          !sum value one dimension is enough 
        integer wtlevelindex
        ! Plant transport
        real PlaCH4(nlayers),Pla_sum
        real LAIMIN,LAIMAX
        real Tveg,Tgr,Tmat,fgrow,Pox,Kpla
        ! Yuan added for soil temp  
        logical doSoilphy
        real testout(11), tsoil_layer(11)
        ! MEMCMC=0   ! note here, any changes here result unexpected bug 
        Rh_h        = Rh_pools(1)+Rh_pools(2)+Rh_pools(3)+Rh_pools(4)+Rh_pools(5)  !hourly Rh_f + Rh_c + Rh_Micr + Rh_Slow + Rh_Pass
        tsoil_layer = testout
        FRLEN       = (/0.75,0.2,0.02,0.015,0.005,0.0,0.0,0.0,0.0,0.0/)             
        thksl       = (/10.,10.,10.,10.,10.,20.,20.,20.,20.,20./)
        simuCH4     = 0.0                 ! v1.2 
        do i = 1, nlayers   ! put it out of the subroutine
            ! Rh weighed according to the distribution of root
            if (i .LE. 3) then                                 ! the empirical method used here is from CLM4.5
                Rh(i)= 0.5*Rh_h*FRLEN(i)+((0.5*Rh_h)/0.3)*0.1  ! Rh(h,i)Rh produced by each layer per hour  unit should be g C m-2 h-1 
            else                                               ! i*10: depth of ith soil layers
                Rh(i)= 0.5*Rh_h*FRLEN(i)
            endif
        enddo   
        Tmax_me=45.0
        do i = 1,nlayers          
            if (doSoilphy) then
                if (tsoil_layer(i+1) .lt. 0.0) then
                    fSTP(i) = 0.0
                else if (tsoil_layer(i+1) .gt. Tmax_me) then
                    fSTP(i) = 0.0
                else if (tsoil_layer(i+1) .ge. 0.0 .and. tsoil_layer(i) .le. Tmax_me) then
                    fSTP(i) = Q10pro**((tsoil_layer(i+1)-Tpro_me)/10)        !Tsoil is the only variable
                endif
            else 
                if (Tsoil .lt. 0.0) then
                    fSTP(i) = 0.0
                else if (Tsoil .gt. Tmax_me) then
                    fSTP(i) = 0.0
                else if (Tsoil .ge. 0.0 .and. Tsoil .le. Tmax_me) then
                    fSTP(i) = Q10pro**((Tsoil-Tpro_me)/10)        !Tsoil is the only variable
                endif
            endif
        enddo
        fpH      = 1.0
        fEhP     = 1.0
        depth(1) = 10.0                                  !calculate soil depth unit cm
        do i=2,nlayers
            depth(i)=depth(i-1)+THKSL(i)
        enddo
        Pro_sum = 0.0
        do i = 1,nlayers
          if ((depth(i)*10) .le. -zwt) then
                ProCH4(i)=0.0
          else
              if (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then
                  ProCH4(i) = Rh(i)*r_me*fSTP(i)*fpH*fEhP*(((depth(i)*10.0)-(-zwt))/(THKSL(i)*10.0))     ! *percent
              elseif (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then
                  ProCH4(i) = Rh(i)*r_me*fSTP(i)*fpH*fEhP
              endif
          endif
          Pro_sum = Pro_sum+ProCH4(i)
        enddo
        ! Add CH4 production to CH4 pool    (gC layer -1)=(gC m-2)
        do i=1,nlayers
            CH4(i) = CH4(i) + ProCH4(i)
        enddo
        ! END OF METHANE PRODUCTION
        ! B. methane oxidation      hourly  unit gC m-2 h-1     !!!!!!!!!!!method of CLM and Zhuang!!!!!!!!!!!!
        ! Methane oxidation is modeled as an aerobic process that occurs in the unsaturated zone of the soil profile ZHUANG
        ! fSTO assignment      
        Q10oxi  = 2.0      !Zhu 2014 results from previous studies  unit 1  also used by zhang
        do i=1,nlayers
            if (doSoilphy) then
                fSTO(i)=Q10oxi**((tsoil_layer(i+1)-Toxi)/10.0)
            else
                fSTO(i)=Q10oxi**((Tsoil-Toxi)/10.0)
            endif
        enddo
        ! fEhO assignment
        fEhO    = 1.0        !Walter 2000  did not consider it, equal to value of 1
        Oxi_sum = 0.0
        do i = 1,nlayers
            Omax_layers(i) = (Omax/(1000000))*12*1000*(wsc(i)*0.001)     !convert the unit of Omax from μmol L-1 h-1 to gC m-2 h-1
            kCH4_layers(i) = (kCH4/(1000000))*12*1000*(wsc(i)*0.001)    !convert the unit of kCH4 from μmol L-1 to gC m-2
            ! then calculate fCH4 with CH4(i) and kCH4_layers(i) 
            fCH4(i)=CH4(i)/(kCH4_layers(i)+CH4(i))   !  CH4 concentration factor
            if ((depth(i)*10.0) .le. -zwt) then                !unit of Omax: gC m-2 h-1
                OxiCH4(i) = Omax_layers(i)*fCH4(i)*fSTO(i)*fEhO!*0.1      !wrong:*(THKSL(i)/1000)!mm to m account for the thickness
            else
                if (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then
                    if (i .eq. 1) then
                        OxiCH4(i)=Omax_layers(i)*fCH4(i)*fSTO(i)*fEhO*((-zwt)/(THKSL(i)*10.0))
                    else
                        OxiCH4(i)=Omax_layers(i)*fCH4(i)*fSTO(i)*fEhO*(((-zwt)-(depth(i-1)*10.0))/(THKSL(i)*10.0))      !  *percent
                    endif
                else if (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then
                    OxiCH4(i)= 0.0
                endif
            endif        
            if (OxiCH4(i) .gt. CH4(i)) then
                OxiCH4(i) = CH4(i)
            endif
            Oxi_sum=Oxi_sum+OxiCH4(i)  
        enddo 
        ! minus CH4 oxidation from CH4 pool
        do i=1,nlayers
            CH4(i)   = CH4(i) - OxiCH4(i)               !minus CH4 oxidation from CH4 pool
            CH4_V(i) = CH4(i)/(wsc(i)*0.001)          !convert concentration from gC/m2 to gC/m3; !CH4_V(i) can be used for DA with observation data in soil layers                                             
        enddo
        ! END OF METHANE OXIDATION
        ! C. methane diffusion
        ! Parameters assignment 
        D_CH4_a = 0.2            !unit cm2 s-1   D_CH4_a is the molecular diffusion coefficient of methane in air
        D_CH4_a = (D_CH4_a/10000.0)*3600.0        !unit m2 h-1
        D_CH4_w = 0.00002        !unit cm2 s-1   D_CH4_a is the molecular diffusion coefficient of methane in water
        D_CH4_w = (D_CH4_w/10000.0)*3600.0        !unit m2 h-1          
        ftort   = 0.66        !tortuosity coefficient with a value of 0.66    Walter and Heimann 2000
        ! parameters for fcoarse algorithm      
        SAND    = 0.4             !   %   SPRUCE site value    0.4
        SILT    = 0.4             !   %   SPRUCE site value   0.4
        CLAY    = 0.2             !   %   SPRUCE site value   0.2
        PVSAND  = 0.45       !relative volume of coarse pores in sandy soils       set to 0.45     value from Walter 2001 zhuang
        PVSILT  = 0.20       !relative volume of coarse pores in silty soils       set to 0.20     value from Walter 2001 zhuang
        PVCLAY  = 0.14       !relative volume of coarse pores in clayish soils     set to 0.14     value from Walter 2001 zhuang  
        fcoarse = SAND*PVSAND+SILT*PVSILT+CLAY*PVCLAY
        CH4_atm = 0.076       !unit umol L-1
        ! Peat soil solution for diffusion coefficient: Equations for D_CH4_soil *         v1.2    Millington and Quirk Model
        do i=1,nlayers
            fwater(i)       = wsc(i)/(THKSL(i)*10)      
            fair(i)         = phi-fwater(i)
            D_CH4_soil_a(i) = (((fair(i))**(10/3))/((phi)**2))*D_CH4_a
            D_CH4_soil_b(i) = D_CH4_W
            if (fair(i) .ge. 0.05) then
                D_CH4_soil(i) = D_CH4_soil_a(i)
            else
                D_CH4_soil(i) = D_CH4_soil_b(i)
            endif      
            Deff(i) = D_CH4_soil(i)
        enddo 
        ! Mineral soil solution for diffusion coefficient: Equations for D_CH4_soil v1.1   Three-porosity-model
        CH4_atm  = (CH4_atm/1000000)*12*1000   
        kH_CH4   = 714.29
        CHinv    = 1600.0
        Tsta     = 298.15
        Ppartial = 1.7E-20 
        ScCH4    = 1898 - 110.1*Tsoil + 2.834*Tsoil**2 - 0.02791*Tsoil**3
        pistonv  = 2.07 * (ScCH4/600)**(-1/2)
        kHinv    = kH_CH4 /((exp(CHinv*(1/(Tsoil+273.15)-1/Tsta))))
        Ceq      = Ppartial / kHinv    ! Ceq: mol L-1   p_partial: atm  kHinv：L atm mol-1    
        Fdifu(1) =  pistonv * (CH4_V(1) - Ceq) !refer to the interface of methane flux from layer 1 to atmosphere  cm to m   switch on/off
        do i = 2,nlayers                                  !refer to flux from layer ii to ii-1 
            Fdifu(i)= Deff(i)*(CH4_V(i)-CH4_V(i-1))/(THKSL(i)*0.01)      !the unit of Fdifu is gC/m-2/h
        enddo     
        ! below I try to keep the CH4 flux no larger than the amount of CH4 that exist at the moment   V1.1 V1.2
        do i=1,nlayers+1
            if (Fdifu(i) .gt. 0.0 .and. (Fdifu(i)) .gt. CH4(i)) then
                Fdifu(i)=CH4(i)
            else if (Fdifu(i) .lt. 0.0 .and. (abs(Fdifu(i))) .gt. CH4(i-1)) then
                Fdifu(i)=-CH4(i-1)    
            endif
        enddo
        do i = 1,nlayers-1                                  !loop of time
            CH4(i) = CH4(i) + (Fdifu(i+1)-Fdifu(i))*1 ! *1   * 1 hour /hour   /15min  *0.25h 
            if (CH4(i) .lt. 0.0) then                     ! this part need to be improved until deleted   V1.2
                CH4(i) = 0.0
            endif
        enddo    
        CH4(10) = CH4(10) - Fdifu(10)                                   !MODIFIED ON 07/25/2016
        if (CH4(10) .lt. 0.0) then                                    !defined the Fdifu(11) to be 0.0
            CH4(10)= 0.0                                              ! switch on/off
        endif
        simuCH4 = simuCH4 + (Fdifu(1)-0.0) 
        !     ********************************************************************************************************************      
        ! D. methane ebullition     !assume bubbles can reach the water table within 1 h&
                                    !& the bubbles is added to the methane concentration in the soil layer just above the wt
                                    !& and then diffused through layers   ??not correct
        ! this subroutine is modified on 02132017 by deleting the unsat from bubble and add unsat to concentration so as to increase diffusion
        ! just by searching "switch" you can switch from old to new mode by adding or deleting "!"
        ! modified threshold value to 100 for testing
        !     ********************************************************************************************************************
        Kebu          = 1.0                    !unit  h-1   rate constant          
        Ebu_sum_unsat = 0.0
        Ebu_sum_sat   = 0.0                                      !initial value
        
        do i=1,nlayers
            CH4_thre       = 1000.0  !!find in parafile  !unit  umol L-1 according to Walter's 500-1000
            CH4_thre_ly(i) = (CH4_thre*1.0e-6)*12*1000*(wsc(i)*0.001)    !convert the unit of CH4_thre from µmol L-1 to gC m-2
        enddo
        if (zwt .ge. 0.0) then                                  !when water table is above the soil surface
            do i=1,nlayers
                if (CH4(i) .gt. CH4_thre_ly(i)) then                  
                        EbuCH4(i)=Kebu*(CH4(i)-CH4_thre_ly(i))     !only if the concentration is larger than threshold
                else !if (CH4(i) .le. CH4_thre_ly(i)) then
                        EbuCH4(i)=0.0
                endif
                Ebu_sum_sat = Ebu_sum_sat+EbuCH4(i)               !& the bubbles are directly added into CH4 efflux into atmosphere
                CH4(i)      = CH4(i)- EbuCH4(i)                        !& update the concentration at the end of this hour in each layers
            enddo
        endif
        if (zwt .lt. 0.0) then                                  !when water table is below the soil surface
            do i=1,nlayers
                if ((depth(i)*10.0) .le. -zwt) then               !acrotelm layers
                    EbuCH4(i)     = 0.0
                    Ebu_sum_unsat = Ebu_sum_unsat+EbuCH4(i)         
                    CH4(i)        = CH4(i)- EbuCH4(i) 
                else
                    if (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then       !partly acrotelm layer
                        wtlevelindex = i
                        if (CH4(i) .gt. CH4_thre_ly(i)) then                  
                        EbuCH4(i) = Kebu*(CH4(i)-CH4_thre_ly(i))!*(((depth(i)*10.0)-(-zwt))/(THKSL(i)*10.0))        ! * percent
                        else !if (CH4(i) .le. CH4_thre_ly(i)) then                     ??????????,??????????????????
                        EbuCH4(i) = 0.0
                        endif 
                        CH4(i)              = CH4(i)- EbuCH4(i)
                        Ebu_sum_unsat       = Ebu_sum_unsat+EbuCH4(i)                ! !  modified by Mary on 02132017
                        CH4(wtlevelindex-1) = CH4(wtlevelindex-1)+EbuCH4(i)    !!!!!-1-!!!! !switch on in new mode should be added add burst bubbles below surface to diffusion modified by Mary on 02132017
                        ! ：the problem is the resolution of soil layer is 10cm and EbuCH4(i) is directly added to the upper layer of boundary layer   02152017  
                    else if (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then   !catotelm layers
                        if (CH4(i) .gt. CH4_thre_ly(i)) then                  
                            EbuCH4(i) = Kebu*(CH4(i)-CH4_thre_ly(i))
                        else !if (CH4(i) .le. CH4_thre_ly(i)) then
                            EbuCH4(i) = 0.0
                        endif 
                        CH4(i)              = CH4(i)- EbuCH4(i)
                        CH4(wtlevelindex-1) = CH4(wtlevelindex-1)+EbuCH4(i)     !!!!!-2-!!!! !switch on in new mode should be added     modified by Mary on 02152017
                        Ebu_sum_unsat       = Ebu_sum_unsat+EbuCH4(i)                  ! modified by Mary on 02132017
                    endif
                endif
            enddo
        endif
        Ebu_sum = Ebu_sum_sat
        ! E. plant mediated methane transportation      totoally used Walter's model also used by Zhuang et. al
        Kpla=0.01         !unit h-1
        ! the Tsoil used here would be better if refer to the 20cm soil temperature after &
        ! & the accomplishment of soil heat dynamics module. according to Zhuang. however Walter used 50cm soil temp.
        Tgr  = 2.0               !unit degree Celsius if annual mean temp is below 5 (otherwise 7)
        Tmat = Tgr+10.0         !unit degree Celsius
        Pox  = 0.5               !50% of mediated methane are oxidised 
        if (Tsoil .lt. Tgr) then ! define fgrow
            fgrow = LAIMIN
        else if (Tsoil .ge. Tgr .and. Tsoil .le. Tmat) then
            fgrow = LAIMIN+LAIMAX*(1-((Tmat-Tsoil)/(Tmat-Tgr))**2)
        else if (Tsoil .gt. Tmat) then
            fgrow = LAIMAX
        endif
        Pla_sum=0.0
        do i=1,nlayers
            PlaCH4(i) = Kpla*Tveg*FRLEN(i)*fgrow*CH4(i)*(1-Pox)         !not sensitive at all to this change, but better
            Pla_sum   = Pla_sum+PlaCH4(i)
            CH4(i)    = CH4(i)-PlaCH4(i)
            CH4_V(i)  = CH4(i)/(wsc(i)*0.001)          !convert concentration from gC/m2 to gC/m3
        enddo
        simuCH4 = simuCH4+Pla_sum
        consum  = simuCH4+OxiCH4(1)+OxiCH4(2)+OxiCH4(3)+OxiCH4(4)+OxiCH4(5)+OxiCH4(6)+OxiCH4(7)+OxiCH4(8)+OxiCH4(9)+OxiCH4(10)
      return
    end subroutine methane

end module mod_methane