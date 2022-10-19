!!
module mod_soil
    use mod_dataTypes
    ! use mod_vegetation
    implicit none
    
    contains
    subroutine soil_thermal(Rsoilab1, Rsoilab2, QLleaf, QLair, FLAIT, raero, & ! input from canopy
                            & Esoil, Hsoil) ! may be out
        implicit none
        real Rsoilab1, Rsoilab2, QLleaf, QLair ! input from canopy module
        real FLAIT, raero                      ! input from canopy module
        ! real rhoS(3), Tsnow, Twater, Tice, ice_tw, diff_s, diff_snow, liq_water, dcount_soil        ! this is const paramter and update in canopy module
        ! ice is used to calculate, then it is updated. zhou: must be test. 
        ! must be define in the driver module
        ! real snow_depth, water_tw    ! initial in main, update in this part
        ! Esoil, Hsoil is updated in this part, but seem not to be inputed?
        ! real zwt ! is updated in soil water processes
        ! ----- local var -----!
        real,dimension(11):: testout    ! maybe just for test?    

        integer n_layers, i
        real ice_density, thkns1, shcap_ice, condu_ice, condu_water, shcap_water
        real condu_soil, shcap_soil, condu_s, thd_t, latent_heat_fusion, condu_air
        real shcap_air, diff_air, water_table_depth, snow_depth_t
        real albedo_water, WILTPT, FILDCP, flux_snow
        real ufw(10),frac_ice1, frac_ice2, QLsoil, Rsoilab3
        real Rsoilabs 
        real rhocp,slope,psyc,Cmolar, H2OLv,fw1
        real Rsoil, Esoil
        real rLAI, resht_lai, hsoil
        real condu(10), shcap(10), difsv1, resdh
        real dnr,dsh,dgh,dle,drsdh, delta, sftmp_pre, tsoill_0
        real tsoill_pre, thkns2, difsv2,temph1,temph2 
        real heat_adjust, heat_excess, ice_incr, inter_var, temph_snow
        real temph_water
        
        real, allocatable ::depth_z(:) 
        n_layers = 10
        allocate(depth_z(n_layers))
        ice_density = 916.                  !   916. soil thermal conductivity W m-2 K-1
        thkns1      = thksl(1)/4.
        shcap_ice   = 2117.27*ice_density
        condu_ice   = 2.29
        condu_water = 0.56!0.56
        shcap_water = 4188000.
        condu_soil  = 0.25
        shcap_soil  = 2600000.
        condu_s     = 0.25
        thd_t       = -1.0
        diff_snow   = 3600.*condu_snow/shcap_snow*10000.
        diff_s      = 3600.*condu_b/shcap_soil*10000.
        latent_heat_fusion = 333700.   ! j kg-1
        condu_air          = 0.023
        shcap_air          = 1255.8
        diff_air    = 3600.*condu_air/shcap_air*10000. 
        water_tw    = zwt*0.001 - ice_tw         ! might means total water that is liquid, add up all layers
        water_table_depth = zwt*0.1
        snow_depth_t = snow_depth - 0.46*0.0     ! warming in Tair impact on snow_depth
                                                 ! in unit cm 0.46 based on snow_depth vs. tair regression     
        if (snow_depth_t .lt. thd_snow_depth) snow_depth_t =0.0
        if (snow_depth_t .gt. 0.) then
            dcount_soil = dcount_soil +1./24.
        else 
            dcount_soil = 0.
        endif

        if (water_table_depth .lt. 4. .and. water_table_depth .gt. 0.0) water_table_depth =0.    ! avoid numerical issues when      
        albedo_water = 0.1      
        ! soil water conditions
        WILTPT       = wsmin/100.
        FILDCP       = wsmax/100. 
        flux_snow    = 0.0 
        depth_z      = (/0., 0., 0., 0., 0., 0., 0.,0.,0.,0./) 
        ufw          = (/0.0163,0.0263,0.0563,0.0563,0.0563,0.1162,0.1162,0.1162,0.1162,0.1162/)   !  ..int add unfrozen water ratio
        frac_ice1    = 0.01         !0.015
        frac_ice2    = 0.001        !0.01
      
        QLsoil       = emsoil*sigma*((sftmp+273.2)**4)
        Rsoilab3     = (QLair+QLleaf)*(1.0-rhoS(3))-QLsoil       

        ! Total radiation absorbed by soil
        if (snow_depth_t .gt. 0.0) then 
            Rsoilabs = (Rsoilab1+Rsoilab2)*(1-albedo_snow)/(1-0.1)+Rsoilab3  
        elseif (water_table_depth .gt. 0.0) then 
            Rsoilabs = (Rsoilab1+Rsoilab2)*(1-albedo_water)/(1-0.1)+Rsoilab3  
        else
            Rsoilabs = Rsoilab1+Rsoilab2+Rsoilab3
        endif
        ! thermodynamic parameters for air
        rhocp  = cpair*Patm*AirMa/(Rconst*TairK)      
        H2OLv  = H2oLv0-2.365e3*Tair
        slope  = (esat(Tair+0.01)-esat(Tair))/0.01   
        psyc   = Patm*cpair*AirMa/(H2OLv*H2OMw)
        Cmolar = Patm/(Rconst*TairK)
        fw1    = AMIN1(AMAX1((FILDCP-wcl(1))/(FILDCP-WILTPT),0.3),1.0)      
        if (water_table_depth .gt. 0.0) then 
            Rsoil = 0. 
        else 
            Rsoil = 30.*exp(0.2/fw1)
        endif 
        rLAI      = exp(FLAIT)     
        ! latent heat flux into air from soil
        Esoil     = (slope*(Rsoilabs-G)+rhocp*Dair/(raero+rLAI))/       &
                    &   (slope+psyc*(Rsoil/(raero+rLAI)+1.))
        resht_lai = resht*FLAIT   
        Hsoil     = rhocp*(sftmp-Tair)/resht_lai  
    
        i = 1 ! zhou: this is for what?
        condu(i) = (FILDCP-wcl(i))*condu_air+liq_water(i)/(thksl(i)*0.01)*condu_water+ &
                    &   ice(i)/(thksl(i)*0.01)*condu_ice +(1-FILDCP)*condu_soil
        shcap(i) = (FILDCP-wcl(i))*shcap_air+liq_water(i)/(thksl(i)*0.01)*shcap_water+ &
                    &   ice(i)/(thksl(i)*0.01)*shcap_ice +(1-FILDCP)*shcap_soil
        difsv1   = 3600.*condu(i)/shcap(i)*10000.
            
        G        = condu(1)*(sftmp-tsoill(1))/(thksl(1)/2.*0.01)
        if (snow_depth_t .gt. 0.0) then 
            G    = condu_snow*(sftmp-Tsnow)/(snow_depth_t/2.*0.01)
        endif
        RESDH = Rsoilabs-Hsoil-Esoil-G    ! Residual heat energy.
        ! First derivative of net radiation; sensible heat; ground heat;
        DNR       = 4.*emsoil*sigma*(sftmp+273.2)**3
        DSH       = rhocp/resht_lai 
        DGH       = condu_s/(thksl(1)/2.*0.01)
        DLE       = (DNR+DGH)*slope/(slope+psyc*(Rsoil/(raero+rLAI)+1.))      
        drsdh     = -DNR-DSH-DGH-DLE
        DELTA     = resdh/drsdh ! Calculate increment DELTA.
        sftmp_pre = sftmp
        sftmp     = sftmp-DELTA
        if (ABS(sftmp_pre -sftmp) .gt. 20. ) sftmp=sftmp_pre
        tsoill_0  = sftmp
        ! zhou: for what?
        do i=1,10
            Tsoill_pre   = tsoill(i) 
            if (water_table_depth .lt. 0.0 .and. -water_table_depth .lt. depth_z(i)) then
                liq_water(i) = FILDCP*thksl(i)*0.01-ice(i)
            else
                liq_water(i) = wcl(i)*thksl(i)*0.01-ice(i)
            endif                       
            if (i .eq. 1) then 
                depth_z(1)=thksl(1)
            else 
                depth_z(i)=depth_z(i-1) + thksl(i)
            endif
            thkns2 = (thksl(i)+thksl(i+1))/2.   
            if (i .eq. 10) then
                difsv2=3600.*condu(i)/shcap(i)*10000. 
            else
                condu(i+1) = (FILDCP-wcl(i+1))*condu_air+liq_water(i+1)/(thksl(i+1)*0.01)*condu_water+ &
                                &   ice(i+1)/(thksl(i+1)*0.01)*condu_ice +(1-FILDCP)*condu_soil
                shcap(i+1) = (FILDCP-wcl(i+1))*shcap_air+liq_water(i+1)/(thksl(i+1)*0.01)*shcap_water+ &
                                &   ice(i+1)/(thksl(i+1)*0.01)*shcap_ice +(1-FILDCP)*shcap_soil 

                difsv2=3600.*condu(i+1)/shcap(i+1)*10000.
            endif
            
            temph2=(difsv1+difsv2)*(Tsoill(i)-Tsoill(i+1))/thkns2 

            if(i.eq.1) then ! adjust if there are snow or water layer above
                if (snow_depth_t .gt. 0.) then   
                    temph_snow = Amin1(diff_snow,difsv1)*(Tsnow-Tsoill(1))/((snow_depth_t+thksl(1))/2.)
                    Tsnow      = Tsnow+(exp(-depth_ex*snow_depth_t)*diff_snow*(sftmp-Tsnow)/(snow_depth_t/2.) &
                                    &   -temph_snow)/(snow_depth_t/2.+(snow_depth_t+thksl(1))/2.) 
                    Tsoill(1)  = Tsoill(1)+(temph_snow &
                                    &   -temph2)/((snow_depth_t+thksl(1))/2.+thkns2) 
                    if (Tsnow .gt.0.0) then 
                        Tsnow =0.0   
                        Tsoill(1)=0.
                    endif
                    drsdh     = 0.0    ! temporarily set drsdh =0 for heat adjustment of soil when  
                    tsoill_0  = (Tsoill(1)+Tsnow)/2.
                elseif (water_table_depth .gt. 0.) then  
                    temph_water = (3600.*condu_water/shcap_water*10000.+difsv1)*(Twater-Tsoill(1))/((water_table_depth+thksl(1))/2.)! there is snow layer 
                    Twater      = Twater+(2.*3600.*condu_water/shcap_water*10000.*(sftmp-Twater)/(water_table_depth/2.) &
                                    &   -temph_water)/(water_table_depth/2.+(water_table_depth+thksl(1))/2.) 
                    ! ----------- Phase change surface water 
                    if (Twater .lt. 0.0 .and. water_tw .gt. 0.0) then  ! freeze 
                        heat_excess = -(shcap_water/360000.*water_tw*100.-drsdh)*Twater
                        ice_incr    = heat_excess*3600./latent_heat_fusion/ice_density
                        if (ice_incr .lt. water_tw) then
                            ice_tw=ice_tw +ice_incr
                            water_tw = water_tw-ice_incr
                            Twater   = 0.0
                            Tice     = 0.0
                        else
                            ice_tw   = ice_tw +water_tw
                            water_tw = 0.0
                            Tice     = Tice - latent_heat_fusion*(ice_incr-water_tw)*ice_density/(shcap_ice*ice_tw)
                        endif     
                    elseif (Twater .gt. 0.0 .and. ice_tw .gt. 0.0) then    ! thraw              
                        heat_excess = (shcap_water/360000.*ice_tw*100.-drsdh)*Twater
                        ice_incr    = heat_excess*3600./latent_heat_fusion/ice_density                 
                        if (ice_incr .lt. ice_tw) then
                            ice_tw   = ice_tw -ice_incr
                            water_tw = water_tw+ice_incr
                            Twater   = 0.0
                            Tice     = 0.0
                        else
                            water_tw = water_tw +ice_tw
                            ice_tw   = 0.0
                            Twater   = Twater + latent_heat_fusion*(ice_incr-ice_tw)*ice_density/(shcap_water*water_tw)
                        endif
                    endif     ! ----------- end of phase change for surface layer --------------!  
                    temph2=(difsv1+3600.*condu_water/shcap_water*10000.)*(Tsoill(i)-Tsoill(i+1))/thkns2 
                    if (water_tw .eq. 0.0 .and. ice_tw .gt. 0.0) then 
                        Tsoill(1) = Tsoill(1)+(2.*3600.*condu_ice/shcap_ice*10000.*(Tice-Tsoill(1))/thkns1 &
                                    &   -temph2)/(thkns1+thkns2) 
                    else 
                        Tsoill(1) = Tsoill(1)+(2.*3600.*condu_water/shcap_water*10000.*(Twater-Tsoill(1))/thkns1 &
                                    &   -temph2)/(thkns1+thkns2) 
                    endif
                    drsdh = 0.0    ! temporarily set drsdh =0 for heat adjustment of soil       
                else                   
                    Tsoill(1) = Tsoill(1)+(diff_s*(sftmp-Tsoill(1))/thkns1 &
                                &   -temph2)/(thkns1+thkns2)
                endif
                !!!!!  phase change in top soil       
                heat_excess = drsdh*(thd_t-Tsoill(i))+shcap(i)*thksl(i)*(Tsoill(i)-thd_t)/360000.         
                ice_incr    = heat_excess*3600./latent_heat_fusion/ice_density               
                inter_var   = ice(i)   
                if (ice_incr .lt. 0.) then     ! freeze             
                    ice(i) = Amin1(liq_water(i)+inter_var,ice(i)-ice_incr)            
                else 
                    ice(i) = Amax1(ice(i)-ice_incr,0.0)              
                endif
                !! readjust energy and temp 
                heat_adjust = heat_excess-latent_heat_fusion*(inter_var-ice(i))*ice_density/3600.
                Tsoill(i)   = thd_t+heat_adjust/(shcap(i)*thksl(i)/360000.-drsdh)      
            else
                if ( i .gt. 9) then 
                    temph2 = 0.00003
                    thkns2 = 500  ! boundary conditions, rethink
                endif            
                Tsoill(i)   = Tsoill(i)+(temph1-temph2)/(thkns1+thkns2)    
                heat_excess = shcap(i)*thksl(i)*(Tsoill(i)-thd_t)/360000.        
                ice_incr    = heat_excess*3600./latent_heat_fusion/ice_density         
                inter_var   = ice(i) 
                if (ice_incr .lt. 0.) then     ! freeze             
                    ice(i)  = Amin1(liq_water(i)+inter_var,ice(i)-ice_incr)             
                else 
                    ice(i)  = Amax1(ice(i)-ice_incr,0.0)              
                endif         
                !! readjust energy and temp 
                heat_adjust = heat_excess-latent_heat_fusion*(inter_var-ice(i))*ice_density/3600.
                Tsoill(i)   = thd_t+heat_adjust/(shcap(i)/360000.*thksl(i))
            endif
            if (ABS(tsoill_pre -tsoill(i)) .gt. 5. ) Tsoill(i)=tsoill_pre
            TEMPH1 = TEMPH2
            THKNS1 = THKNS2
            DIFSV1 = DIFSV2        
        enddo
        testout(1)=tsoill_0
        testout(2:11)=tsoill(1:10)
        deallocate(depth_z)
    return 
    end subroutine soil_thermal

    ! --------------------------------------------
    subroutine soil_water(transp, evap, liq_water, ice, accumTa, &
                        & wsc, phi)
        ! All of inputs, the unit of water is 'mm', soil moisture or soil water content is a ratio
        implicit none
        real transp, evap   ! from canopy?
        real liq_water(10), ice(10)  ! from soil them?
        real melt           ! from driver module of snow_d
        real runoff         ! calculate in this part, and be used in update part and transferCpool module
        ! real infilt, fwsoil, topfws, omega, zwt  ! has the initial data, and update in this part, maybe used to calculate the water processes
        real accumTa        ! from driver
        ! outputs ----------------------!
        real wsc(10), phi   ! phi used in the methane module
       
        !--------------local vars ------- !
        integer i,j,k
        real    infilt_max
        real(KIND=8) FLDCAP, WILTPT ! soil traits: ie. 0.xx
        real DWCL(10), evapl(10), WUPL(10),SRDT(10),depth(10)
        integer nfr                              ! !   plant traits
        real rain_new,rain_t,snowDepth  ! int added from soil thermal
        real infilt_dbmemo, twtadd, wtadd
        real exchangeL,supply,demand,omegaL(10)
        real Tsrdt
        real vtot     ! characters annotation for water table module  -MS
        real zmax,thetasmin,zthetasmin,az
        real zwt1,zwt2,zwt3
        real fw(10),ome(10)
        ! ------- output? -----------------------
        real  tr_allo, tr_ratio(10), plantup(10)

        infilt_max = 15.
        WILTPT     = wsmin/100.000
        FLDCAP     = wsmax/100.000
        do i=1,10
            dwcl(i)  = 0.0
            evapl(i) = 0.0
            WUPL(i)  = 0.0
            SRDT(i)  = 0.0
            DEPTH(i) = 0.0
        enddo
        ! Determine which layers are reached by the root system. 
        ! Layer volume (cm3)
        DEPTH(1) = THKSL(1)
        DO i=2,10
            DEPTH(i) = DEPTH(i-1) + THKSL(i)
        enddo
        do i=1,10
            IF(rdepth.GT.DEPTH(i)) nfr=i+1
        enddo
        IF (nfr.GT.10) nfr=10
        !   *** ..int    
        !   ******** added for soil thermal    

        if (doSoilphy) then 
            rain_new = rain
            if (accumTa .lt. -4.) rain_new =0.           !dbice   !tuneice
            ! here it defines how the melt water is added to water input   
            ! add melted water hourly
            rain_t = melt/24+rain_new
            ! add melted water daily all at once
            infilt = infilt+rain_t
            if (ice(1) .gt. 0.0) then
                !infilt = 0.0
            endif
        else
            infilt = infilt+rain
        endif     
        infilt_dbmemo = infilt
        ! Loop over all soil layers.
        TWTADD = 0
        IF(infilt.GE.0.0)THEN
            ! Add water to this layer, pass extra water to the next.
            ! cc = wcl(1)                            !dbmemo
            WTADD  = AMIN1(INFILT,infilt_max,AMAX1((FLDCAP-wcl(1))*thksl(1)*10.0,0.0)) ! from cm to mm
            ! change water content of this layer
            ! write(*,*) 'before  update',wcl(1)    !dbmemo
            WCL(1) = (WCL(1)*(thksl(1)*10.0)+WTADD)/(thksl(1)*10.0)
            ! ddd = (FLDCAP-cc)*thksl(1)*10.0        !dbmemo
            ! dd = (FLDCAP-0.564999998)*thksl(1)*10.0   !dbmemo
            ! write (*,*) 'wcl(1)',wcl(1), 'WTADD',WTADD,'INFILT',INFILT,'infilt_max',infilt_max,'ddd',ddd    !dbmemo
            ! FWCLN(I)=WCL(I)       !  /VOLUM(I)! update fwcln of this layer
            TWTADD = TWTADD+WTADD   ! calculating total added water to soil layers (mm)
            INFILT = INFILT-WTADD   ! update infilt
        ENDIF
        ! runoff method 1 !dbice
        if (doSoilphy) then 
            runoff = INFILT*0.005
        else
            runoff = INFILT*0.001              ! Shuang added this elseif line
        endif   
        infilt = infilt-runoff
        ! -----------------------------------------------------------------------
        if (transp .gt. 0.2 .and. transp .le. 0.22) then
            infilt = infilt+transp*0.4
        else if (transp .gt. 0.22) then
            ! infilt = infilt+infilt*0.0165
            ! infilt = infilt+0.22*0.4+(transp-0.22)*0.9
            infilt = infilt+transp*0.8
        else
            infilt = infilt+transp*0.001
        endif
        ! !    
        if (evap .ge. 0.1 .and. evap .le. 0.15) then
            infilt = infilt+evap*0.4
        else if (evap .gt. 0.15) then
            infilt = infilt+evap*0.8
        else
            infilt = infilt+evap*0.001
        endif
        ! --------------------------------------------------------------------------- 
        do i=1,10 !   water redistribution among soil layers
            wsc(i) = Amax1(0.00,(wcl(i)-wiltpt)*THKSL(i)*10.0)
            ! ..int commented lines for soil thermal        
            ! omegaL(i)=Amax1(0.001,(wcl(i)-WILTPT)/(FLDCAP-WILTPT))
            if (doSoilphy) then 
                omegaL(i) = Amax1(0.001,(liq_water(i)*100./thksl(i)-WILTPT)/(FLDCAP-WILTPT))
            else
                omegaL(i) = Amax1(0.001,(wcl(i)-WILTPT)/(FLDCAP-WILTPT))
            endif        
        enddo
        supply = 0.0
        demand = 0.0

        do i=1,9
            if(omegaL(i).gt.0.3)then
                supply    = wsc(i)*(omegaL(i)-0.3)
                demand    = (FLDCAP-wcl(i+1))*THKSL(i+1)*10.0      &
                            &  *(1.0-omegaL(i+1))
                exchangeL = AMIN1(supply,demand)
                wsc(i)    = wsc(i)- exchangeL
                wsc(i+1)  = wsc(i+1)+ exchangeL
                wcl(i)    = wsc(i)/(THKSL(i)*10.0)+wiltpt
                wcl(i+1)  = wsc(i+1)/(THKSL(i+1)*10.0)+wiltpt
            endif
        enddo
        wsc(10) = wsc(10)-wsc(10)*0.00001     ! Shuang modifed
        runoff  = runoff+wsc(10)*0.00001     ! Shuang modifed
        wcl(10) = wsc(10)/(THKSL(10)*10.0)+wiltpt
        Tsrdt   = 0.0
        DO i = 1,10  ! Fraction of SEVAP supplied by each soil layer
            SRDT(I) = EXP(-6.73*(DEPTH(I)-THKSL(I)/2.0)/100.0) !/1.987
            ! SRDT(I)=AMAX1(0.0,SRDT(I)*(wcl(i)-wiltpt)) !*THKSL(I))
            Tsrdt   = Tsrdt+SRDT(i)  ! to normalize SRDT(i)
        enddo
        do i=1,10
            EVAPL(I) = Amax1(AMIN1(evap*SRDT(i)/Tsrdt,wsc(i)),0.0)  !mm
            DWCL(I)  = EVAPL(I)/(THKSL(I)*10.0) !ratio
            wcl(i)   = wcl(i)-DWCL(i)
        enddo
        evap = 0.0       
        do i=1,10
            evap = evap+EVAPL(I)
        enddo
        ! Redistribute transpiration according to root biomass
        ! and available water in each layer
        tr_allo = 0.0
        do i=1,nfr
            tr_ratio(i) = FRLEN(i)*wsc(i) !*(wcl(i)-wiltpt)) !*THKSL(I))
            tr_allo     = tr_allo+tr_ratio(i)
        enddo
        do i=1,nfr
            plantup(i) = AMIN1(transp*tr_ratio(i)/tr_allo, wsc(i)) !mm              
            wupl(i)    = plantup(i)/(thksl(i)*10.0)
            wcl(i)     = wcl(i)-wupl(i)
        enddo
        transp = 0.0
        do i=1,nfr
            transp = transp + plantup(i)
        enddo
        !******************************************************    
        !   water table module starts here
        if (doSoilphy) then
            vtot = (liq_water(1)+liq_water(2)+liq_water(3))*1000+(ice(1)+ice(2)+ice(3))*1000+infilt
        else 
            vtot = wsc(1)+wsc(2)+wsc(3)+infilt
        endif

        !   infilt means standing water according to jiangjiang
        !    vtot = MAX(145.,vtot+145.+rain-evap-transp-runoff)         ! vtot should not be smaller than 145, which is the water content when wt is at -300mm
        phi        = 0.56                           ! soil porosity   mm3/mm3   the same unit with theta
        zmax       = 300                            ! maximum water table depth   mm
        thetasmin  = 0.25                           ! minimum volumetric water content at the soil surface   cm3/cm3
        zthetasmin = 100                            ! maximum depth where evaporation influences soil moisture   mm
        az         = (phi-thetasmin)/zthetasmin     ! gradient in soil moisture resulting from evaporation at the soil surface    mm-1

        zwt1       = -sqrt(3.0*(phi*zmax-vtot)/(2.0*az))
        zwt2       = -(3.0*(phi*zmax-vtot)/(2.0*(phi-thetasmin)))
        zwt3       = vtot-phi*zmax                                   
        if ((zwt1 .ge. -100) .and. (zwt1 .le. 0))   zwt = zwt1  !the non-linear part of the water table changing line
        if (zwt2 .lt. -100)                         zwt = zwt2  !the linear part of the water table changing line
        if (phi*zmax .lt. vtot)                     zwt = zwt3  !the linear part when the water table is above the soil surface   
        !   ..int new lines added for soil thermal module 
        do i=1,nfr       
            if (doSoilphy) then 
                ome(i) = (liq_water(i)*100./thksl(i)-WILTPT)/(FLDCAP-WILTPT)
            else 
                ome(i) = (wcl(i)-WILTPT)/(FLDCAP-WILTPT)
                ome(i) = AMIN1(1.0,AMAX1(0.0,ome(i)))
            endif 
            fw(i)=amin1(1.0,3.333*ome(i))
        enddo

        if (doSoilphy) then 
            topfws=amax1(0.0,topfws)
        else 
            topfws=amin1(1.0,(wcl(1)-WILTPT)/((FLDCAP-WILTPT)))
        endif     

        fwsoil = 0.0
        omega  = 0.0
        do i=1,nfr
            fwsoil = fwsoil+fw(i)*frlen(i)
            omega  = omega+ome(i)*frlen(i)
        enddo   
        write(*,*)"test-omega:", omega, ome
        return
    end subroutine soil_water

    ! ----------------------------------------------------------
    subroutine snow_daily(rain_d,lat,days,ta,snow_dsim,fa,fsub,rho_snow,melt,dcount,decay_m)
                        ! rain_d,lat,iday,accumTa,snow_dsim,fa,fsub,rho_snow,melt,dcount,decay_m
        ! subroutines from methane and soil thermal modules
        implicit none
        real lat,tr,daylength,dec,melt,fa,sublim,dsnow,snow_in,decay_m,fsub
        real rain_d,snow_dsim,rho_snow,dcount,ta
        integer days
        real snow_dsim_pre
        tr        = 0.0174532925
        dec       = sin(((real(days)-70.)/365.)*360.*tr)*23.44
        daylength = acos(-tan(lat*tr)*tan(dec*tr))/7.5 
        daylength = daylength/tr/24.
        if (snow_dsim .ge. 0.) then
            dcount = dcount +1.
        else 
            dcount = 0.
        endif
        sublim = 0.
        if (ta .gt. 0. .and. snow_dsim .gt. 0.) sublim = fsub*715.5*daylength*esat(ta)/(ta+273.2)*0.001   ! 0.001 from Pa to kPa
        melt   = 0.
        if (ta .gt. 1.0e-10 .and. snow_dsim .gt. 0.) melt=fa*(2.63+2.55*ta+0.0912*ta*rain_d)   !dbmemo updated version
        if (dcount .gt.0. .and. ta .lt.5.) then
            melt = melt*EXP(-decay_m*dcount/365.)  !dbmemo dbice
        endif
        if (ta .le. 0.) then         ! dbmemo second bug in dbmemo
            snow_in = rain_d
        else
            snow_in = 0.
        endif
        dsnow         = snow_in-sublim-melt 
        snow_dsim_pre = snow_dsim
        snow_dsim     = snow_dsim + dsnow/rho_snow 
        if (snow_dsim .le. 0.0) then 
            snow_dsim = 0.0 
            melt      = snow_dsim_pre*rho_snow +snow_in-sublim    !! for water part
        endif 
        melt = AMAX1(melt, 0.)
    return
    end subroutine snow_daily
    
    ! ----------------------------------------------------------
    real function esat(T)
        real T
        ! returns saturation vapour pressure in Pa
        esat = 610.78*exp(17.27*T/(T+237.3))
    return
    end
end module mod_soil