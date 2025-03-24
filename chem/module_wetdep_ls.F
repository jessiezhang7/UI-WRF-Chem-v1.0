MODULE module_wetdep_ls
USE module_state_description, only:p_qv,p_qc,p_so2,p_sulf, &
                              p_bc1,p_bc2,p_oc1,p_oc2,p_seas_1,p_seas_2,     &
                              p_seas_3,p_seas_4,p_dms,                       &
                              !HZ, adding some new for MADE/SORGAM   
						      p_so4ai,p_so4aj,p_no3ai,p_no3aj,p_nh4ai,       &
						      p_nh4aj,p_p25i,p_p25j,p_dust_1,p_dust_2,         &
						      p_dust_3,p_dust_4,p_dust_5,                      &
						      p_eci,p_ecj,p_orgaro1j,p_orgaro1i,p_orgaro2j,p_orgaro2i,&
                              p_orgalk1j,p_orgalk1i,p_orgole1j,p_orgole1i,p_orgba1j,&
                              p_orgba1i,p_orgba2j,p_orgba2i,p_orgba3j,p_orgba3i,p_orgba4j,&
                              p_orgba4i,p_orgpaj,p_orgpai
						      
CONTAINS
subroutine wetdep_ls(dt,var,rain,moist,rho,num_moist, &
         num_chem,numgas,dz8w,vvel,chem_opt,             &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
  IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: num_chem,numgas,num_moist,               &
                                  chem_opt,                                &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   real, INTENT(IN ) :: dt
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(IN ) ::                                   moist
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) :: rho,dz8w,vvel        
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ,1:num_chem),                        &
          INTENT(INOUT) :: var        
   REAL,  DIMENSION( ims:ime , jms:jme ),                                  &
          INTENT(IN   ) :: rain
!  REAL,  DIMENSION( ims:ime ,  jms:jme ),                                  &
!         INTENT(INOUT   ) :: var_rmv
   REAL,  DIMENSION( its:ite ,  jts:jte ) :: var_sum,var_rmv
   REAL,  DIMENSION( its:ite ,  kts:kte, jts:jte ) :: var_rmvl
   REAL,  DIMENSION( its:ite ,  jts:jte ) :: frc,var_sum_clw,rain_clw     
    real :: dvar,factor,clsum,alpha,rho_water
   integer :: nv,i,j,k,km,kb,kbeg
    rho_water = 1000.
!   write(0,*) 'in wetdepls, numgas,num_chem = ',numgas,num_chem,chem_opt
!   nv=p_bc1
    do nv=1,num_chem
       if(nv.le. numgas .and. nv.ne.p_sulf)cycle
       alpha = .5    ! scavenging factor
       if(chem_opt >= 300)then
          if(nv.eq.p_bc1 .or. nv.eq.p_oc1 .or. nv.eq.p_dms)cycle
          if(nv.eq.p_bc2 .or. nv.eq.p_oc2)alpha=0.8
       endif
       !HZ, the orginal, needs to be commented out
       !if(nv.eq.p_sulf .or. nv.eq.p_seas_1 .or. nv.eq.p_seas_2    &
        !               .or. nv.eq.p_seas_3 .or. nv.eq.p_seas_4)alpha=1.
       
       !HZ-Jessie, 12/10/2024, reset the scavenging efficiency, decrease from convective scavenging
!		if (nv.eq.p_so4ai .or. nv.eq.p_nh4ai .or. nv.eq.p_no3ai .or. nv.eq.p_p25i) alpha = 0.
!		if (nv.eq.p_so4aj .or. nv.eq.p_nh4aj .or. nv.eq.p_no3aj .or. nv.eq.p_p25j) alpha = 0.5
!		if (nv.eq.p_dust_1 .or. nv.eq.p_dust_2 ) alpha=0.1
!		if (nv.eq.p_dust_3 ) alpha=0.2
!		if (nv.eq.p_dust_4 ) alpha=0.6
!		if (nv.eq.p_dust_5 ) alpha=0.8

		!HZ-Jessie, 12/17/2025, readjusting the efficiency
		if (nv.eq.p_so4ai .or. nv.eq.p_nh4ai .or. nv.eq.p_no3ai .or. nv.eq.p_p25i)  alpha = 0.
		if (nv.eq.p_p25j) alpha = 0.4
		if (nv.eq.p_so4aj) alpha = 0.8
		if (nv.eq.p_no3aj) alpha = 0.6
		if (nv.eq.p_nh4aj) alpha = 0.8
		
		if (nv.eq.p_dust_1 .or. nv.eq.p_dust_2 )  alpha=0.1
		if (nv.eq.p_dust_3 )  alpha=0.2
		if (nv.eq.p_dust_4 )  alpha=0.6
		if (nv.eq.p_dust_5 )  alpha=0.8
		
		!HZ-Jessie, 12/13/2024, need adjust EC,OC and seasalt as well 
		if (nv.eq.p_eci ) alpha = 0.
		if (nv.eq.p_ecj ) alpha = 0.4
		
		if (nv.eq.nv.eq.p_orgole1i .or. nv.eq.p_orgpai .or. nv.eq.p_orgba1i) alpha = 0.
		if (nv.eq.p_orgaro1i .or. nv.eq.p_orgaro2i .or. nv.eq.p_orgalk1i ) alpha = 0.
		if (nv.eq.p_orgba2i .or. nv.eq.p_orgba3i .or. nv.eq.p_orgba4i) alpha = 0.
		
		if (nv.eq.nv.eq.p_orgole1j .or. nv.eq.p_orgpaj .or. nv.eq.p_orgba1j) alpha = 0.4
		if (nv.eq.p_orgaro1j .or. nv.eq.p_orgaro2j .or. nv.eq.p_orgalk1j) alpha = 0.4
		if (nv.eq.p_orgba2j .or. nv.eq.p_orgba3j .or. nv.eq.p_orgba4j) alpha = 0.4
	
		if (nv.eq.p_seas_1 ) alpha=0.4
		if (nv.eq.p_seas_2 ) alpha=0.4
		if (nv.eq.p_seas_3 ) alpha=0.8
		if (nv.eq.p_seas_4 ) alpha=1.0
		
       do i=its,ite
       do j=jts,jte
        var_sum_clw(i,j)=0.
        var_sum(i,j)=0.
        var_rmvl(i,:,j)=0.
        frc(i,j)=0.
        rain_clw(i,j)=0.
        if(rain(i,j).gt.1.e-10)then
!          write(0,*)i,j,rain(i,j)
! convert rain back to rate
!
           rain_clw(i,j)=rain(i,j)/dt
! total cloud water
!
           do k=1,kte-1
              dvar=max(0.,moist(i,k,j,p_qc)*rho(i,k,j)*vvel(i,k,j)*dz8w(i,k,j))
              var_sum_clw(i,j)=var_sum_clw(i,j)+dvar
              var_sum(i,j)=var_sum(i,j)+var(i,k,j,nv)*rho(i,k,j)
           enddo
           if(var_sum(i,j).gt.1.e-10 .and. var_sum_clw(i,j).gt.1.e-10 ) then
   !        assuming that frc is onstant, it is my conversion factor 
!       (just like in convec. parameterization
              frc(i,j)=rain_clw(i,j)/var_sum_clw(i,j)
!    write(0,*)'frc ', frc(i,j),var_sum_clw(i,j),var_sum(i,j)
              frc(i,j)=max(1.e-6,min(frc(i,j),.005))
           endif
        endif
      enddo
    enddo
!
! get rid of it
!
    do i=its,ite
    do j=jts,jte
     if(rain(i,j).gt.1.e-10 .and. var_sum(i,j).gt.1.e-10 .and. var_sum_clw(i,j).gt.1.e-10)then
       do k=kts,kte-2
        if(var(i,k,j,nv).gt.1.e-16 .and. moist(i,k,j,p_qc).gt.0.)then
        factor = max(0.,frc(i,j)*rho(i,k,j)*dz8w(i,k,j)*vvel(i,k,j))
!       write(0,*)'var before ',k,km,var(i,k,j,nv),factor
        dvar=alpha*factor/(1+factor)*var(i,k,j,nv)
        var(i,k,j,nv)=max(1.e-16,var(i,k,j,nv)-dvar)
!       write(0,*)'var after ',km,var(i,k,j,nv),dvar
        endif
       enddo
!      var_rmv(i,j)=var_rmv(i,j)+var_rmvl(i,j)
    endif
    enddo
    enddo
    enddo
END SUBROUTINE WETDEP_LS

END MODULE module_wetdep_ls


