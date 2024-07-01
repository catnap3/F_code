! initialize =========================================================
! initialize =========================================================
      subroutine initial
        use globals ! import "globals" module
		use load_wna_lat_mod
		implicit none 
		include 'mpif.h'

! Set parameters --------------------------------------------------------
        aqm     = abs(qm)
        cs      = cv*cv;	tcs = 2.0d0*cs;               ! set (light speed)**2
		mu0     = 1;                                      ! Vacuum permeability
		eps0    = 1/(mu0*cs);                             ! Vacuum permittivity
		tpi     = 2.0d0*pi;                               ! set parameter 2*pi		
		zunit   = real_cv/real_wce;                       ! set real number of unit length [m]
		REr     = RE/zunit*cv/wceq;                       ! Earth radii unit dh
		radeq   = Lv*RE*real_wce/real_cv*cv/wceq             ! distance at the equator of the L-value [dh] 
		w_ph    = wpeq*dsqrt(N_ph); w_phs = w_ph**2       ! plasma frequency of source electron	
		w_phc   = wpeq*dsqrt(N_phc); w_phcs = w_phc**2       ! plasma frequency of source electron	

		! Nomalized -----------------------------------------------------		
        vpro   = vpro*cv;                                 ! average perpendicular velocity	
		Vtpara = Vtpara*cv; 
		
! Setting ----------------------------------------------------------------
		if (test==1) then
		  jobng  = 0
		  jobnp  = 0  
		  jobne  = 0
		  joblo  = 1 
		  jobnk  = 0
		  nphase = 3600
		end if

! Observed ---------------------------------------------------------------
		!trajectory @test==1
		tra_time   = kidint(lifetime/real_dt) !observed time period [sec]
		tra_start1 = 0.0d0
		tra_start  = kidint(tra_start1/real_dt)
		
		!calculation time period [sec]
		if(test==1)then
		   ntime = tra_time+tra_start	   
		elseif(waveonly .ne. 0)then
		   ntime = kidint(3/real_dt)	
		else
		   ntime = kidint(lifetime/real_dt)		   
		end if

		pskip_trajectory=3000
		if(waveonly==1)then
		   pskip_wave = kidint(0.0004/real_dt)
		else 
		   pskip_wave = kidint(0.0001/real_dt) !waveonly=2
		end if

		pskip_h      = 4
		p_trajectory = int8(tra_time/pskip_trajectory)+1
		p_wave       = int8(ntime/pskip_wave)+1	

! Set numbers associate with Green's function -----------------------------
        if(keido_del==0.0d0)then                               ! number if longitudinal grid
		   nkeido=1     
		else
		   nkeido=int((keido_max-keido_min)/keido_del)
		   if(nkeido==0)then
			  nkeido=1
		   end if
		end if

        nenergy = kidint((energy_max-energy_min)/energy_del)+1    ! number of kinetic energy grid	
		nenergy1= nenergy
		npitch  = kidint((pitcheq_max-pitcheq_min)/pitcheq_del)+1 ! number of equatorial pitch angle grid
		if(test==1 .or. waveonly==1  .or. waveonly==2)then
		    keido_min   = 0.0d0 
		    keido_del   = 0.0d0
		    nkeido      = 1
			nenergy     = 1
			npitch      = 1
        end if
		totalp  = nenergy*npitch*nphase*nkeido	! number of total particles  
		
! Particles in process---------------------------------------------------
		totalp=nenergy*npitch*nphase*nkeido ! total particle
		if(mod(totalp, nprocs).ne.0)then
		   np = int8(totalp/nprocs)+1  !number of particles per process
		else
		   np = int8(totalp/nprocs)
		end if
		istart= myrank*np+1
		iend  = (myrank+1)*np
		if (myrank .eq. 0) then
			print *,'myrank=',myrank,'totalp=',totalp
			print *,'np',np,'nprocs',nprocs
		end if		
		
! Allocate for particles--------------------------------------------------------------
		allocate( px(np), py(np), pz(np), hh(np), vx(np), vy(np), vz(np), &
			vpara(np), vperp(np), phi(np), KE(np), PA(np), check(np))	
allocate( trajectoryrx(np), trajectoryry(np), trajectoryrz(np), Bmoment(np)	)		
		check(:)=2 !  0:move 1:precipitation 2:non-memory		

! Allocate for graph------------------------------------------------------
		
		allocate(wfBp(p_wave,1),BwBp(p_wave,1),opt(p_wave), thre(p_wave), &
					  Keq_first(nenergy1,nx),Keq(nenergy1,nx))
		
			  if(test==1)then
				 allocate(trajectory(np),trajectorypitch(np), &
						  trajectorylambda(np), & 
						  trajectoryz(np),trajectoryst(np),trajectorye(np))
				 allocate(tp_Green(np))
			  else
				 allocate(tp_Green(np))
			  end if
				
	
! First --------------------------------------------------------------
		if(waveonly==1 .or. waveonly==2)then
	      if(waveonly==1)then
              Keq_first(:,:)=0 !S parameter check
              Keq(:,:)=-2	
		  else
              opt(:)=-1
              thre(:)=-1
          end if			  
          BwBp(:,:)=0
          wfBp(:,:)=0  
		  tp_wave=0
		else 
		  if(test==1)then
			trajectory=-1.0d0
	        trajectorypitch=-1.0d0
			trajectorylambda=0.0d0
	        trajectoryz=1000.0d0
	        trajectoryst=1.0d0
			trajectorye=-1.0d0
			tp=0
	      else
		    tp_Green=0
		  end if
        end if
		
	      if(waveonly==1)then
              Keq_first(:,:)=0 !S parameter check
              Keq(:,:)=-2	
		  else
              opt(:)=-1
              thre(:)=-1
          end if			  
          BwBp(:,:)=0
          wfBp(:,:)=0  
		  tp_wave=0
	 
		  if(test==1)then
			trajectory=-1.0d0
	        trajectorypitch=-1.0d0
			trajectorylambda=0.0d0
	        trajectoryz=1000.0d0
	        trajectoryst=1.0d0
			trajectorye=-1.0d0
			tp=0
	      else
		    tp_Green=0
		  end if
      
      
		aaa=0 !whistler wave packet duration at the equator
		check(:)=2 !  0:move 1:precipitation 2:non-memory	  
		swave=0	  		
		
		
!#################################################################################
		rhw = int(dh/dhw) ! ratio of dh/dhw
		index_w(1) = 1
		do i=2,nhh1
			index_w(i) = index_w(i-1)+rhw  
		end do
		
! Grid of the dipole--------------------------------------------------------------
        ! for grids
		Beq        = abs(wceq/qm)
		radeq3     = radeq*radeq*radeq  !R^3
		threei     = 1.0d0/3.0d0
		cet        = nhh1 
		xbh        = radeq		
		zbh        = 0.5d0         !half move
		xb         = radeq         !R_EQ
		zb         = 0.0d0
		bbb(1)     = Beq     
		xg(cet)    = xb
		zg(cet)    = zb
		hg(cet)    = 0.0d0
		lambdag(1) = 0.0d0
		Bcheck     = Bcheck*Beq;
		Bcut       = Bcut*Beq; 
		iconvective= 1
		iv2=0
		do i = 1,nhh
		   r2       = xbh*xbh + zbh*zbh
		   wratio   = -3.0d0*Beq*radeq3/(r2*r2*dsqrt(r2))
		   bx0      = wratio*xbh*zbh
		   bz0      = wratio*(zbh*zbh - r2*threei)
		   bab      = dsqrt(bx0*bx0 + bz0*bz0)
		   xb       = xb + bx0/bab                    !bx0/bab*1 (deltah=1)
		   zb       = zb + bz0/bab
		   r2       = xb*xb + zb*zb
		   wratio   = -3.0d0*Beq*radeq3/(r2*r2*dsqrt(r2))
		   bx0      = wratio*xb*zb
		   bz0      = wratio*(zb*zb - r2*threei)
		   bab      = dsqrt(bx0*bx0 + bz0*bz0)
		   bbb(1+i) = bab   ! total magnetic field
		   xg(cet+i) = xb
		   xg(cet-i) = xb
		   zg(cet+i) = zb
		   zg(cet-i) = -zb
		   hg(cet+i) = dble(i)
		   hg(cet-i) = dble(-i)
		   lambdag(1+i)= datand(zg(cet+i)/xg(cet+i))  ! Latitude in degrees
		   if(lambdag(1+i)<2)then
				iconvective=i+1
		   end if
		   xbh = xb + bx0/bab*0.5d0
		   zbh = zb + bz0/bab*0.5d0
		end do
		
		iconvective = (iconvective-1)*rhw+1
		nhh2   = nhh
		polar_latitude=datand(zg(nhp1)/xg(nhp1))
		
		! for wave equations		
		xbh        = radeq		
		zbh        = 0.5d0         !half move
		xb         = radeq         !R_EQ
		zb         = 0.0d0
		w_bbb(1)   = Beq
		w_wc(1)    = aqm*Beq
		w_xg(1)    = xb
		w_zg(1)    = zb
		w_hg(1)    = 0.0d0
		do i = 1,nhhw
		   r2       = xbh*xbh + zbh*zbh
		   wratio   = -3.0d0*Beq*radeq3/(r2*r2*dsqrt(r2))
		   bx0      = wratio*xbh*zbh
		   bz0      = wratio*(zbh*zbh - r2*threei)
		   bab      = dsqrt(bx0*bx0 + bz0*bz0)
		   xb       = xb + bx0/bab*dhww                    !bx0/bab*1 (deltah=1)
		   zb       = zb + bz0/bab*dhww
		   r2       = xb*xb + zb*zb
		   wratio   = -3.0d0*Beq*radeq3/(r2*r2*dsqrt(r2))
		   bx0      = wratio*xb*zb
		   bz0      = wratio*(zb*zb - r2*threei)
		   bab      = dsqrt(bx0*bx0 + bz0*bz0)
		   w_bbb(1+i) = bab   ! total magnetic field
		   w_wc(1+i)  = aqm*bab
		   w_xg(1+i) = xb
		   w_zg(1+i) = zb
		   w_hg(1+i) = dble(i)
		   w_lambdag(1+i)= datand(w_zg(1+i)/w_xg(1+i))  ! Latitude in degrees
		   
		   xbh = xb + bx0/bab*0.5d0*dhww
		   zbh = zb + bz0/bab*0.5d0*dhww
		end do	
!##########################################################################################		
! Sbb(R,alpha_eq) given by (10) of Omura et al. 2015----------------------
		Sbb = 0
		do bi = 1, nx-1
			a1 = dsqrt(1.0d0-1.0d0*(dsind(dble(bi-1)*dx))**2) !Equator
            Sbb(bi) = Sbb(bi)+1.0/a1
            do bj = 2,nhh
				a1=dsqrt(1.0-bbb(bj)/beq*(dsind(dble(bi-1)*dx))**2)
				if(isnan(a1))then
				else
					Sbb(bi) = Sbb(bi)+2.0d0*1.0/a1
				end if
			end do
		end do		

!#################################################################################
! background gyrofrquency & plasma frequency**************************************
		wc(1)  = wceq
		do bi = 2, nhh1
			wc(bi) = aqm*bbb(bi)
		end do
		wpe2= wpeq*wpeq                   ! (plasma frequency)**2
!#################################################################################
! Oblique whistler mode wave******************************************************  
!*********************************************************************************
		! variables
		Bw(:)  = 0;
		Bw(1)  = Beq*1;	
		w_Bw(:)= 0;
		w_Bw(1)= Beq*1;
		ww_min = ww_min*wceq
		ww_max = ww_max*wceq
		wf     = ww_min		
		w_wf   = ww_min			
		AAAeq  =2.0d0/pi*vpro**2/Vtpara**2-1 ! assume vpro=sqrt(pi/2)*Vtperp	
        
		
		
		
		! Set wave normal angles		
	    fpl   = 90-asind(1/dsqrt(Lv)); ! Latitute of the footpoint of a field line   	
		Lat   = (/ (ii, ii = 0, nhh) /); 
		Lat   = Lat*(fpl+0.1)/nhh1;          ! Set Latitude array		
		idx_L1            = maxloc(Lat , Lat <= 2);
		thetaL(1:idx_L1(1))  = 0;                    ! Lat <= 2: wave normal angle is 0
		idx_L2            = minloc(Lat , Lat >= 45);
		thetaL(idx_L2(1):nhh1)= wna;                  ! Lat >=45: wave normal angle is theta(set in input_para) [unit:deg]
		do ai = (idx_L1(1)+1),(idx_L2(1)-1)               ! 2 < Lat < 45 : linear
			thetaL(ai) = wna/(idx_L2(1)-idx_L1(1))*(ai-idx_L1(1)); 
		end do			
		! !====for wna close to resonance cone  angle==========================
		! fpl      = 90-dasind(1/dsqrt(Lv)); ! Latitute of the footpoint of a field line   		
		! Lat      = (/ (ii, ii = 0, nhh) /); 
		! Lat      = Lat*(fpl+0.1)/nhh1;          ! Set Latitude array	
		! allocate(XX(nhh1), YY(nhh1), RR(nhh1), LL(nhh1), PP(nhh1), SS(nhh1))	
		! XX       = wpe2/(ww_max*ww_max);
	    ! YY       = wc/ww_max;
		! PP       = 1-XX;
		! RR       = 1 + XX/(YY-1);
		! LL       = 1 - XX/(YY+1);
		! SS       = 0.5d0*(RR+LL);		
		! thetaL = dasind(dsqrt(-PP/(SS-PP)))*0.9d0 ! wna = resonance cone angle 90%		
		! idx_L1            = maxloc(Lat , Lat <= 2);
		! thetaL(1:idx_L1(1))  = 0;                    ! Lat <= 2: wave normal angle is 0
		! idx_L2            = minloc(Lat , Lat >= 10); 	
		! do ai = (idx_L1(1)+1),(idx_L2(1)-1)               ! 2 < Lat < 10 : linear
			! thetaL(ai) = thetaL(idx_L2(1))/(idx_L2(1)-idx_L1(1))*(ai-idx_L1(1)); 
		! end do	
		
		! deallocate(XX, YY, RR, LL, PP, SS)
		! !======================================================
	
		!$OMP PARALLEL DO PRIVATE(bi,a1,dhL,hjj,hii,ii1,inth)
		do bi = 1, nhh1       !! linear interpolation for thetag 
			a1        = lambdag(bi);
			if (a1 <= 2) then
				thetag(bi) = 0
			else if (a1 >= 45) then
				thetag(bi) = wna
			else	
				dhL        = a1-Lat;
				hjj        = minloc(dhL, dhL .ge. 0);
				hii        = hjj(1);
				ii1        = hii+1;
				inth       = abs(Lat(hii)-Lat(ii1));  		
				thetag(bi) = (thetaL(ii1)*abs(a1-Lat(hii))+thetaL(hii)*abs(a1-Lat(ii1)))/inth;  				
			end if	
        end do
		!$OMP END PARALLEL DO

		
		
		!$OMP PARALLEL DO PRIVATE(bi,hii,ii1,ddh)
		do bi = 1, nhhw1
			hii = int((bi-1)*dhww)+1
			ii1 = hii+1
			ddh = (bi-1)*dhww+1-hii
			w_thetag(bi) = thetag(hii)*(1.0d0-ddh) + thetag(ii1)*ddh
		end do	
		!$OMP END PARALLEL DO		
				
		Bp = cshift(bbb, 1)         ! calculate dB/dh
		Bn = cshift(bbb, -1)
		dBz = 0.5*(Bp-Bn)
		dBz(1)   = 2*dBz(2)-dBz(3)	
		dBz(nhh1) = 2*dBz(nhh1-1)-dBz(nhh1-2)
		
! make table for rk and Vg interpolation
        hgrid    = (/ (ii, ii = 1, nhhw1) /)
		index_a  = (/ (ii, ii = 0, ngrid-1) /)	
		index_a  = index_a*int(nhhw1/ngrid)+1
		index_a(ngrid)=nhhw1
		wcg      = w_wc(index_a)
		tg       = w_thetag(index_a)
		wpe2g(:) = wpe2(1) !! if plasma frequency is not uniform, modify this part
		dhg      = index_a(2)-index_a(1) 
		dw1      = (ww_max-ww_min)/(nwf-2)	
		if (ww_max .eq. ww_min) dw1=0.1  ! to avoid problem on constant frequency cases
			
		allocate(sint(ngrid), cost(ngrid), sint2(ngrid), cost2(ngrid), sint4(ngrid), &
				 XX(ngrid), YY(ngrid), Xw(ngrid), Yw(ngrid), Y2(ngrid), Q1(ngrid), Q2(ngrid), rk(ngrid), &
				 RR(ngrid), LL(ngrid), PP(ngrid), SS(ngrid), DD(ngrid))		
		
		do bi = 1, nwf
			w1    = ww_min+dw1*(bi-1); 
			ww_array(bi) = w1
			wf2   = w1*w1
			cost  =  dcosd(tg); cost2 = cost*cost;
			sint  =  dsind(tg); sint2 = sint*sint; sint4 = sint2*sint2;  
		    XX    = wpe2g/(wf2);
			YY    = wcg/w1;
			Xw    = -2*wpe2g/(w1**3);
			Yw    = -wcg/wf2;
			RR    = 1 + XX/(YY-1);
			LL    = 1 - XX/(YY+1);
			PP    = 1 - XX;
			SS    =  0.5d0*(RR+LL);    
			DD    =  0.5d0*(RR-LL); 		

			Y2    = YY*YY;
			Q2    = Y2*sint4+4*(1-XX)**2*cost2;
			Q1    = 2*(1-XX)-Y2*sint2+YY*sqrt(Q2);
			rk_array(bi,:) = sqrt(1-(2*XX*(1-XX)/Q1))*w1/cv;
			kpara_array(bi,:) = rk_array(bi,:)*cost
			kperp_array(bi,:) = rk_array(bi,:)*sint
			
			!----------------- Calculate Vg_para ---------------------------------		
			Kappa  = sint*cost*kperp_array(bi,:)/(kpara_array(bi,:)*kpara_array(bi,:))/(1+(kperp_array(bi,:)/kpara_array(bi,:))**2);
			Kappar = Kappa*kpara_array(bi,:)/kperp_array(bi,:); 
			do ai  = 1,ngrid
			   if (kperp_array(bi,ai) .EQ. 0) then
				 Kappar(ai) = 0;
			   end if
			end do

			dDkpara   = 2*kpara_array(bi,:)*cs/wf2 &
					 - 2*XX*(1-XX)*(Q1**(-2))* &
					 (2*Y2*Kappa + 0.5*YY*(Q2**(-0.5))* &
					 (-4*Y2*sint2*Kappa+8*Kappa*(1-XX)**2));

			dDkperp   = 2*kperp_array(bi,:)*cs/wf2 - 2*XX*(1-XX)*(Q1**(-2))* &
					 (-2*Y2*Kappar + 0.5*YY*(Q2**(-0.5))* &
					 (4*Y2*sint2*Kappar - 8*Kappar*(1-XX)**2));

			dQ2    = 2*YY*Yw*sint4 - 8*(1-XX)*Xw*cost2;      ! if have time, modify this part
			dQ1    =-2*Xw-2*YY*Yw*sint2+Yw*sqrt(Q2)+ &       ! if have time, modify this part
					 0.5*YY/sqrt(Q2)*dQ2;                    ! if have time, modify this part

			dDw    =-2*cs*rk_array(bi,:)*rk_array(bi,:)/(w1**3)+(2*Xw-4*XX*Xw)/Q1-2*XX*(1-XX)/(Q1*Q1)*dQ1; ! if have time, modify this part
			
			Vgpara_array(bi,:) = -dDkpara/dDw;
			Vgperp_array(bi,:) = -dDkperp/dDw;
		end do		
		

		deallocate(sint, cost, sint2, cost2, sint4, &
				   XX, YY, Xw, Yw, Y2, Q1, Q2, rk, &
				   RR, LL, PP, SS, DD)
		allocate(sint(nhh1), cost(nhh1), sint2(nhh1), cost2(nhh1), sint4(nhh1), &
				   rk(nhh1), RR(nhh1), LL(nhh1), PP(nhh1), SS(nhh1), DD(nhh1), &
				   XX(1), YY(1), Xw(1), Yw(1), Y2(1), Q1(1), Q2(1))		   
			
!intepolate initial wave number and group velocity for wave equations 	
	    w1 = ww_max
		w_kpara_max(1)  = kpara_array(nwf-1,1)  		
		w_kpara_max(nhhw1)  = kpara_array(nwf-1,ngrid) 		
	
		!$OMP PARALLEL DO PRIVATE(bi,hgi,dhgh,hjj,hii,ii1) 		
		do bi = 2,int(index_a(ngrid-1))
			hgi        = hgrid(bi);
			dhgh       = hgi-index_a;
			hjj        = minloc(dhgh, dhgh .ge. 0);
			hii        = hjj(1);
			ii1        = hii+1;
			w_kpara_max(bi)  = (kpara_array(nwf-1,ii1)*(hgi-index_a(hii))+kpara_array(nwf-1,hii)*(index_a(ii1)-hgi))/dhg;  		
		end do
		!$OMP END PARALLEL DO
		dhg_end = index_a(ngrid)-index_a(ngrid-1)

		!$OMP PARALLEL DO PRIVATE(bi,hgi,dhgh,hjj,hii,ii1) 	
		do bi = index_a(ngrid-1)+1, nhhw1-1
			hgi        = hgrid(bi);
			dhgh       = hgi-index_a; 
			hjj        = minloc(dhgh, dhgh .ge. 0);
			hii        = hjj(1);
			ii1        = hii+1;
			w_kpara_max(bi)  = (kpara_array(nwf-1,ii1)*abs(hgi-index_a(hii))+kpara_array(nwf-1,hii)*abs(hgi-index_a(ii1)))/dhg_end; 	
		end do		
		!$OMP END PARALLEL DO
		
		w1 = ww_min
		w_kpara(1)  = kpara_array(1,1)  
		w_Vgpara(1) = Vgpara_array(1,1)  		
		w_kpara(nhhw1)  = kpara_array(1,ngrid) 		
		w_Vgpara(nhhw1) = Vgpara_array(1,ngrid)  			
	
		!$OMP PARALLEL DO PRIVATE(bi,hgi,dhgh,hjj,hii,ii1) 		
		do bi = 2,int(index_a(ngrid-1))
			hgi        = hgrid(bi);
			dhgh       = hgi-index_a;
			hjj        = minloc(dhgh, dhgh .ge. 0);
			hii        = hjj(1);
			ii1        = hii+1;
			w_kpara(bi)  = (kpara_array(1,ii1)*(hgi-index_a(hii))+kpara_array(1,hii)*(index_a(ii1)-hgi))/dhg;  	
			w_Vgpara(bi) = (Vgpara_array(1,ii1)*(hgi-index_a(hii))+Vgpara_array(1,hii)*(index_a(ii1)-hgi))/dhg;  	
		end do
		!$OMP END PARALLEL DO
		dhg_end = index_a(ngrid)-index_a(ngrid-1)

		!$OMP PARALLEL DO PRIVATE(bi,hgi,dhgh,hjj,hii,ii1) 	
		do bi = index_a(ngrid-1)+1, nhhw1-1
			hgi        = hgrid(bi);
			dhgh       = hgi-index_a; 
			hjj        = minloc(dhgh, dhgh .ge. 0);
			hii        = hjj(1);
			ii1        = hii+1;
			w_kpara(bi)  = (kpara_array(1,ii1)*abs(hgi-index_a(hii))+kpara_array(1,hii)*abs(hgi-index_a(ii1)))/dhg_end; 
			w_Vgpara(bi) = (Vgpara_array(1,ii1)*abs(hgi-index_a(hii))+Vgpara_array(1,hii)*abs(hgi-index_a(ii1)))/dhg_end;  			
		end do		
		!$OMP END PARALLEL DO
		
! resample the wave number and group velocity bk file
		kpara  = w_kpara(index_w) 
		Vgpara = w_Vgpara(index_w) 
        kpara_max  = w_kpara_max(index_w) 		

! Initial wave profile------------------------------------------------
        ! wave phase
        psi(1) = 0;   
        do pii = 2,nhh1  
			psi(pii) = mod(psi(pii-1)-(kpara(pii-1)+kpara(pii))*0.5, tpi);
		end do 	

! (2.16) of Yuko's thesis ----------------------------------------------------------
        Hje = 0.0d0
		Hjb = 0.0d0
		do pii = 1, nSbb
			 sss    = 1.0d0/dble(nSbb)*dble(pii-1.0d0) ! Assuming that S>0
			 zeta1 =-dasin(sss)+tpi
			 a1    = dcos(zeta1)-sss*zeta1
			 do hii = 1, nzeta
				  a2 =a1-(dcos(zeta1-zeta1/dble(nzeta)*dble(hii))-sss*(zeta1-zeta1/dble(nzeta)*dble(hii)))
				  if (dsign(1.0d0, a2) < 0.0d0)exit
			 end do
			 zeta2 = zeta1-zeta1/dble(nzeta)*dble(hii-1.0d0)
			 do hii = 1, nzeta
				  zeta = zeta1+(zeta2-zeta1)/dble(nzeta)*dble(hii)
				  Hje(pii)=Hje(pii)+dsin(zeta)*dsqrt(dcos(zeta1)-dcos(zeta)+sss*(zeta-zeta1))*(zeta2-zeta1)/dble(nzeta)
				  Hjb(pii)=Hjb(pii)+dcos(zeta)*dsqrt(dcos(zeta1)-dcos(zeta)+sss*(zeta-zeta1))*(zeta2-zeta1)/dble(nzeta)
			 end do
		end do						
!====================================================================
! Particles
!====================================================================		

		npi = istart-1
		npp = 0
		npj = 0
	
		call random_seed(size=seedsize)
        allocate(seed(seedsize))
		allocate(ranran(6,np), wave_ran(200))			 
		seed(:)=(myrank+1)*seednum
        CALL random_seed(put=seed)
        CALL random_number(ranran)
		
		do ienergy = 1, nenergy
			do ipitch = 1, npitch
				do ikeido = 1, nkeido
					do iphase = 1, nphase	  
					
						if(npj+1 > npi) then
							npp = npp+1
							npj = npj+1
							if(npp > np) exit
							check(npp) = 0
							
							if(jobne == 1) then ! initial kinetic energy
								KE0 = energy_min+energy_del*dble(nenergy-1)*ranran(2,npp)		 
							else 
								KE0 = energy_min+energy_del*dble(ienergy-1)
							end if
							if(jobnp == 1) then ! initial equatorial pitch angle
								PA0 = pitcheq_min+(pitcheq_max-pitcheq_min)*ranran(3,npp)
							else
								PA0 = pitcheq_min+pitcheq_del*dble(ipitch-1)
							end if
							if(jobnk==1)then
								Long0 = keido_min+(keido_max-keido_min)*ranran(4,npp)
							else
								Long0 = keido_min+keido_del*dble(ikeido-1)
							end if
							if(joblo == 1) then ! initial location within -1~1 
								Locat0  = -1+2*ranran(5,npp)
							else
								Locat0  = 0; 
							end if 										
							if(jobng == 1) then ! initial gyro phase
								Pha0    = tpi*ranran(6,npp)
							else
								dphi    = tpi/dble(nphase)
								Pha0    = dphi*dble(iphase)
							end if
 	
							KEr     = KE0 *1.956d0*1e-3
							Vt      = sqrt(cs*(1.0d0-1.0d0/(KEr+1.0d0)**2))
							! find mirror point
							Vperp0  = Vt*sind(PA0)
							mpb     = Beq/(dsind(PA0)**2); ! background magnetic field at the mirror point
							idx_L1  = minloc( bbb, bbb >= mpb); mp_idx  =   idx_L1(1)-1
							if (mp_idx < 1) then
								mp_idx = 1;
							end if
							
							hloc    = nint(dble(mp_idx-1)*Locat0) 													
							px(npp) = xg(cet+hloc)*dcosd(Long0) 
							py(npp) = xg(cet+hloc)*dsind(Long0)
							pz(npp) = zg(cet+hloc)
							hh(npp) = hg(cet+hloc)
							
							vperp(npp) = dsqrt(bbb(1+iabs(hloc))*Vperp0**2/Beq)
							vpara(npp) = dsqrt(Vt**2-vperp(npp)**2)   !times -1 if south hemistphere  		
													
							if( isnan(vpara(npp)))then
								vperp(npp) = v0
								vpara(npp) = 0.0d0
							end if
 	 							
							r2 = px(npp)*px(npp) +py(npp)*py(npp) + pz(npp)*pz(npp) 							
							wratio = -3.0d0*Beq*radeq3/(r2*r2*dsqrt(r2))
 							
							bx0 = wratio*px(npp)*pz(npp)
							by0 = wratio*py(npp)*pz(npp)
							bz0 = wratio*(pz(npp)*pz(npp)- r2*threei)

							t_eta = dsqrt(bx0*bx0+by0*by0)/abs(bz0) ! tan(eta)
							c_eta = 1.0/dsqrt(1.0+ t_eta*t_eta)       ! cos(eta)
							s_eta = c_eta*t_eta				        ! sin(eta)

							vx(npp) = (-dsign(1.0d0, pz(npp))*vpara(npp)*s_eta+ dsign(1.0d0, bz0)*vperp(npp)*dcos(Pha0)*c_eta)*dcosd(Long0)-vperp(npp)*dsin(Pha0)*dsind(Long0)
							vy(npp) = (-dsign(1.0d0, pz(npp))*vpara(npp)*s_eta+ dsign(1.0d0, bz0)*vperp(npp)*dcos(Pha0)*c_eta)*dsind(Long0)+vperp(npp)*dsin(Pha0)*dcosd(Long0)
							vz(npp) = dsign(1.0d0, bz0)*vpara(npp)*c_eta + vperp(npp)*dcos(Pha0)*dsign(1.0d0, pz(npp))*s_eta
							ccc = 1.0d0/dsqrt(bx0*bx0 + by0*by0 + bz0*bz0)
							vpara(npp) = (vx(npp)*bx0 + vy(npp)*by0 + vz(npp)*bz0)*ccc
							vperp(npp) = dsqrt(vx(npp)**2 + vy(npp)**2 + vz(npp)**2 - vpara(npp)**2) 				
							if( isnan(vperp(npp)))then
								vperp(npp)=0.0d0
							end if									
				! check if vperp/vpara are the same
							KE(npp) = KE0	
							PA(npp) = asin(dsqrt(bbb(hloc)/Beq)*dsind(PA0))      ! local pitch angle
							
						else
							npj = npj+1
						end if						
					end do
					if(npp > np) exit
				end do	
				if(npp > np) exit
			end do
			if(npp > np) exit
		end do		
 
!wave*************************************************************************	  
	  seed(:)=100
      call random_seed(put=seed)
      call random_number(wave_ran)
	  wsign = 1.0d0
	  thranind = 2
	  
	  
	  ! Generate Gaussian-shape wave packet
	  wave_eq_end  = 1.4*tlength/real_dt ! convert wave end time from second to the time step  
	  ta           = (tlength/2+ta)/real_dt 
	  sigma        = sigma/real_dt
	  Bwmax        = Bwmax*Beq
	  
!diagnostic----------------------------------------------------
     jt_wave=0
     jt_trajectory=0
     iv=0 
	 iv2=0
 
        RETURN	
      end