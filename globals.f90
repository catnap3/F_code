!--- Module for declaring------------ 
module globals
	use omp_lib	
	implicit none

!====================================================
! simulation control
!====================================================	
    integer(kind=4),parameter :: test      = 0 ! 0 = ordinary; 1 = detailed trajectory
	integer(kind=4),parameter :: waveonly  = 1 ! 0 = with particles; 1 = waveonly; 2= equator waveonly
	integer(kind=4),parameter :: waveskip  = 1

!====================================================
! input parameters
!====================================================	
	! simulation parameters
	integer(kind=4), parameter :: nhnum    = 9836   ! number of grid along the geomagnetic field line. The numbere is according to L-value
	integer(kind=4), parameter :: sgn      = -1     ! sign of particles : -1 = electron; 1 = ion
	integer(kind=4)            :: nphase   = 36     ! number of gyrophases
	
    double precision, parameter :: Lv       = 4.7d0  ! L-values
    double precision, parameter :: dt       = 0.01d0 ! time step [1/wceq]
	double precision, parameter :: dh       = 1.0d0  ! grid size for others [c/wceq]  	
	double precision, parameter :: dhw      = 0.5d0  ! grid size for wave subpackets [c/wceq]  	
	double precision, parameter :: lifetime = 0.6944d0! wave duration time [sec] @waveonly==0
	                                          ! please set up the parameter according as the result@waveonly==1	
											  
	double precision :: wna    = 60.0d0                 ! maximum wave normal angle [deg]
	double precision :: wceq   = 1.0d0                  ! normalized electron gyro frequency at equator
	double precision :: wpeq   = 3.7503d0               ! electron plasma frequency [wceq]
	double precision :: qm     = -1.0d0                 ! normalized charge/mass of particle
	double precision :: m0     = 1.0d0                  ! normalized mass of particle
	double precision :: e0     = 1.0d0                  ! normalized charge of particle
    double precision :: cv     = 1.0d0                  ! light speed   	
	double precision :: Bw_ini = 0.000001d0             ! initial wave amplitude of B [Beq]
	double precision :: ww_min = 0.5d0                 ! wave frequency range [wceq]
	double precision :: ww_max = 0.5d0                  ! wave frequency range [wceq]
	
	! wave generation
	double precision :: vpro      = 0.1d0  			  !averaged perpendicular velocity for hot electrons [c]
	double precision :: Vtpara    = 0.0692d0 		      !thermal parallel velocity [c]
	double precision :: Q         = 0.5d0
    double precision :: tau       = 0.5d0
    double precision :: threshold = 1.001d0
    double precision :: bwmod     = 0.4d0
	double precision :: N_ph      = 0.0018d0            ! source hot electron density at equator normalized by cold electron density [ne]
	double precision :: N_phc     = 0.0003d0           ! source hot electron density at convective growth region normalized by cold electron density [ne]
    double precision :: Bcheck    = 1.0d-8             ! check for convective growth [Beq]
    double precision :: Bcut      = 0.015d0            ![Beq]
	
	double precision :: tlength = 0.12     ! duration of a subpacket [sec]		
    double precision :: Bwmax   = 5*1e-4   ! maximum wave amplitude [B0eq] 
	double precision :: ta    = 0.02d0     ! time shift from the beginning of the simulation [sec]
	double precision :: sigma = 0.015d0    ! standard deviation [sec]
	integer(kind=4)  :: wave_eq_end 
	
	
    ! set up the RADIATION BELT ELECTRON parameter
	integer(kind=4) :: jobng = 1                 ! gyro phase :                0 = interval, 1=random
	integer(kind=4) :: jobnp = 1                 ! pitch angle :               0 = interval, 1=random
	integer(kind=4) :: joblo = 1                 ! location along field line : 0 = equator, 1=random
	integer(kind=4) :: jobne = 1                 ! kinetic energy :            0 = interval, 1=random
	integer(kind=4) :: jobnk = 1                 ! longitudinal direction :    0 = interval, 1=random

	! Green's function parameters
	integer(kind=4), parameter :: nSbb    = 100   ! mesh of whistler Sbb parameter !Sbb(R,alpha_eq) given by (10) of Omura et al. 2015
	integer(kind=4), parameter :: nzeta   = 300   ! total integral number of zeta
	integer(kind=4), parameter :: nwf     = 2     ! mesh of wave vectors & group velocities array in wave frequencies
	integer(kind=4), parameter :: ngrid   = 1000  ! mesh of wave vectors & group velocities array along grids
	integer(kind=4), parameter :: seednum = 404  ! seed number for generating randam number for different runs
	integer(kind=4), parameter :: strank  = 0     ! filename contral number. Use it if you have many runs in a folder. The strank should be the number of the largest filenumber of the last run.
	double precision :: energy_min  = 10.0d0         ! minimum energy [keV]
	double precision :: energy_max  = 20.0d0       ! maximum energy [keV]
	double precision :: energy_del  = 10.0d0         ! energy interval for Green's function [keV]
	double precision :: pitcheq_min = 1.0d0          ! minimum equatorial pitch angle [deg]
	double precision :: pitcheq_max = 10.0d0         ! maximum equatorial pitch angle [deg]
	double precision :: pitcheq_del = 1.0d0          ! pitch angle interval for Green's function [deg]
	double precision :: keido_min   = 0.0d0          ! minimum longitude [deg]
    double precision :: keido_max   = 360.0d0        ! mzximum longitude [deg]
    double precision :: keido_del   = 36.0d0         ! logitude interval for Green's function [deg] @jobnk==0 
	double precision :: whistler_kakudo_min  = 0.0d0  
    double precision :: whistler_kakudo_max  = 360.0d0
	
	! constants--------------------------------------
	double precision, parameter   :: pi      = 4.0*atan(1.0d0)
	double precision, parameter   :: qe      = 1.602E-19        ! Electric charge [A*s]
	double precision, parameter   :: me      = 9.1E-31          ! Mass of an electron [kg]
	double precision, parameter   :: B0E     = 3.12E-5          ! Magnetic field at the euqator on the Earth's surface [T]
	double precision, parameter   :: RE      = 6378000          ! Earth radius [m]
	double precision, parameter   :: real_cv = 3.0d8            ! Light speed [m/s]
	double precision, parameter   :: real_Beq= 2.4012d-7!B0E/(Lv*Lv*Lv)        ! set real magnetic field at the equator of the L-value [T] 
	double precision, parameter   :: real_wce= qe*real_Beq/me        ! set real electron gyro-frequency at the equator with L-value [rad]
	
	! relatived numbers
	integer(kind=4), parameter :: nhp1 = nhnum+1       !nhnum +1(grid number +1)
	integer(kind=4), parameter :: nhh  = int(nhnum/2)  !nhnum half
	integer(kind=4), parameter :: nhh1 = nhh+1
	integer(kind=4), parameter :: nhw  = nhnum*int(dh/dhw)+1 ! for wave equations	
	integer(kind=4), parameter :: nhhw = int(nhw/2)
	integer(kind=4), parameter :: nhhw1= nhhw+1
	double precision, parameter :: dx  = 1.0d0
    integer(kind=4), parameter :: nx   = int(90/dx)+1
	

!====================================================
! main.f90 & MPI & random number
!====================================================	
    ! loop ------------------------------------------
	integer(kind=8) :: itime, ai, bi, bj, npi, mi, pii, i, j, ii, jj, m
	integer(kind=8) :: aaa, swave
	double precision :: a1, a2, b1, bb, time
	! MPI -------------------------------------------
	integer(kind=8) :: ierr, nprocs, myrank
	integer(kind=8) :: np, istart,iend 
	! random number ---------------------------------
	integer(kind=8) :: seedsize
	integer(kind=8), allocatable, dimension(:) :: seed
	double precision, allocatable, dimension(:,:) :: ranran	
	! character--------------------------------------
	character(len = 10) :: cols ! type the digit
	character(len = 10) :: nzs		
	character(len = 10) :: folder
	! output
	double precision :: pitchangle, zz

!====================================================
! rescaling.f90
!====================================================
	double precision, save  :: rnh, rnt, rnv, rne, rnb, rnj	
	double precision, dimension(nhh1) :: rkpara, rkperp, rkpara_max, rwc, rB0, rVgpara
	
!====================================================
! initial.f90
!====================================================		
	! general used ----------------------------------
	integer(kind=8) :: npj, npp, iv, iv2 !important
	integer(kind=8) :: ntime
	integer(kind=8) :: cet, nhh2
	integer(kind=8) :: nenergy, nenergy1, npitch, nkeido, totalp
	integer(kind=8) :: ii1, hii, hjj(1)
	integer(kind=8) :: pskip_trajectory, pskip_h, pskip_wave, p_trajectory, p_wave
	integer(kind=8) :: jt_wave, jt_trajectory, iconvective
	integer(kind=8) :: tra_time, tra_start
	double precision    :: tra_start1
    double precision    :: cs, tcs, tpi, mu0, eps0, wpes, aqm, REr
	double precision    :: Beq, ab, aa
	double precision    :: ddz, ddz1, dhww
	!double precision    :: Bweq, wweq, PAeq, Aeq	
	double precision    :: radeq, g_pe, g_pe1
	double precision    :: w_ph, w_phs, w_phc, w_phcs, AAAeq
	double precision    :: xbh, zbh, xb, zb, r2, wratio, bx0, by0, bz0, bab, &
					zeta,zeta1,zeta2, &
					Kinetic,Long, v0, pit, Vperp0,  dphi, dphim
	double precision    :: polar_latitude, atomos, wsign, real_dt, real_dh
	
	integer(kind=8) :: index_w(nhh1), rhw
	integer(kind=8) :: ienergy, ipitch, iphase, ilocat, ikeido
	
! for wave equations	
	double precision, dimension(nhhw1) :: w_wf, w_wfn, w_Bw, w_Bwn, w_wc, w_bbb !half grids for wave equations
	double precision, dimension(nhhw1) :: w_kpara, w_kpara_max, w_Vgpara, w_thetag !half grids for wave equations
	double precision, dimension(nhhw1) :: w_zg, w_xg, w_hg, w_lambdag !half grids for wave equations
	double precision :: w_rk
	
	double precision, dimension(nhh1) :: wf, Bw, wc, wpe2, bbb, lambdag !half grids
	double precision, dimension(nhh1) :: thetag
	double precision, dimension(nhh1) :: kperp, kpara, kpara_max, psi, cospsi2, sinpsi2
	double precision, dimension(nhh1) :: dBz, wcyclo, dwc
	double precision, dimension(nhh1) :: Vgpara, As, Ap, rn2
	double precision, dimension(nhh1) :: Bxw, Byw, Bzw, Exw, Eyw, Ezw
	
	double precision, dimension(nhp1) :: xg, zg, hg ! All grids

	double precision, dimension(nSbb) :: Hje, Hjb
	double precision, dimension(nx) :: Sbb
! for interpolation of wave number and group velocity
	integer(kind=8) :: wii, wjj(1), wi1
	double precision :: w1, wf2, dw1, dhg, dhg_end, dhgg, hgi
	double precision :: rk_a, rk_b, Vgpara_a, Vgpara_b, Vgperp_a, Vgperp_b
	double precision :: w_a, w_b, hg_a, hg_b
	integer(kind=4), dimension(ngrid) :: index_a
	double precision, dimension(ngrid) :: tg, wcg, wpe2g, dhgh
	double precision, dimension(ngrid) :: dQ1, dQ2
	double precision, dimension(ngrid) :: Kappa, Kappar, dDkpara, dDkperp, dDw 
	double precision, dimension(nwf,ngrid) :: rk_array, kperp_array, kpara_array
	double precision, dimension(nwf,ngrid) :: Vgpara_array, Vgperp_array
	double precision, dimension(nwf) :: ww_array, dwf
	double precision, dimension(nhhw1) :: hgrid
	double precision, allocatable, dimension(:) :: sint, cost, sint2, cost2, sint4, &
												XX, YY, Xw, Yw, Y2, Q1, Q2, rk, &
												RR, LL, PP, SS, DD
														
	double precision, allocatable, dimension(:) :: wave_ran
	
! for output
    double precision :: tp, tp_wave
	double precision, allocatable, dimension(:) :: px, py, pz, hh, vx, vy, vz,&   ! px: position x
										 vpara, vperp, energy, phi, trajectoryrx, trajectoryry,trajectoryrz, Bmoment
	double precision, allocatable, dimension(:,:) :: wfBp, BwBp, Keq_first, Keq
    double precision, allocatable, dimension(:)   :: opt, thre      
    double precision, allocatable, dimension(:) :: trajectory, trajectorypitch, &
                                           trajectorylambda, trajectoryz,trajectoryst,trajectorye
    double precision, allocatable, dimension(:) :: tp_Green, check	
	double precision, allocatable, dimension(:) :: PA, KE
	
	! used only in initial.f90 ----------------------
    integer(kind=8) :: Li
	integer(kind=8) :: idx_L1(1), idx_L2(1), hloc, mp_idx
	double precision, dimension(nhh1) :: Bp, Bn, thetaL, Lat, dhL
	double precision  :: zunit, radeq3, threei
    double precision  :: Bx00, By00, Ex00, Ey00, Ez00
    double precision  :: kzeq, Vt, KEr, pangle, ddh, inth
	double precision  :: KE0, PA0, Locat0, Long0, Pha0
	double precision  :: Vxy0, mpb, mp
	double precision  :: wnaL, fpl, ccc

!====================================================
! diagno
!====================================================
    integer(kind=8) :: nh_end
    
!====================================================
! vdist
!====================================================
    character(len = 5) :: cy_rank
    integer(kind=4) :: fi
	double precision    :: Vperpeq, Vparaeq
    double precision    :: rw_min, rw_max, rc, rwp
	
!====================================================
! velocity.f90
!====================================================		
	double precision :: gamma, signal, bh, rL, cosphip	
	double precision :: vpara_x, vpara_y, vpara_z, vperp_x, vperp_y, vperp_z
	double precision :: Rgc_x, Rgc_y, Rgc_z
	double precision, dimension(3) :: gca, xpp1, xpp
	double precision :: wf1, wc1, kpara1, kperp1
	double precision :: eta, c_eta, s_eta, t_eta
	double precision :: long_phi, c_phi, s_phi, &
                    psi1, cpsi1, spsi1, As1, Ap1,&
					Bw1, Bxw1, Byw1, Bzw1, Exw1, Eyw1, Ezw1, &
                    bwr, bwt, bwr0, bwt0, ewr, ewt, &
					bxwx, bywx, bzwx, bxwy, bywy, bzwy, &
					bx2, by2, bz2, ex1, ey1,ez1, &
					ux, uy, uz, bx1, by1, bz1, uxt, uyt, uzt, boris
	double precision :: bx0r, by0r, bz0r
	
!====================================================
! oblique_whistler.f90
!====================================================	
    integer(kind=8) :: wai, hi, ig, ig2, thranind
	double precision :: Wh, vp, VR1, s0, s1, s2, sss, wtr2, vtr, &
						Hje1, Hjb1, J0, Je, Jb, GG, &
						chi, xi, gamma_wave, w_i, w_in,&
				        dwt, dBt, opt1, thre1, thrand
    double precision :: Bop, Bth, Bthrand, ddx							  
    double precision :: wfg2
	
end module globals