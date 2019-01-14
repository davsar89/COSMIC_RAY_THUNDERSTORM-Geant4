
use lib_array, only : logspace, linspace


!  Calculate cosmic-ray fluxes in the atmosphere based on PARMA model


parameter(npart=33) ! number of applicable particle
implicit real*8 (a-h, o-z)
dimension IangPart(0:npart)
data IangPart/1,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,5,5,6/ ! Angular particle ID
integer, parameter :: NB_ENERS=128
integer, parameter :: NB_ANGLS=64
integer,parameter :: idp = selected_int_kind(13)
integer,parameter :: sp = selected_real_kind(p=6,r=37)
integer,parameter :: dp = selected_real_kind(p=15,r=307)

real(dp) :: energies(NB_ENERS)
real(dp) :: xmin=0.01
real(dp) :: xmax=1000000.0
real(dp) :: cos_angles(NB_ANGLS)
real(dp) :: ang_min=-1.0
real(dp) :: ang_max=1.0


call logspace(xmin,xmax,energies)
call linspace(ang_min,ang_max,cos_angles)

	do jj = 1, size(energies)
		do kk = 1, size(cos_angles)

			ip=33 ! Particle ID (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
			e=energies(jj)    ! Energy (MeV/n)
			iyear=2016  ! Year
			imonth=1    ! Month
			iday=1      ! Date
			glat=-13.0   ! Latitude (deg), -90 =< glat =< 90
			glong=130.0 ! Longitude (deg), -180 =< glong =< 180
			Alti=12.0   ! Altitude (km)
			g=0.0      ! Local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin
			ang=cos_angles(kk)     ! cosine of zenith angle (e.g. ang=1.0 for vertical direction, ang=0.0 for holizontal direction)

			s=getHP(iyear,imonth,iday,idummy) ! Solar activity (W index) of the day
			r=getr(glat,glong)                ! Vertical cut-off rigidity (GV)
			d=getd(alti,glat)                 ! Atmospheric depth (g/cm2), set glat = 100 for use US Standard Atmosphere 1976.

			Flux=getSpec(ip,s,r,d,e,g)
			!write(6,*) 'Angular Integrated Flux(/cm2/s/(MeV/n))=',Flux
			!if(IangPart(ip).ne.0) then ! Angular distribution is available for the particle
			 DifFlux=Flux*getSpecAngFinal(iangpart(ip),s,r,d,e,g,ang)
			 write(6,*) 'Angular Differential Flux(/cm2/s/(MeV/n)/sr)=',DifFlux
			!endif
		end do

end do

end
