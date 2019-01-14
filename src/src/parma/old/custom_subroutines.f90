!  Calculate cosmic-ray fluxes in the atmosphere based on PARMA model
subroutine getParmaSpec(ip,ener,cosang,alti,glat,glong,iyear,imonth,iday,DifFlux)
        USE parma, ONLY : getSpec, getSpecAngFinal,getHP,getr,getd
	implicit none
	integer, parameter :: npart=33 ! number of applicable particle
	integer :: IangPart(0:npart)

	data IangPart/1,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,5,5,6/ ! Angular particle ID

	real*8, intent(in) :: ener,glat,glong,alti,cosang
	integer*4, intent(in) :: ip,iyear,imonth,iday
	real*8, intent(out) :: DifFlux
	real*8 :: g=0.,Flux=0.,s=0.,r=0.,d=0.
	integer :: idummy=0


	!ip=33 ! Particle ID (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
	!e=energies(jj)    ! Energy (MeV/n)

	!iyear=2016  ! Year
	!imonth=3    ! Month
	!iday=1      ! Date
	!glat=-13.0   ! Latitude (deg), -90 =< glat =< 90
	!glong=130.0 ! Longitude (deg), -180 =< glong =< 180
	!Alti=12.0   ! Altitude (km)
	g=10.0      ! Local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin
	!ang=cos_angles(kk)     ! cosine of zenith angle (e.g. ang=1.0 for vertical direction, ang=0.0 for holizontal direction)

	s=getHP(iyear,imonth,iday,idummy) ! Solar activity (W index) of the day
	r=getr(glat,glong)                ! Vertical cut-off rigidity (GV)
	d=getd(alti,glat)                 ! Atmospheric depth (g/cm2), set glat = 100 for use US Standard Atmosphere 1976.

        if (isnan(r)) r = 14.0
        write(*,*) ip,s,r,d,ener,g, getSpec(ip,s,r,d,ener,g)
	!write(6,*) 'Angular Integrated Flux(/cm2/s/(MeV/n))=',Flux
	!if(IangPart(ip).ne.0) then ! Angular distribution is available for the particle
	 DifFlux=getSpec(ip,s,r,d,ener,g)*getSpecAngFinal(iangpart(ip),s,r,d,ener,g,cosang)

         write(*,*) 'ERROR : DifFlux is nan. Aborting'
         if (isnan(DifFlux)) call abort
	 !write(6,*) 'Angular Differential Flux(/cm2/s/(MeV/n)/sr)=',DifFlux

end subroutine getParmaSpec
