!  Calculate cosmic-ray fluxes in the atmosphere based on PARMA model


parameter(npart=33) ! number of applicable particle
implicit real*8 (a-h, o-z)
dimension IangPart(0:npart)
data IangPart/1,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,5,5,6/ ! Angular particle ID


			ip=31 ! Particle ID (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
			e=1.0    ! Energy (MeV/n)
			iyear=2016  ! Year
			imonth=3    ! Month
			iday=1      ! Date
			glat=-13.0   ! Latitude (deg), -90 =< glat =< 90
			glong=130.0 ! Longitude (deg), -180 =< glong =< 180
			Alti=12.0   ! Altitude (km)
			g=0.0      ! Local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin
			ang=-1.0     ! cosine of zenith angle (e.g. ang=1.0 for vertical direction, ang=0.0 for holizontal direction)

			s=getHP(iyear,imonth,iday,idummy) ! Solar activity (W index) of the day
			r=getr(glat,glong)                ! Vertical cut-off rigidity (GV)
			d=getd(alti,glat)                 ! Atmospheric depth (g/cm2), set glat = 100 for use US Standard Atmosphere 1976.

			Flux=getSpec(ip,s,r,d,e,g)
			!write(6,*) 'Angular Integrated Flux(/cm2/s/(MeV/n))=',Flux
			!if(IangPart(ip).ne.0) then ! Angular distribution is available for the particle
			 DifFlux=Flux*getSpecAngFinal(iangpart(ip),s,r,d,e,g,ang)
			 write(6,*) 'Angular Differential Flux(/cm2/s/(MeV/n)/sr)=',DifFlux
			!endif

end
