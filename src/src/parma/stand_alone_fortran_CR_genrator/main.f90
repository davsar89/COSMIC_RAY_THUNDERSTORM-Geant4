!  Calculate cosmic-ray fluxes in the atmosphere based on PARMA model
implicit real*8 (a-h, o-z)
parameter(nebin=140)
parameter(npart=33) ! number of particle
parameter(nang=101)  ! from -1 to 1, 0.2 step
dimension emid(nebin+1),ewid(nebin+1) ! nebin+1 for energy integrated
dimension dcc(0:npart,nebin) ! dose conversion coefficient (pSv*cm^2)
dimension Flux(0:npart,nebin+1) ! energy differential flux, ie=nebin+1 for energy integrated
dimension dose(0:npart+1)    ! dose from each particle, npart+1 for total
dimension ang(nang) ! angle, 0 for angular integrated
dimension FluxAngEdif(0:npart,nebin+1,nang) ! Energy differntial & Angular differential flux
dimension IangPart(0:npart)
character fname*40,pname(0:npart)*6
character chatmp*1,chatmp1*1,chatmp2*2
character condition*40,chatmp40*40
character dccname*7
character uname*6

data IangPart/1,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,5,5,6/ ! Angular particle ID

data pname/'neutro','proton','he---4','li---7','be---9','b---11','c---12','n---14','o---16','f---19','ne--20', & 
    &      'na--23','mg--24','al--27','si--28','p---31','s---32','cl--35','ar--40', & 
    &      'k---39','ca--40','sc--45','ti--48','v---51','cr--52','mn--55','fe--56','co--59','ni--59', &
    &      'muplus','mumins','electr','positr','photon'/
data condition/'Tokyo-Smin'/ ! Default input condition name
 data dccname/'ICRP116'/     ! Effective dose for ISO irradation (uSv/h)
! data dccname/'NM06tub'/     ! Count rate of 6 tube neutron monitor (100 count per hour)
! data dccname/'Air-pGy'/     ! Dose rate in air (uGy/h)
! data dccname/'h10ICRP'/     ! H*(10) (uSv/h)

data ie511/78/ ! energy bin for 511keV photon

!    Read Condition Name
write(6,*) 'Input Condition Input file name'
write(6,*) 'Default= ',condition
read(5,'(a40)') chatmp40
if(chatmp40.ne.'                                        ') condition=chatmp40
do i=1,40
 if(condition(i:i).eq.' ') exit
enddo
length=i-1  ! number of character

open(unit=1,file='Condition/'//condition(1:length)//'.inp',status='old')
read(1,*) isout,istyle 
! isout=0: No output for flux data for each condition, 
!      =1: Output Angular integrated flux in "SpecOut" folder
!      =2: Output Angular differential flux in "Angout" folder as a function of energy
!      =3: Output Angular differential flux in "Angout" folder as a function of cos(theta)
!      =4: Output Angular differential flux in "Angout" folder as a function of cos(theta) (absolute value)
!istyle =-1:Direct input (s****-r***-d*****-g***)
!       = 0:Direct input (s****-r***-d****-g***)
!       = 1:year, month, day, latitude, longitude, atmospheric detph(g/cm^2), g(direct input)
!       = 2:year, month, day, latitude, longitude, altitude (m), g(direct input) 
!       = 3:year, month, day, latitude, longitude, altitude (ft), g(direct input)
!       = 4:year, month, day, cutoff rigidity, altitude (m), g(direct input)

if(istyle.ge.1) itimedep=1 ! time-dependence output
if(isout.eq.1.and.istyle.ne.0) then ! output all spectrum in one file
 open(unit=2,file='SpecOut/'//condition(1:length)//'.out')
 write(2,*) '-s6 -dline'
endif

if(isout.ge.2) then
 do ip=0,npart
  if(iangpart(ip).ne.0) open(100+ip,file='AngOut/'//condition(1:length)//'-'//pname(ip)//'.out')
 enddo
endif

open(unit=4,file='DoseOut/'//condition(1:length)//'.'//dccname)
write(4,*) '-s3'
write(4,'(2i8,i12,i8,35i13)') (ip,ip=0,npart+5)
write(4,1000) (pname(ip),ip=0,npart),uname
1000	format(' W-index COR(GV) Dep(g/cm^2)   Local',34(7x,a6),'        Total ',a6)

!     Read Dose Conversion Coefficient Data
open(unit=3,file='dcc/'//dccname//'.inp',status='old')
read(3,*) chatmp1
read(3,*) chatmp1
do ie=1,nebin
 read(3,*) emid(ie),ewid(ie),(dcc(ip,ie),ip=0,npart)
enddo
close(unit=3)
emid(nebin+1)=0.0d0 ! energy integrated
ewid(nebin+1)=1.0d0 ! energy integrated

if(dccname.eq.'NM06tub') then ! for neutron monitor count rate
 uname='100cph'
 unitconv=3600.0/100.0
else ! for dose
 uname='uSv/h '
 unitconv=1.0e-6*3600
endif

! Set angle
ang(1)=-1.0d0
do ia=2,nang
 ang(ia)=ang(ia-1)+2.0/(nang-1)
 if(abs(ang(ia)).lt.1.0e-10) ang(ia)=0.0d0
enddo

dose(:)=0.0 ! initialization

!     Main Calculation Start!!
do i=1,100000
 10 if(istyle.le.0) then ! classical style
  read(1,'(a40)',end=999) fname
  do ii=1,40
   if(fname(ii:ii).eq.' ') exit
  enddo
  lengthF=ii-1  ! number of character
  read(fname,'(a)') chatmp1
  if(chatmp1.ne.'s') goto 10 ! skip data
  if(istyle.eq.0) then
   read(fname,1010) chatmp1,is,chatmp2,ir,chatmp2,id,chatmp,chatmp1,ig
   dtmp=1.0
  else
   read(fname,1012) chatmp1,is,chatmp2,ir,chatmp2,id,chatmp,chatmp1,ig
   dtmp=0.1
  endif
  s=is*1.0 ! in MV
  r=max(0.0,ir*0.1)
  d=id*dtmp
  if(chatmp.ne.'-') then ! use 0.1 digit
   read(chatmp,'(i1)') idtmp
   d=d+idtmp*0.1
  endif
  g=ig*0.01d0
  if(g.ge.-0.02d0.and.g.lt.0.0) g=10.0d0    ! for semi-infinite atmosphere, e.g. g-01
  if(g.lt.-0.02d0) g=100.0d0 ! for Black hole, e.g. g-10
  if(chatmp1.eq.'p') g=-g        ! for pilot location
  if(chatmp1.eq.'c') g=-g-10.0d0 ! for cabin location 
 else ! specify year, latitude, longitude directly
  if(istyle.ne.4) then
   read(1,*,end=999) iyear,imonth,iday,cido,ckei,d,g ! time, latitude, longitude, alti(or depth), local
  else
   read(1,*,end=999) iyear,imonth,iday,r,d,g ! time, COR, altitude, local
  endif
  s=getHP(iyear,imonth,iday,idummy) ! solar modulation potential in MV for the day
  if(idummy.ge.4) then
   write(6,*) 'Error in getHP, error code =',idummy
   write(6,*) 'yy/mm/dd',iyear,imonth,iday
   goto 5
  endif
  if(istyle.ne.4) then
   r=getr(cido,ckei)
  else
   cido=-100.0 ! identification that latitude is not specified
  endif
  if(istyle.ge.2) then ! convert from m to g/cm^2
   if(istyle.eq.3) d=d*0.3048 ! convert from ft to m
   dtmp=d*0.001 ! altitude in km
   d=getd(dtmp,cido)
  endif
 endif
 do ie=1,nebin 
  e=emid(ie)
  do ip=0,npart
   Flux(ip,ie)=getSpec(ip,s,r,d,e,g)
   if(ip.eq.npart.and.ie.eq.ie511) Flux(ip,ie)=Flux(ip,ie)+get511flux(s,r,d)/ewid(ie) ! add 511 keV flux
   Flux(ip,nebin+1)=Flux(ip,nebin+1)+Flux(ip,ie) ! energy integrated flux
   dose(ip)=dose(ip)+Flux(ip,ie)*dcc(ip,ie)*ewid(ie)  ! Each Particle Dose (pSv/s)
   dose(npart+1)=dose(npart+1)+Flux(ip,ie)*dcc(ip,ie)*ewid(ie)  ! Total Dose (pSv/s)
  enddo
 enddo
 if(isout.eq.1) then
  if(istyle.eq.0) then
   open(unit=2,file='SpecOut/'//fname(1:lengthF)//'.out')
   write(2,*) '-s5 -dline'
  else
   write(2,1015) s,r,d,g 
  endif
  write(2,1020) dose(npart+1)*unitconv,uname
  write(2,1025) (ip,ip=0,npart+1)
  write(2,1030) (dose(ip)*unitconv,ip=0,npart),uname
  write(2,1035) (dose(ip)/dose(npart+1)*100.0,ip=0,npart)
  write(2,1040) (pname(ip),ip=0,npart)
  do ie=1,nebin
   write(2,1050) emid(ie),(Flux(ip,ie),ip=0,npart)
  enddo
 elseif(isout.ge.2) then ! Angular differential calculation
  do ip=0,npart ! currently only neutron & proton
   if(iangpart(ip).ne.0) then
    iemin=61 ! basically down to 10 keV
    if(ip.eq.0) iemin=1 ! only neutron, down to 1e-8 MeV
    do ie=iemin,nebin+1
     do ia=1,nang
      FluxAngEdif(ip,ie,ia)=Flux(ip,ie)*getSpecAngFinal(iangpart(ip),s,r,d,emid(ie),g,ang(ia)) ! get flux in (/cm2/s/sr)
     enddo
    enddo
    write(100+ip,1015) s,r,d,g 
    if(isout.eq.2) then ! as a function of energy
     write(100+ip,'(202i13)') (ia,ia=1,nang+1)
     write(100+ip,'(''       E(MeV)'',201f13.2,'' (/(MeV/n)/cm2/s/sr)'')') (ang(ia),ia=1,nang)
     do ie=iemin,nebin ! output result
      write(100+ip,'(202es13.5)') emid(ie),(FluxAngEdif(ip,ie,ia),ia=1,nang)
     enddo
     if(iangpart(ip).eq.4) then ! for muon
      do ie=1,10 ! extra one order
       etmp=emid(ie)*1.0e14
       write(100+ip,'(202es13.5)') etmp,(getSpec(ip,s,r,d,etmp,g)*getSpecAngFinal(iangpart(ip),s,r,d,etmp,g,ang(ia)),ia=1,nang)
      enddo
     endif
    elseif(isout.eq.3) then ! as a function of cos(theta), relative value
     write(100+ip,'(200i13)') (ie,ie=1,nebin+3-iemin)
     write(100+ip,'(''   cos(theta)'',200es13.5,'' (/sr)'')') (emid(ie),ie=iemin,nebin+1)
     do ia=1,nang ! output result
      write(100+ip,'(200es13.5)') ang(ia),(FluxAngEdif(ip,ie,ia)/Flux(ip,ie),ie=iemin,nebin+1)
     enddo
    else ! as a function of cos(theta), absolute value
     write(100+ip,'(200i13)') (ie,ie=1,nebin+3-iemin)
     write(100+ip,'(''   cos(theta)'',200es13.5,'' (/(MeV/n)/cm2/s/sr)'')') (emid(ie),ie=iemin,nebin+1)
     do ia=1,nang ! output result
      write(100+ip,'(200es13.5)') ang(ia),(FluxAngEdif(ip,ie,ia),ie=iemin,nebin+1)
     enddo
    endif
   endif
  enddo
 endif
 if(itimedep.ne.1) then ! normal output
  write(4,1060) s,r,d,g,(dose(ip)*unitconv,ip=0,npart+1)
 else  ! output timedependence
  tmp=iyear+(imonth-1)/12.0+(iday-1)/365.0  ! rough year
  write(4,1065) s,r,d,g,(dose(ip)*unitconv,ip=0,npart+1),tmp,imonth,iday  
 endif
 do ip=0,npart+1  ! initialized
  dose(ip)=0.0
 enddo
 5 continue  
enddo
1010	format(a1,i4,a2,i3,a2,i4,a1,a1,i3)
1012	format(a1,i4,a2,i3,a2,i5,a1,a1,i3)
1015	format('W=',f6.1,'    , r=',f5.1,' (GV), d=',f6.1,' (g/cm^2), g=',f6.2)
1020	format('Total Dose Rate =',f15.7,' ',a6)
1025    format(100i13)
1030	format('Dose Rate   =',34es13.5,' ',a6)
1035	format('Dose Contri =',34es13.5,' %')
1040	format('Energy(MeV/n)',34(7x,a6),' (/(MeV/n)/cm^2/s)')
1050	format(100es13.5)
1060	format(f8.2,f8.2,f12.2,f8.2,35es13.5)
1065	format(f8.2,f8.2,f12.2,f8.2,35es13.5,f10.4,2i3)

999	continue

end

