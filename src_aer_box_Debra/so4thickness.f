      Subroutine SO4THICKNESS(isz,alm,also4,wpp,dens,Rgm,sulfthick,Rcurve,totweight,sadone)
      include 'aerosol.i'
      COMMON /SIZESOL/Rg(NSZAL),dRg(nszal),cores(NSZAL), &
     &    VOLal(NSZAL),SADAL1(NSZAL)
      COMMON /SIZEWET/Rgwet(NSZAL),dRgwet(nszal),SAD1wet(NSZAL)
      COMMON/FLAGS/IFL(90)
      PARAMETER (pi=3.1415926)
      PARAMETER (onethird=1./3.)
!
!        WRITE PROGRAM VERSION NUMBER
      IF(IFL(58).EQ.0) THEN
         WRITE(16,*) '$Id: so4thickness.f,v 1.1 2014/05/14  dkweis Exp $'
         IFL(58)=1
      END IF
!      
!  This routine calculates sulfate thickness on a mixed alumina-sulfate particles
!           for one grid point and one size bin
!      totweight and sadone are weight and surface area density for a single particle
!
!##      write(16,*) 'start so4thick: ',isz,alm,also4,wpp,dens
!
!  Default values when alm = 0.
         Rgm=Rgwet(isz)*1.e-4
         sulfthick=0.
         Rcurve=Rgwet(1)*1.e-4
         totweight=volal(isz)*rhoal
         sadone=sad1wet(isz)
         if(alm.le.0.) return
!
! sulfmono is thickness of an H2SO4 monolayer in cm 
!    ie. H2SO4 molecule liquid diameter ~ 0.6 nm
         sulfmono=2.*(98./((4./3.)*PI*dens*wpp*6.02E23))**onethird
!
!  Calculations for alm > 0.
         so4perpart = also4/alm ! molec H2SO4 per part
         so4weight = so4perpart/6.02E23*98./wpp ! weight H2SO4+H2O per part
         so4vol = so4weight/dens ! volume H2SO4+H2O per part
         alumweight = rhoal*volal(isz) ! weight Al2O3 per part
         totweight = so4weight + alumweight ! total weight per part
         totvol = volal(isz) + so4vol ! particle volume in cm^3
         Raddry=Rgwet(isz)*1.e-4
         RSphere=(totvol*3./(4.*pi))**onethird
         Rcore = ((totvol/cores(isz))*3./(4.*pi))**onethird
         sulfthick = Rcore - Rgwet(1)*1.e-4
!!            sulfthick = max(sulfmono,sulfthick)
!!            sulfthick = min(sulfmono,sulfthick)
!!            sulfthick = min(10.*sulfmono,sulfthick)
         Rgm=Raddry + sulfthick
         Rcurve = Rgwet(1)*1.E-4 + sulfthick
         sadone=(4.*pi*Rcore**2) * (cores(isz))**(2./Dfhwet)
         if(RSphere.ge.Raddry+sulfmono) then
            Rgm = RSphere
            Rcurve=RSphere
            RAlSphere=(volal(isz)*3./(4.*pi))**onethird
            sulfthick=RSphere - RAlSphere
            sadone=4.*pi*RSphere**2
         end if
!##      write(16,*) 'end so4thickness: ',Rgm,sulfthick,Rcurve,totweight,sadone
         if(sulfthick.lt.0.) then
!            write(16,*) 'negative so4thickness = ',sulfthick,' for size ',isz,alm,also4
!            write(16,*) 'alm=',alm,'  also4=',also4,' totvol=',totvol,' Rcore=',Rcore
            sulfthick=0.
         end if
      return
      end
