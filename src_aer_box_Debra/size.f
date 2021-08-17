      SUBROUTINE SIZE
!
!           $Id: size.f,v 1.9 2012/12/16 23:48:53 dkweis Exp $
!
      include 'aerosol.i'
      real Rsphere(nszal),Raleff(nszal)
      COMMON /SIZE1/R(NSIZE),DR(NSIZE),VOL(NSIZE),SAD1(NSIZE)
      COMMON /SIZESOL/Rg(NSZAL),dRg(nszal),cores(NSZAL), &
     &    VOLal(NSZAL),SADAL1(NSZAL)
      COMMON /SIZEWET/Rgwet(NSZAL),dRgwet(nszal),SAD1wet(NSZAL)
      DATA PI/3.1415926/
      dimension dlogr(nsize)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCC                                                 CCCCC
!CCCC  THIS SUBROUTINE FILLS IN THE VOLUME, RADIUS    CCCCC
!CCCC  AND BIN WIDTH ARRAYS WITH THE APPROPRIATE      CCCCC
!CCCC  VALUES.  THE INITIAL RADIUS IS "RZERO" AND     CCCCC
!CCCC  BY CHANGING THAT, ONE CAN OBTAIN A DIFFERENT   CCCCC
!CCCC  SIZE SPREAD IN THE TWENTY BINS.                CCCCC
!CCCC                                                 CCCCC
!CCCC  NOTE THAT UNITS ARE MICRONS                    CCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!-----------------------------------------------------------
!          sizes for spherical H2SO4-H2O particles
      write(16,554) Rzero, Vrat, nsize
  554 format(/'H2SO4 SIZE PARAMETERS:  &
                        R0=',F7.5,'  Vrat=',F4.1,'  NSIZE=',i3)
      write(16,555)
555   FORMAT(4x,'BIN#',10X,'R(um)',9x,'R1(um)',9x,'R2(um)',9x, &
     &  'dR(um)',6x,'dlogR(um)',7x,'Vol(cm^3)',6x,'SAD(cm^-1)')
      VOL0=(4./3.)*PI*(Rzero**3)
      DO 20 K=1,NSIZE
         VOL(K)=VOL0*Vrat**(K-1)
         R(K)=(VOL(K)*3.0/(4.0*PI))**(1./3.)
         r1=r(k)*(2./(1.+Vrat))**(1./3.)
         r2=r(k)*(2.*Vrat/(1.+Vrat))**(1./3.)
         dr(k)=r2-r1
!##         dlogr(k)=log10(r2)-log10(r1)
         dlogr(k)=log(r2)-log(r1)
         VOL(K)=VOL(K)*1.E-12
         SAD1(K)=4.*PI*R(K)*R(K)*1.E-8
         write(16,666) K,R(K),r1,r2,DR(k),dlogr(k),VOL(K),SAD1(K)
666   FORMAT(4X,I4,5X,f10.7,5X,f10.7,5x,f10.7,5x,F10.7,5x,F10.7,5X,1PE11.5,5x,1PE11.5)
20    CONTINUE

!-----------------------------------------------------------
!          sizes for solid Al2O3 particles
      write(16,557) Rzal, Vrat, nszal,rhoal,Df,Dfh
  557 format(/'ALUMINA SIZE PARAMETERS:  R0=',F6.3,'  Vrat=',F4.1,  &
     & '  NSIZE=',i3,'  DENSITY=',f6.2,'  Df=',f4.1,'  Dh=',f4.1)
      write(16,556)
  556 FORMAT(4x,'BIN#',7X,'#CORES',9X,'Rg(um)',8x, &
     &  'R1g(um)',8x,'R2g(um)',7x,'Reff(um)',7X,'Rsph(um)',7x, &
     &  'Vol(cm^3)',6x,'SAD(cm^-1)')
      VOL0=(4./3.)*PI*(RzAl**3)
      DO 40 K=1,NSZAL
         cores(k) = Vrat**(k-1)  ! # of primary particles in aggregate K
         Rg(k) = RzAl * cores(k)**(1./Df)
         r1=Rg(k)*(2./(1.+Vrat))**(1./Df)
         r2=Rg(k)*(2.*Vrat/(1.+Vrat))**(1./Df)
         DRg(k)=r2-r1
         Raleff(k) = RzAl * cores(k)**(1./Dfh)
         SADal1(k)=4.*pi*Raleff(k)*Raleff(k)*1.e-8
         Volal(K)=VOL0*cores(k)
         Rsphere(K)=(Volal(K)*3.0/(4.0*PI))**(1./3.)
         volal(k)=volal(k)*1.e-12
         write(16,667) K,cores(k),Rg(K),r1,r2,Raleff(k),Rsphere(k),Volal(K),SADal1(k)
  667    FORMAT(4X,I4,5X,f8.0,5x,f10.6,5X,f10.6,5x,f10.6,5x,f10.6,5X,f10.6,5x,1pE11.5,5x,1pE11.5)
40    CONTINUE

!-----------------------------------------------------------
!          sizes for wetted Al2O3 particles
      write(16,558) Rzal, Vrat, nszal,rhoal,Dfwet,Dfhwet
  558 format(/'WETTED ALUMINA SIZE PARAMETERS:  R0=',F6.3,'  Vrat=',F4.1,  &
     & '  NSIZE=',i3,'  DENSITY=',f6.2,'  Dfwet=',f4.1,'  Dhwet=',f4.1)
      write(16,559)
  559 FORMAT(4x,'BIN#',7X,'#CORES',6X,'Rgwet(um)',8x, &
     &  'R1g(um)',8x,'R2g(um)',7x,'Reff(um)',7X,'Rsph(um)',7x, &
     &  'Vol(cm^3)',6x,'SAD(cm^-1)')
      DO K=1,NSZAL
         cores(k) = Vrat**(k-1)  ! # of primary particles in aggregate K
         Rgwet(k) = RzAl * cores(k)**(1./Dfwet)
         r1=Rg(k)*(2./(1.+Vrat))**(1./Dfwet)
         r2=Rg(k)*(2.*Vrat/(1.+Vrat))**(1./Dfwet)
         DRgwet(k)=r2-r1
         Ralef = RzAl * cores(k)**(1./Dfhwet)
         SAD1wet(k)=4.*pi*Ralef*Ralef*1.e-8
         write(16,667) K,cores(k),Rgwet(K),r1,r2,Ralef,Rsphere(k),Volal(K),SAD1wet(k)
      END DO

      RETURN
      END
