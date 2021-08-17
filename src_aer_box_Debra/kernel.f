!*********************************************************************
!
      SUBROUTINE COAGKERNEL
!
!*********************************************************************
!
      include 'aerosol.i'
      COMMON /SIZE1/R(NSIZE),DR(NSIZE),VOL(NSIZE),SAD1(NSIZE)
      COMMON /SIZESOL/Rg(NSZAL),dRg(nszal),cores(NSZAL), &
     &    VOLal(NSZAL),SADAL1(NSZAL)
      COMMON /SIZEWET/Rgwet(NSZAL),dRgwet(nszal),SAD1wet(NSZAL)
      COMMON /CaCO3/Rgca(nszal),Volca(nszal),SADca1(nszal),rhoca
      COMMON /KERNEL/AKSS(NSIZE,NSIZE),AKALAL(NSZAL,NSZAL), &
     &  AKSAL(NSIZE,NSZAL),AKSM(NSIZE,NSZAL),  &
     &  AKALM(NSZAL,NSZAL),AKMM(NSZAL,NSZAL)
      COMMON /INFOR/WP,DEN,BVP,ST
      COMMON /AIR/TMP,PR,AIRVD
      DIMENSION RGmix(nszal)
      DIMENSION D(NSIZE),G(NSIZE),DEL(NSIZE),AR(NSIZE)
      DIMENSION DAL(NSZAL),GAL(NSZAL),DELAL(NSZAL)
      DIMENSION DAM(NSZAL),GAM(NSZAL),DELAM(NSZAL)
      REAL LB,rad(nszal),XKNAL(NSZAL),cslip(nszal),jcores,kcores
      PARAMETER (PI=3.1415926, AKB=1.38054E-16, AV=6.0225E23)
      PARAMETER (FP0=6.6E-6, P0=1013., T0=293.15)
      parameter (onethird=1./3.)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                          cccccccc
!  This subroutine calculates the change of N due to       cccccccc
!      coagulation                                         cccccccc
!                                                          cccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      write(6,*) 'coagkernel'
      AR = R*1.E-4  ! CONVERT RADIUS FROM MICRONS TO CM

         wpp=wp*0.01
         write(6,*) 'coagkernel,  wpp=',wpp
         call flush(16)
!       thickness of an H2SO4 monolayer in cm (ie. H2SO4 molecule liquid diameter ~ 0.6 nm)
         sulfmono=2.*(98./((4./3.)*PI*den*wpp*6.02E23))**onethird
         TTE=TMP
!              TMP IN DEGREES CELCIUS
         TPC=TTE-273.15
         write(6,*) 'coagkernel,  sulfmono=',sulfmono,'  tte=',tte
!         call flush(16)
!              VISCOSITY
!##         IF(TPC.GE.0.) THEN
!##            VIS=(1.718+0.0049*TPC)*1.E-4
!##         ELSE
!##            VIS=(1.718+0.0049*TPC-1.2E-5*TPC*TPC)*1.E-4
!##         END IF
!        vis [kg/m s]-> [10 g/cm s]
         vis=1.8325E-4*(416.16/(TMP+120.))*(TMP/296.16)**1.5

!              MEAN FREE PATH OF AIR MOLECULES
         FP=FP0*(P0/PR)*(TTE/T0)
!
!  COAGULATION KERNEL FOR SULFATE - SULFATE COLLISIONS
!
         DO J=1,NSIZE
            AMASS=VOL(J)*DEN
            XKN=FP/AR(J)
            B=(1./(6.*PI*VIS*AR(J)))* &
     &           (1.+1.249*XKN+.42*XKN*EXP(-.87/XKN)) ! slip-flow (Kasten 1968)
!!     &           (1.+1.257*XKN+0.4*XKN*EXP(-1.1/XKN)) ! (Seinfeld and Pandis 2006)
            D(J)=AKB*TTE*B
            G(J)=SQRT(8.*AKB*TTE/(PI*AMASS))
            LB=(2./pi)*D(J)/G(J)  ! (Jacobson 1999)
!!            LB=(8./pi)*D(J)/G(J)  ! (Seinfeld and Pandis 2006)
            AAA=(2.*AR(J)+LB)**3.
            BBB=(4.*AR(J)*AR(J)+LB*LB)**(3./2.)
            DEL(J)=(1./(6.*AR(J)*LB))*(AAA-BBB)-2.*AR(J)
         END DO
         DO J=1,NSIZE
            DO K=1,NSIZE
               RR=AR(J)+AR(K)
               DD=D(J)+D(K)
               GG=SQRT(G(J)*G(J)+G(K)*G(K))
               DDEL=SQRT(DEL(J)*DEL(J)+DEL(K)*DEL(K))
               FAC1=RR/(RR+DDEL)
               FAC2=4.*DD/(GG*RR)
               FAC3=1./(FAC1+FAC2)
               AKSS(J,K)=4.*PI*RR*DD*FAC3
            END DO
         END DO
!
!  COAGULATION KERNEL FOR AL2O3 - AL2O3 COLLISIONS
!
         DO J=1,NSZAL  ! for alumina
            AMASS=VOLAL(J)*RHOAL
            RGcm=RG(J)*1.E-4
            XKN=FP/RGcm
            B=(1./(6.*PI*VIS*RGcm))* &
     &           (1.+1.249*XKN+.42*XKN*EXP(-.87/XKN)) ! slip-flow  (Kasten 1968)
!!     &           (1.+1.257*XKN+0.4*XKN*EXP(-1.1/XKN)) ! (Seinfeld and Pandis 2006)
            DAL(J)=AKB*TTE*B
            GAL(J)=SQRT(8.*AKB*TTE/(PI*AMASS))
            LB=(2./pi)*DAL(J)/GAL(J)
!!            LB=(8./pi)*D(J)/G(J)  ! (Seinfeld and Pandis 2006)
            AAA=(2.*RGcm+LB)**3.
            BBB=(4.*RGcm*RGcm+LB*LB)**(3./2.)
            DELAL(J)=(1./(6.*RGcm*LB))*(AAA-BBB)-2.*RGcm
         END DO
         DO J=1,NSZAL   ! ALUMINA
            DO K=1,NSZAL  ! alumina
               RR=(RG(J)+RG(K))*1.E-4
               DD=DAL(J)+DAL(K)
               GG=SQRT(GAL(J)*GAL(J)+GAL(K)*GAL(K))
               DDEL=SQRT(DELAL(J)*DELAL(J)+DELAL(K)*DELAL(K))
               FAC1=RR/(RR+DDEL)
               FAC2=4.*DD/(GG*RR)
               FAC3=1./(FAC1+FAC2)
               AKALAL(J,K)=4.*PI*RR*DD*FAC3
            END DO
         END DO
!
!  COAGULATION KERNEL FOR SULFATE - ALUMINA COLLISIONS
!
         DO J=1,NSIZE   ! sulfate
            DO K=1,NSZAL  ! alumina
               RR=AR(J)+RG(K)*1.E-4
               DD=D(J)+DAL(K)
               GG=SQRT(G(J)*G(J)+GAL(K)*GAL(K))
               DDEL=SQRT(DEL(J)*DEL(J)+DELAL(K)*DELAL(K))
               FAC1=RR/(RR+DDEL)
               FAC2=4.*DD/(GG*RR)
               FAC3=1./(FAC1+FAC2)
               AKSAL(J,K)=4.*PI*RR*DD*FAC3
            END DO
         END DO
!

! optional output
      write(6,*)
         WRITE(6,*) 'SULFATE COAG KERNELS'
         write(6,410) (k,k=1,nsize)
  410    format(3x,40i8)
         do j=1,nsize
            WRITE(6,420) J,(AKSS(J,K),k=1,nsize)
  420       format(i3,40(1pe9.2))
         END DO

      write(6,*)
         WRITE(6,*) 'AL2O3-AL2O3 COAG KERNELS'
         write(6,410) (k,k=1,nszal)
         do j=1,nszal
               WRITE(6,420) J,(AKALAL(J,K),k=1,nszal)
         END DO


      write(6,*)
         WRITE(6,*) 'SULFATE-AL2O3 COAG KERNELS'
         write(6,410) (k,k=1,nszal)
         do j=1,nsize
               WRITE(6,420) J,(AKSAL(J,K),k=1,nszal)
         END DO



      return
      end
