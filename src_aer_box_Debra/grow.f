     SUBROUTINE GROW(dgas,tstep)
!
!         $Id: grow.f,v 1.16 2012/12/15 14:29:26 dkweis Exp $
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC                                                              CCC
!CCC This subroutine calculate the growth of sulfate particles    CCC
!CCC through hetermolecular condensation process                  CCC
!CCC                                                              CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      include 'aerosol.i'
      DIMENSION DN(NSIZE),DC(NSIZE),iout(nsize),DALSO4(NSZAL)
      character*24 ident
      SAVE IOUT
      COMMON /SIZE1/R(NSIZE),DR(NSIZE),VOL(NSIZE),SAD1(NSIZE)
      COMMON /SIZESOL/Rg(NSZAL),dRg(nszal),cores(NSZAL), &
     &    VOLal(NSZAL),SADAL1(NSZAL)
      COMMON /SIZEWET/Rgwet(NSZAL),dRgwet(nszal),SAD1wet(NSZAL)
      COMMON /NUMBER1/ AN(NSIZE)
      COMMON /SOLIDS/ALN(NSZAL),ALM(NSZAL), &
     &  ALSO4(NSZAL)
      COMMON /CHEMA/H2O,H2SO4  ! number density 
      COMMON /INFOR/WP,DEN,BVP,ST
      COMMON /AIR/TMP,PR,AIRVD
      COMMON /FLAGS/IFL(90)
      PARAMETER (PI=3.1415926, AV=6.0225E23)
      PARAMETER (AKB=1.38054E-16, Rgas=8.314E7)
      parameter (onethird=1./3.)
      data iout/nsize*1/
!
      if(ifl(55).eq.0) then
!         write(16,*) '$Id: grow.f,v 1.16 2012/12/15 14:29:26 dkweis Exp $'
         ifl(55)=1
      end if
      write(6,*) 'Grow: H2SO4=',H2SO4
!
      DEN1=0.2*1.14+0.8*0.808
!  AIR IS DIAMETER OF AN AIR MOLECULE, SULFATE IS DIAMETER OF AN H2SO4
      AIR=2.*(28.9/((4./3.)*PI*DEN1*AV))**onethird
      SULFATE=2.*(98./((4./3.)*PI*1.8409*AV))**onethird
      D=0.5*(AIR+SULFATE)
!  HAMU IS AMU TO THE POWER OF 0.5
      HAMU=(28.9/(28.9+98.))**0.5
      SULF=98./av
!
         wpp=wp*0.01
!      thickness of an H2SO4 monolayer in cm (ie. H2SO4 molecule liquid diameter ~ 0.6 nm)
!         sulfmono=2.*(98./((4./3.)*PI*den*wpp*6.02E23))**onethird

         DGAS=0.
         DN(:)=0.
         DALSO4(:)=0.
         do kk=1,nszal
            if(also4(kk).ge.0. .or. also4(kk).le.0.) then
               continue
            else
               write(16,*) 'also4=NaN start of grow at ',kk
               write(16,*) 'dalso4=',dalso4(kk),'  also4=',also4(kk)
               stop 'error start of grow'
            end if
         end do
!
      AD=AIRVD
      EFF=1./(PI*AD*D**2)*HAMU
!   WPP IS WEIGHT % OF H2SO4, FRA IS MOLE FRACTION
!      WPP=wp*0.01  ! calculated above
      FRA=18.*WPP/(18.*WPP+98.*(1.-WPP))
      avgmolwt=fra*98.+(1.-fra)*18.
      VOLUME=1./DEN
      AVE=VOLUME/(((WPP/98.)+(1.-WPP)/18.)*AV)
!   DIFF IS THE DIFFUSION COEFFICIENT
      DIFF=3./(8.*D**2*AD)*SQRT(AKB*TMP/(2.*PI*SULF)* &
     &     ((28.9+98.)/28.9)) ! shw

      write(*,*)'shw, DIFF:', DIFF ! shw
!
      DO 100 K=1,NSIZE
         RAD=R(K)*1.E-4
         VKNUD=EFF/RAD
         CORR=(1.3333+0.71/VKNUD)/(1.+1./VKNUD)
!     CURVATURE CORRECTION (KELVIN EFFECT)
         curvature=exp(2.*avgmolwt*st/(rad*Rgas*tmp*den))
         PREDIFF=H2SO4-BVP*curvature
!     CHANGE IN PARTICLE VOLUME PER SECOND PER CM**2 OF SURFACE AREA
         GROWTH=AVE*DIFF*PREDIFF/(RAD*FRA*(1.+CORR*VKNUD))

         write(6,*)'shw check GROWTH:', GROWTH
         write(6,*)'bin, DIFF, PREDIFF, (1.+CORR*VKNUD), AVE/(RAD*FRA):'
         write(6,*)K, DIFF, PREDIFF, (1.+CORR*VKNUD), AVE/(RAD*FRA)

!
!     The depletion of H2SO4 is given by
!
!$$$         DC(K)=-4.*PI*RAD*DIFF*PREDIFF/(1.+CORR*VKNUD)
!
!     Change in number of particles in bin k corresponding to volume change = S/V*growth*AN*tstep
         GROWTH=GROWTH*TSTEP
         CHANGE=3.*(GROWTH/RAD)*AN(K)
         IF (CHANGE.LT.0.) THEN
            IF (K.GT.1) then
               IF (Vrat/(Vrat-1.)*abs(CHANGE).GT.AN(K)) &
     &              CHANGE=-0.95*AN(K)*(Vrat-1.)/Vrat
               DN(K)=DN(K)+CHANGE*Vrat/(Vrat-1.)
               DN(K-1)=DN(K-1)-CHANGE*Vrat/(Vrat-1.)
            else
               if(abs(change).gt.AN(K)) change=-0.95*AN(K)
               dn(k)=dn(k)+change
            end if
         ELSE IF (CHANGE.GT.0.) THEN
            IF (K.LT.NSIZE) THEN
               IF (CHANGE/(Vrat-1.).GT.AN(K)) &
     &              CHANGE=0.95*AN(K)*(Vrat-1.)
               DN(K)=DN(K)-CHANGE/(Vrat-1.)
               DN(K+1)=DN(K+1)+CHANGE/(Vrat-1.)
            else
               dn(k)=dn(k)+change
            END IF
         END IF
  100 continue
      DO 160 K=1,NSIZE
         vmol=vol(k)*den*wp*0.01/98.*Av
         dgas=dgas-dn(k)*vmol
  160 CONTINUE

!      if(lalum) then
!         if(.not.nomixed) then
!  condensation to mixed alumina-H2SO4 particles
      DO 200 K=1,NSZAL
         if(alm(k).gt.0. .and. also4(k).gt.0.) then
            call so4thickness(k,alm(k),also4(k),wpp,den, &
     &                    Rad,sthick,Rcurve,partmass,sadone)
            VKNUD=EFF/RAD
            CORR=(1.3333+0.71/VKNUD)/(1.+1./VKNUD)
!     CURVATURE CORRECTION (KELVIN EFFECT)
            curvature=exp(2.*avgmolwt*st/(Rcurve*Rgas*tmp*den))
            PREDIFF=H2SO4-BVP*curvature
!     CHANGE IN H2SO4 on mixed particles, molecules/sec
            RATE=SADONE*DIFF*PREDIFF/(RAD*(1.+CORR*VKNUD))
         else
            RATE=0.
         end if
!
         RATE=RATE*ALM(K)*TSTEP
         IF (RATE.LT.0.) THEN  ! EVAPORATION
            IF (ABS(RATE).GT.ALSO4(K)) THEN  ! ALL H2SO4 EVAPORATED
               DGAS=dgas+ALSO4(K)
               DALSO4(K)=DALSO4(K)-ALSO4(K)
!                 MOVE MIXED PARTICLES TO PURE ALUMINA PARTICLES
            ELSE
               DALSO4(K)=DALSO4(K) + RATE
               DGAS=dgas-RATE
            END IF
         ELSE IF (RATE.GT.0.) THEN  ! CONDENSATION
            DALSO4(K)=DALSO4(K) + RATE
            DGAS=dgas-RATE
         END IF
  200 continue
!      end if  ! if(.not.nomixed)
!      end if  ! if(lalum)
!
      ratio=1.
      if(dgas.lt.0.) then
!           H2SO4 should not go below BVP, since sign of DN would change
         if( (h2so4+dgas) .lt. bvp) then
            ratio=(h2so4-bvp)/abs(dgas)
!            write(16,*) 'ratio=',ratio,h2so4,dgas,bvp
!            call flush(16)
         end if
      end if
      if(ratio.lt.0. .or. ratio.gt.1.0) then
         write(16,*) 'bad ratio in grow at point ',' ratio=',ratio
         write(16,*) 'H2SO4=',H2SO4,'  dgas=',dgas,'  bvp=',bvp
         write(16,*) 'h2so4-bvp=',h2so4-bvp,'  dgas=',dgas
         stop 'bad ratio in grow'
      end if
      dgas=dgas*ratio
      H2SO4=H2SO4+dgas
      if(h2so4.ge.0.) then
         continue
      else if (h2so4.le.0.) then
         continue
      else
         write(16,*) 'h2so4=nan in grow for point ',i,j
         write(16,*) 'dgas=',dgas,'  ratio=',ratio,'  bvp=',bvp
         stop 'error in grow'
      end if
      if(h2so4.lt.0.) then
         h2so4=0.
      end if
      do kk=1,nsize
         AN(KK)=AN(KK)+DN(KK)*ratio
         dn(kk)=dn(kk)*ratio/tstep
      end do
      totdalso4=0.
!      if(lalum) then
      do kk=1,nszal
         ALSO4(KK)=ALSO4(KK)+DALSO4(KK)*ratio
         DALSO4(kk)=DALSO4(kk)*ratio/tstep
         totdalso4=totdalso4+dalso4(kk)
! MOVE MIXED PARTICLE TO PURE AL if fractal dimension of wet and dry are same
         IF(ALSO4(KK).LE.0.) THEN
!            write(16,*) 'Grow, also4<0 at ',kk,'  alm=',alm(kk),aln(kk)
            ALSO4(KK)=0.
            IF(ALM(KK).GT.0. .AND. abs(Df-Dfwet).lt..001) THEN  
               ALN(KK)=ALN(KK) + ALM(KK)
               ALM(KK)=0.
!               write(16,*) 'New alm=',alm(kk),'  new aln=',aln(kk)
            END IF
!            if(also4(kk).le.0. .and. alm(kk).gt.0.) then
!               write(16,*) 'Dry aged particles for ',kk,'  alm=',alm(kk)
!            end if
         end if
         if(also4(kk).ge.0. .or. also4(kk).le.0.) then
            continue
         else
            write(16,*) 'also4=NaN in grow at ',kk
            write(16,*) 'dalso4=',dalso4(kk),'  also4=',also4(kk)
            stop 'error in grow'
         end if
      end do
!##      write(16,*) 'grow, point ','  totdalso4=',totdalso4
!      end if  ! if(lalum)
  300 continue

      RETURN
      END
