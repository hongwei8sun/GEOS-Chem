      SUBROUTINE WPDEN(H2O)
!
!         $Id: wpden.f,v 1.22 2012/12/16 23:55:24 dkweis Exp $
!
      include 'aerosol.i'
      DIMENSION Rgm(NSZAL),Rsphere(nszal)
      DIMENSION VOLALM(NSZAL),DENSALM(NSZAL)
      CHARACTER*24 IDENT,AERID,ALUMID,ALMIXID,ALSO4ID
      COMMON /AIR/TMP,PR,AIRVD
      common /AEROID/AERID(nsize),ALUMID(NSZAL),almixid(nszal),also4id(nszal)
      COMMON /INFOR/WP,DEN,BVP,ST
      COMMON /NUMBER1/ AN(NSIZE)
      COMMON /SOLIDS/ALN(NSZAL),ALM(NSZAL),ALSO4(NSZAL)
      COMMON /SIZESOL/Rg(NSZAL),dRg(nszal),cores(NSZAL), &
     &    VOLal(NSZAL),SADAL1(NSZAL)
      COMMON /SIZEWET/Rgwet(NSZAL),dRgwet(nszal),SAD1wet(NSZAL)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC                                                             CCCCC
!CCC This subroutine calculates the density of sulfate particles CCCCC
!CCC                                                             CCCCC
!CCC Input parameter: WP Weight percentage of H2SO4 in the       CCCCC
!CCC                     sulfate aerosol particle                CCCCC
!CCC                  TMP Temperature in K                       CCCCC
!CCC                                                             CCCCC
!CCC Output parameter: Den density in gram per c.c.              CCCCC
!CCC                                                             CCCCC
!CCC DD0 are 101 values of density at 0 C for weight %=0 to 100  CCCCC
!CCC DD10 are 101 values of density at 10 C for weight %=0 to 100 CCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER (PI=3.1415926,AV=6.0225E23,AKB=1.38054E-16)
!
!
!
!               CALCULATE WEIGHT % H2SO4 IN AEROSOLS AS FN OF TMP
!               AND H2O PARTIAL PRESSURE USING PARAMETERIZATION OF
!               TABAZADEH, TOON, CLEGG, AND HAMILL, GRL, 24, 1931, 1997
!
!               partial pressure of ambient h2o
            pph2o=h2o*pr
      write(6,*) 'entered wpden with h2o=',h2o,'  pr=',pr, &
  &      '  h2opp=',pph2o
!               saturation water vapor partial pressure
            pph2osat=18.452406985-3505.1578807/tmp &
     &           -330918.55082/tmp/tmp &
     &           +12725068.262/(tmp**3)
            pph2osat=exp(pph2osat)
!               water activity
            aw=pph2o/pph2osat
            if(aw.gt.1.0) aw=1.0
            if(aw.le.0.05) then
               y1=12.372089320*aw**(-0.16125516114) &
     &              -30.490657554*aw -2.1133114241
               y2=13.455394705*aw**(-0.19213122550) &
     &              -34.285174607*aw -1.7620073078
            else if(aw.le.0.85) then
               y1=11.820654354*aw**(-0.20786404244) &
     &              -4.8073063730*aw -5.1727540348
               y2=12.891938068*aw**(-0.23233847708) &
     &              -6.4261237757*aw -4.9005471319
            else
               y1=-180.06541028*aw**(-0.38601102592) &
     &              -93.317846778*aw +273.88132245
               y2=-176.95814097*aw**(-0.36257048154) &
     &              -90.469744201*aw +267.45509988
            end if
            sulfmolal=y1+((tmp-190.)*(y2-y1)/70.)
            wp=9800.*sulfmolal/(98.*sulfmolal+1000.)
            if (wp.lt.15.) then
               wp=15.
            end if
            if (wp.gt.99.) wp=99.
            if (wp.lt.15. .or. wp.gt.100.) then
               write(16,*) 'weight percent= ',wp,' at ',i,j, &
     &              '  T=',tmp
               stop 'Stupid Value of Weight Percent'
            end if
!                mass fraction
            wpp=wp*0.01
!                mole fraction
            xa=18.*wpp/(18.*wpp+98.*(1.-wpp))
            write(6,*) '  wp=',wp,'  wpp=',wpp,'  xa=',xa
!
!               CALCULATE DENSITY OF AEROSOLS (GM/CC) AS FN OF WT %
!               AND TEMPERATURE
            den=density(tmp,wpp)
            write(6,*) '  den=',den
!
!               CALCULATE SURFACE TENSION OF AEROSOLS AGAINST AIR
            st=surftension(tmp,xa)
            write(6,*) '  st=',st
!
!               CALCULATE EQUILIBRIUM VAPOR PRESSURE OF H2SO4
            bvp=solh2so4(tmp,xa)
            write(6,*) '  bvp=',bvp
!
!               CALCLATE DENSITY OF MIXED AL2O3-H2SO4 PARTICLES
!                  AL2O3 density is RHOAL, H2SO4 density is DEN
!                  grid box contains ALM mixed particles/cc and
!                     ALSO4 molecules of H2SO4/cc coating alumina
            do isz=1,nszal
               if(alm(isz).gt.1.e-20) then
                  so4perpart = also4(isz)/alm(isz) ! molec H2SO4 per part
                  so4weight = so4perpart/6.02E23*98./wpp ! weight H2SO4+H2O per part
                  so4vol = so4weight/den ! volume H2SO4+H2O per part
                  alumweight = rhoal*volal(isz) ! weight Al2O3 per part
                  totweight = so4weight + alumweight ! total weight per part
                  totvol = volal(isz) + so4vol ! particle volume in cm^3
                  volalm(isz) = totvol
                  densalm(isz) = totweight/totvol ! density in gm/cc
                  Rsphere(isz)=1.e4*(totvol*3./(4.*pi))**0.333
                  Rgm(isz)=Rgwet(isz) + so4vol/sad1wet(isz)*1.e4
                  Rgm(isz)=max(Rgm(isz),Rsphere(isz))
               else
                  densalm(isz)=rhoal
                  volalm(isz)=volal(isz)
                  Rgm(isz)=Rgwet(isz)
               end if
            end do
      RETURN
      END
!
!*************************************************************
!
      function density(temp,so4mfrac)
!        calculation of particle density
!        requires Temperature (temp) and acid mass fraction (so4mfrac)
!
!
!---->Vehkamaeki et al., 2002  (JGR, doi:10.1029/2002JD002184)

!
      parameter (a1= 0.7681724d0, &
     &           a2= 2.184714d0, &
     &           a3= 7.163002d0, &
     &           a4=-44.31447d0, &
     &           a5= 88.74606d0, &
     &           a6=-75.73729d0, &
     &           a7= 23.43228d0)
      parameter (b1= 1.808225d-3, &
     &           b2=-9.294656d-3, &
     &           b3=-3.742148d-2, &
     &           b4= 2.565321d-1, &
     &           b5=-5.362872d-1, &
     &           b6= 4.857736d-1, &
     &           b7=-1.629592d-1)
      parameter (c1=-3.478524d-6, &
     &           c2= 1.335867d-5, &
     &           c3= 5.195706d-5, &
     &           c4=-3.717636d-4, &
     &           c5= 7.990811d-4, &
     &           c6=-7.458060d-4, &
     &           c7= 2.581390d-4)

      so4m2=so4mfrac*so4mfrac
      so4m3=so4mfrac*so4m2
      so4m4=so4mfrac*so4m3
      so4m5=so4mfrac*so4m4
      so4m6=so4mfrac*so4m5

      a=+a1+a2*so4mfrac+a3*so4m2+a4*so4m3 &
     &        +a5*so4m4+a6*so4m5+a7*so4m6
      b=+b1+b2*so4mfrac+b3*so4m2+b4*so4m3 &
     &        +b5*so4m4+b6*so4m5+b7*so4m6
      c=+c1+c2*so4mfrac+c3*so4m2+c4*so4m3 &
     &        +c5*so4m4+c6*so4m5+c7*so4m6
      density=(a+b*temp+c*temp*temp) ! units are gm/cm**3
      return
      end
!
!************************************************************
!
      function surftension(temp,so4frac)
!        calculation of surface tension
!        requires Temperature (temp) and acid mole fraction (so4frac)
!        from Vehkamaeki et al., 2002 (JGR, doi:10.1029/2002JD002184)
!
!
!---->Vehkamaeki et al. (2002)
!
      parameter (a1= 0.11864d0, &
     &           a2=-0.11651d0, &
     &           a3= 0.76852d0, &
     &           a4=-2.40909d0, &
     &           a5= 2.95434d0, &
     &           a6=-1.25852d0)
      parameter (b1=-1.5709d-4, &
     &           b2= 4.0102d-4, &
     &           b3=-2.3995d-3, &
     &           b4= 7.611235d-3, &
     &           b5=-9.37386d-3, &
     &           b6= 3.89722d-3)
      parameter (convfac=1.d3)  ! convert from newton/m to dyne/cm

!        so4 mass fraction
      so4mfrac=98.*so4frac/(98.*so4frac+18.*(1-so4frac))
      so4m2=so4mfrac*so4mfrac
      so4m3=so4mfrac*so4m2
      so4m4=so4mfrac*so4m3
      so4m5=so4mfrac*so4m4

      a=+a1+a2*so4mfrac+a3*so4m2+a4*so4m3+a5*so4m4+a6*so4m5
      b=+b1+b2*so4mfrac+b3*so4m2+b4*so4m3+b5*so4m4+b6*so4m5
      so4sig=a+b*temp
      surftension=so4sig*convfac
      return
      end
!
!**********************************************************************
!
      function sath2so4(temp)
!
!---->Ayers et.al. (1980), GRL (7) pp 433-436
!     plus corrections for lower temperatures by Kulmala and Laaksonen (1990)
!     and Noppel et al. (1990)
!
      parameter (b1=1.01325d5, &
     &           b2=11.695d0, &
     &           b3=1.0156d4, &
     &           b4=0.38d0/545.d0, &
     &           tref=360.15d0)
      parameter (akb=1.38054d-16)

!         saturation vapor pressure (N/m**2)
      ppacid=b1*exp(-b2+b3*(1.d0/tref-1.d0/temp &
     &                        +b4*(1.d0+log(tref/temp)-tref/temp)))
!         saturation number density (cm**-3)
      sath2so4=ppacid*10.d0/(akb*temp)
      return
      end
!
!**********************************************************************
!
      function sath2o(T)
      parameter (akb=1.38054d-16)
      if(T.gt.229.) then
!         Preining et al., 1981 (from Kulmala et al., 1998)
!         saturation vapor pressure (N/m**2)
         ppwater=exp(77.34491296d0 -7235.424651d0/T -8.2d0*log(T) &
     &        +5.7133d-3*T)
!         saturation number densities (cm**-3)
         sath2o=10.d0/(akb*T)*ppwater
      else
!         Tabazadeh et al., 1997, parameterization for 185<T<260
!               saturation water vapor partial pressure (mb)
         pph2osat=18.452406985d0-3505.1578807d0/T &
     &        -330918.55082d0/T/T &
     &        +12725068.262d0/(T**3)
         pph2osat=exp(pph2osat)
         sath2o=1.d3/(akb*T)*pph2osat
      end if

      return
      end
!
!**********************************************************************
!
      function solh2o(T,xa)
      xw=1.0-xa
!           compute activity of water
      a11=2.989E3 -2.147E6/T + 2.33E8/(T*T)
      b11=0.527
      cactw=10.**(a11*xa*xa/(xa+b11*xw)**2/T)
      solh2o=cactw*xw*sath2o(T)
      return
      end
!
!**********************************************************************
!
      function solh2so4(T,xa)
      xw=1.0-xa
!          compute acitvity of acid
      a12=5.672E3 -4.074E6/T +4.421E8/(T*T)
      b12=1./0.527
      cacta=10.**(a12*xw*xw/(xw+b12*xa)**2/T)
      solh2so4=cacta*xa*sath2so4(T)
      return
      end
