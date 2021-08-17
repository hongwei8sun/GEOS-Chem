      SUBROUTINE NUCLEAR(bn,TSTEP)
!
!         $Id: nuclear.f,v 1.23 2012/12/15 14:43:50 dkweis Exp $
!
      include 'aerosol.i'
      COMMON /AIR/TMP,PR,AIRVD
      COMMON /CHEMA/H2O,H2SO4  ! number density 
      COMMON /SIZE1/RAD(NSIZE),DR(NSIZE),VOL(NSIZE),SAD1(NSIZE)
      COMMON /NUMBER1/ AN(NSIZE)
      COMMON /INFOR/WPCT,DEN,BVP,STS
      DIMENSION H2SO4MIN(7,9)
      DOUBLE PRECISION dtrate
      SAVE H2SO4MIN
      data h2so4min/ &
     & 1.00E+03,3.16E+02,1.00E+02,3.16E+00,1.00E-02,1.00E-02,1.00E-02, &
     & 1.00E+04,3.16E+03,3.16E+02,3.16E+01,3.16E+00,1.00E+00,1.00E+00, &
     & 3.16E+04,1.00E+04,3.16E+03,1.00E+03,3.16E+01,3.16E+00,3.16E+00, &
     & 3.16E+05,1.00E+05,3.16E+04,3.16E+03,3.16E+02,3.16E+01,1.00E+01, &
     & 1.00E+07,3.16E+05,1.00E+05,3.16E+04,3.16E+03,3.16E+02,3.16E+01, &
     & 1.00E+07,1.00E+07,1.00E+07,1.00E+05,3.16E+04,3.16E+03,1.00E+02, &
     & 1.00E+07,1.00E+07,1.00E+07,1.00E+07,1.00E+05,3.16E+04,1.00E+03, &
     & 1.00E+07,1.00E+07,1.00E+07,1.00E+07,1.00E+07,1.00E+05,1.00E+04, &
     & 1.00E+07,1.00E+07,1.00E+07,1.00E+07,1.00E+07,1.00E+07,1.00E+05/
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC
!CCC   This subroutine calculates the nucleation rate for a given
!CCC   ambient condtion of temperature temper, water vapor concentration
!CCC   water (no/cc) and  sulfuric acid vapor concentration H2SO4
!CCC   (number per cc)
!CCC
!CCC   radius is the radius of the embryo, rate is the nucleation rate
!CCC    com is the acidity of the embryo, drate is the rate of
!CCC    depletion of H2SO4
!CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!
      PARAMETER (AV=6.0225E23)
      COMMON /FLAGS/IFL(90)

      bn=0.
      dn=0.
      write(6,*) 'nuclear, H2SO4=',H2SO4

            dtrate=0.d0
            dtradius=0.
               ttmp=TMP
               if(ttmp.ge.232.) go to 600
               itmp=(TTMP-190.)/5. + 1.
               if(itmp.lt.1) itmp=1
               if(itmp.gt.9) itmp=9
               ih2o=(alog10(h2o)-12.)/0.5 + 1.95
               if(ih2o.lt.1) ih2o=1
               if(ih2o.gt.7) ih2o=7
               if(h2so4.le.h2so4min(ih2o,itmp)) go to 600
!
               call vehkamaki(TTMP,H2O,H2SO4,dtrate,dtradius,WPCT,DEN)
               if(dtrate.lt.0.d0) then
                  write(16,*) 'Nucleation trouble'
                  write(16,*) dtrate,dtradius,ttmp
                  STOP 'ERROR IN NUCLEAR'
               end if
               if(dtrate.le.0.d0) then
                  bn=0.
                  dn=0.
                  go to 600
               end if
!
            do ksize=1,NSIZE
               rupper=rad(ksize)*(2.*Vrat/(1.+Vrat))**(1./3.)
               if(dtradius*1.E+4.le.rupper) go to 550
            end do
            stop 'ERROR determining radius for nucleation'
  550       continue
!              bn is loss of H2SO4 gas, dn is addition of particles
            bn=dtrate*tstep
            dn=dtrate*tstep*98./(vol(ksize)*den*wpct*0.01*Av)
            if(bn.lt.h2so4) then
               h2so4=h2so4-bn
            else
               ratio=h2so4/bn
               dn=dn*ratio
               bn=h2so4
               h2so4=0.
            end if
            an(ksize)=an(ksize)+dn
  600    CONTINUE
  700 CONTINUE

!
      RETURN
      END
!
! *********************************************************************
!
      SUBROUTINE VEHKAMAKI(TEMP,WATER,H2SO4C,DRATE,RADIUS,WPBACK, &
     &     DENBACK)
!
      DOUBLE PRECISION drate
      REAL MW,MO
      PARAMETER (Mw=18.,Ma=98.)
      PARAMETER (PI=3.1415926,Av=6.0225E23)
!
!---->define limit values for nucleation rate and mass of nucleus
!
      parameter (rnucmax=1.E12)
      parameter (rmasmin=2.E0,rmasmax=300.E0)
      parameter (so4gmin=1.E1,so4gmax=1.E12)
      parameter (tmin=190.d0)
!
      drate=0.d0
      radius=0.
      if(h2so4c.lt.so4gmin) return
!         saturation number densities (cm**-3)
      h2so4sat=sath2so4(temp)
      h2osat=sath2o(temp)
!         RH=relative humidity, RA=relative acidity
      RH=min(1.,water/h2osat)
      RA=min(1.,h2so4c/h2so4sat)
!##      write(16,*) 'In Vehkamaki, T=',temp,'  h2so4=',h2so4c,
!##     $     '  h2o=',water,'  h2so4sat=',h2so4sat,
!##     $     '  h2osat=',h2osat,'  rh=',rh,'  ra=',ra
!
!---->calculate h2so4 mole/mass fraction in critical nucleus
!
      tempnew=max(tmin,temp)
      so4new=min(so4gmax,h2so4c)
!##      write(16,*) 'tempnew=',tempnew,'  so4new=',so4new,'  rh=',rh
!
!---->calculate h2so4 mole fraction of critical nucleus
      so4nuc=so4fracnuc_f(tempnew,rh,so4new)
      so4nuc=min(1.,max(0.,so4nuc))
      so4mnuc=so4nuc*Ma &
     &     /(so4nuc*Ma+(1.-so4nuc)*Mw)
!##      write(16,*) 'so4nuc=',so4nuc,'  so4mnuc=',so4mnuc
!
!---- >calculate nucleation rate: rnucnum [#particles/cm3]
      rnucnum=rnucnum_f(tempnew,rh,so4new,so4nuc)
!      rnucnum=min(rnucmax,rnucnum)
!##      write(16,*) 'rnucnum=',rnucnum
      if (rnucnum.le.0.) return
!
!---->calculate the number of total and h2so4 molecules in the critical cluster
         rnuctot=rnucmas_f(tempnew,rh,so4new,so4nuc)
         rnucso4=so4nuc*rnuctot   ! H2SO4 molecules per particle
!           H2SO4 molecules per second nucleated
         rnucso4=min(rmasmax,max(rmasmin,rnucso4))
         rnucmas=rnucnum*rnucso4
!##         write(16,*) 'rnuctot=',rnuctot,'  rnucso4=',rnucso4,
!##     $        '  rnucmas=',rnucmas
         if(rnucmas.lt.0.) then
            write(16,*) 'Trouble in Vehkamaeki'
            write(16,*) 'Inputs: T=',TEMP,'  H2O=',WATER,'  H2SO4=', &
     &           H2SO4C
            write(16,*) 'so4nuc=',so4nuc,so4mnuc,'  rnucnum=',rnucnum
            write(16,*) 'rnuctot=',rnuctot,rnucso4,'  rnucmas=',rnucmas
            return
         end if
!---->calculate the radius of the critical cluster (centimeters)
!
         rnucrad=1.E-7*exp(-1.6524245 + 0.42316402*so4nuc + &
     &        0.3346648*log(rnuctot))
         densnuc=density(tempnew,so4mnuc)
         volnuc=rnucso4/Av*98./so4mnuc/densnuc
         radnuc=(volnuc*3./4./pi)**0.3333
!    adjust particle radius for ambient density and wp
         radius=radnuc*(so4mnuc*densnuc/(wpback*0.01)/denback)**0.33333
!##         write(16,*) 'rnucrad=',rnucrad,'  volnuc=',volnuc,
!##     $        '  radnuc=',radnuc,'  densnuc=',
!##     $        densnuc,' radius=',radius,'  Av=',Av,'  pi=',pi,
!##     $        '  so4mnuc=',so4mnuc
      drate=rnucmas
!
      RETURN
      END
!
!-----------------------------------------------------------------------
      function rnucnum_f(temp,relhum,so4gas,so4frac)
!
!---->calculate nucleation rate: rnucnum [#particles/cm3]
!
      implicit real (a-h,o-z), integer (i-n)
!
!---->Vehkamaeki et al. (2002)
!
      parameter (a1=0.14309d0, &
     &           a2=2.21956d0, &
     &           a3=2.73911d-2, &
     &           a4=7.22811d-5, &
     &           a5=5.91822d0)
      parameter (b1=0.117489d0, &
     &           b2=0.462532d0, &
     &           b3=1.18059d-2, &
     &           b4=4.04196d-5, &
     &           b5=15.7963d0)
      parameter (c1=0.215554d0, &
     &           c2=8.10269d-2, &
     &           c3=1.43581d-3, &
     &           c4=4.7758d-6, &
     &           c5=2.91297d0)
      parameter (d1=3.58856d0, &
     &           d2=4.9508d-2, &
     &           d3=2.1382d-4, &
     &           d4=3.10801d-7, &
     &           d5=2.93333d-2)
      parameter (e1=1.14598d0, &
     &           e2=6.00796d-1, &
     &           e3=8.64245d-3, &
     &           e4=2.28947d-5, &
     &           e5=8.44985)
      parameter (f1=2.15855d0, &
     &           f2=8.08121d-2, &
     &           f3=4.07382d-4, &
     &           f4=4.01947d-7, &
     &           f5=7.21326d-1)
      parameter (g1=1.6241d0, &
     &           g2=1.60106d-2, &
     &           g3=3.77124d-5, &
     &           g4=3.21794d-8, &
     &           g5=1.13255d-2)
      parameter (h1=9.71682d0, &
     &           h2=1.15048d-1, &
     &           h3=1.57098d-4, &
     &           h4=4.00914d-7, &
     &           h5=0.71186d0)
      parameter (p1=1.05611d0, &
     &           p2=9.03378d-3, &
     &           p3=1.98417d-5, &
     &           p4=2.46048d-8, &
     &           p5=5.79087d-2)
      parameter (q1=0.148712d0, &
     &           q2=2.83508d-3, &
     &           q3=9.24619d-6, &
     &           q4=5.00427d-9, &
     &           q5=1.27081d-2)
      parameter (so4gmin=1.E04,so4gmax=1.E11)
      parameter (expmax=46.E0)

      sfracinv=1./so4frac
      temp2=temp*temp
      temp3=temp*temp2

      a=+a1+a2*temp-a3*temp2+a4*temp3+a5*sfracinv
      b=+b1+b2*temp-b3*temp2+b4*temp3+b5*sfracinv
      c=-c1-c2*temp+c3*temp2-c4*temp3-c5*sfracinv
      d=-d1+d2*temp-d3*temp2+d4*temp3-d5*sfracinv
      e=+e1-e2*temp+e3*temp2-e4*temp3-e5*sfracinv
      f=+f1+f2*temp-f3*temp2-f4*temp3+f5*sfracinv
      g=+g1-g2*temp+g3*temp2+g4*temp3-g5*sfracinv
      h=+h1-h2*temp+h3*temp2+h4*temp3+h5*sfracinv
      p=-p1+p2*temp-p3*temp2+p4*temp3-p5*sfracinv
      q=-q1+q2*temp-q3*temp2+q4*temp3-q5*sfracinv


      rhln=log(relhum)
      rhln2=rhln*rhln
      rhln3=rhln*rhln2
      so4ln=log(min(so4gmax,max(so4gmin,so4gas)))
      so4ln2=so4ln*so4ln
      so4ln3=so4ln*so4ln2

      expon=min(expmax,a+b*rhln+c*rhln2+d*rhln3 &
     &                +e*so4ln+f*rhln*so4ln+g*rhln2*so4ln &
     &                +h*so4ln2+p*rhln*so4ln2+q*so4ln3)
      rnucnum_f=exp(expon)

      return
      end
!-----------------------------------------------------------------------
      function so4fracnuc_f(temp,relhum,so4gas)
!
!---->calculate h2so4 molfraction of critical nucleus
!
      implicit real(a-h,o-z), integer (i-n)
!
!---->Vehkamaeki et al. (2002)
!
      parameter (a1=7.40997d-1, &
     &           a2=2.66379d-3, &
     &           a3=3.49998d-3, &
     &           a4=5.04022d-5, &
     &           a5=2.01048d-3, &
     &           a6=1.83289d-4, &
     &           a7=1.57407d-3, &
     &           a8=1.79059d-5, &
     &           a9=1.84403d-4, &
     &           a10=1.50345d-6)

      parameter (so4gmin=1.E04,so4gmax=1.E11)

      rhln=log(relhum)
      rhln2=rhln*rhln
      rhln3=rhln*rhln2
      so4ln=log(min(so4gmax,max(so4gmin,so4gas)))
      so4fracnuc_f=a1       -a2*temp &
     &            -a3*so4ln +a4*temp*so4ln &
     &            +a5*rhln  -a6*temp*rhln &
     &            +a7*rhln2 -a8*temp*rhln2 &
     &            +a9*rhln3-a10*temp*rhln3

      return
      end
!-----------------------------------------------------------------------
      function rnucmas_f(temp,relhum,so4gas,so4frac)
!
!---->calculate total number of molecules in critical cluster [#molec/part]
!     Vehkamaeki et al. (2002)
!
      implicit real (a-h,o-z), integer (i-n)

      parameter (a1=2.95413d-3, &
     &           a2=9.76834d-2, &
     &           a3=1.02485d-3, &
     &           a4=2.18646d-6, &
     &           a5=1.01717d-1)
      parameter (b1=2.05064d-3, &
     &           b2=7.58504d-3, &
     &           b3=1.92654d-4, &
     &           b4=6.70430d-7, &
     &           b5=2.55774d-1)
      parameter (c1=3.22308d-3, &
     &           c2=8.52637d-4, &
     &           c3=1.54757d-5, &
     &           c4=5.66661d-8, &
     &           c5=3.38444d-2)
      parameter (d1=4.74323d-2, &
     &           d2=6.25104d-4, &
     &           d3=2.65066d-6, &
     &           d4=3.67471d-9, &
     &           d5=2.67251d-4)
      parameter (e1=1.25211d-2, &
     &           e2=5.80655d-3, &
     &           e3=1.01674d-4, &
     &           e4=2.88195d-7, &
     &           e5=9.42243d-2)
      parameter (f1=3.85460d-2, &
     &           f2=6.72316d-4, &
     &           f3=2.60288d-6, &
     &           f4=1.19416d-8, &
     &           f5=8.51515d-3)
      parameter (g1=1.83749d-2, &
     &           g2=1.72072d-4, &
     &           g3=3.71766d-7, &
     &           g4=5.14875d-10, &
     &           g5=2.68660d-4)
      parameter (h1=6.19974d-2, &
     &           h2=9.06958d-4, &
     &           h3=9.11728d-7, &
     &           h4=5.36796d-9, &
     &           h5=7.74234d-3)
      parameter (p1=1.21827d-2, &
     &           p2=1.06650d-4, &
     &           p3=2.53460d-7, &
     &           p4=3.63519d-10, &
     &           p5=6.10065d-4)
      parameter (q1=3.20184d-4, &
     &           q2=1.74762d-5, &
     &           q3=6.06504d-8, &
     &           q4=1.42177d-11, &
     &           q5=1.35751d-4)


      parameter (so4gmin=1.E04,so4gmax=1.E11)

      sfracinv=1./so4frac
      temp2=temp*temp
      temp3=temp*temp2

      a=-a1-a2*temp+a3*temp2-a4*temp3-a5*sfracinv
      b=-b1-b2*temp+b3*temp2-b4*temp3-b5*sfracinv
      c=+c1+c2*temp-c3*temp2+c4*temp3+c5*sfracinv
      d=+d1-d2*temp+d3*temp2-d4*temp3-d5*sfracinv
      e=-e1+e2*temp-e3*temp2+e4*temp3+e5*sfracinv
      f=-f1-f2*temp+f3*temp2+f4*temp3-f5*sfracinv
      g=-g1+g2*temp-g3*temp2-g4*temp3+g5*sfracinv
      h=-h1+h2*temp-h3*temp2-h4*temp3-h5*sfracinv
      p=+p1-p2*temp+p3*temp2-p4*temp3+p5*sfracinv
      q=+q1-q2*temp+q3*temp2-q4*temp3+q5*sfracinv


      rhln=log(relhum)
      rhln2=rhln*rhln
      rhln3=rhln*rhln2
      so4ln=log(min(so4gmax,max(so4gmin,so4gas)))
      so4ln2=so4ln*so4ln
      so4ln3=so4ln*so4ln2

      rnucmas_f=exp(a+b*rhln+c*rhln2+d*rhln3 &
     &             +e*so4ln+f*rhln*so4ln+g*rhln2*so4ln &
     &             +h*so4ln2+p*rhln*so4ln2+q*so4ln3)

      return
      end
