      include 'aerosol.i'
      COMMON /AIR/TMP,PR,AIRVD
      COMMON /INFOR/WP,DEN,BVP,ST
      COMMON /CHEMA/H2Ond,H2SO4nd
      COMMON /NUMBER1/ AN(NSIZE)
      COMMON /SOLIDS/ALN(NSZAL),ALM(NSZAL),ALSO4(NSZAL)
      character*80 arg

!      data H2SO4/0./
!      data H2O/10.E-6/
!      DATA TMP/213/
!      DATA PR/54.7/
!      data H2SO4/10.E-9/
!      data H2O/2.e-6/
!      data H2O/50.E-9/
!      data TMP/240./
! SPARC Fig 6.2a
      data TMP/200./            ! [K]
      data PR/124./             ! [mb]
      data H2O/3.14E-4/         ! [mol/mol]
      DATA H2SO4/40.E-12/       ! [mol/mol]
!  H2SO4 injection rate of 40 pptv per day
      H2SO4rate = 0.0 ! [(mol/mol)/s]
!        AN is sulfate particle concentration in #/cc
      AN(:) = 1000.0
      write(6,*) 'Initial AN:  ',AN
!      AN(1) = 1.E5
!       ALN, ALM, ALSO4 are solid, mixed particles and SO4 on mixed particles
      ALN(:)=0.
!      ALN(1)=1.E5
      ALM(:)=0.
      ALSO4(:)=0.

!      tstep = 0.01
!      time = 20.
!      ncoag=1
!      ngrow=1
      tstep=3.          ! [s]
      time= 3600*24     ! [s] total run time
      ncoag=1
      ngrow=1

      nargs=iargc()
      if(nargs.lt.1) go to 2
      iarg=1
    1 continue
         CALL GETARG(iarg,arg)
         if(arg(1:1).eq.'-') then
            if(arg(2:4).eq.'h2o' .or. arg(2:4).eq.'H2O') then
               iarg=iarg+1
               call getarg(iarg,arg)
               H2O=ffroma(arg,80)
            else if(arg(2:6).eq.'h2so4' .or. arg(2:4).eq.'H2SO4') then
               iarg=iarg+1
               call getarg(iarg,arg)
               H2SO4=ffroma(arg,80)
            else if(arg(2:4).eq.'tmp' .or. arg(2:4).eq.'TMP') then
               iarg=iarg+1
               call getarg(iarg,arg)
               TMP=ffroma(arg,80)
            else if(arg(2:3).eq.'pr' .or. arg(2:3).eq.'PR') then
               iarg=iarg+1
               call getarg(iarg,arg)
               pr=ffroma(arg,80)
            else if(arg(2:6).eq.'ncoag' .or. arg(2:6).eq.'NCOAG') then
               iarg=iarg+1
               call getarg(iarg,arg)
               ncoag=ifroma(arg,80)
            else if(arg(2:6).eq.'ngrow' .or. arg(2:6).eq.'NGROW') then
               iarg=iarg+1
               call getarg(iarg,arg)
               ngrow=ifroma(arg,80)
            else if(arg(2:3).eq.'dt' .or. arg(2:3).eq.'DT') then
               iarg=iarg+1
               call getarg(iarg,arg)
               tstep=ffroma(arg,80)
            else if(arg(2:5).eq.'time' .or. arg(2:5).eq.'TIME') then
               iarg=iarg+1
               call getarg(iarg,arg)
               time=ffroma(arg,80)
            end if
         endif
         iarg=iarg+1
         if(iarg.le.nargs) go to 1
    2 continue
      open(16,file='boxmod.out')
      write(16,*) 'Aerosol Box Model'
      write(16,*)
      write(6,*) 'Input Parameters:  T=',tmp,' Pr=',pr,'  H2O=',H2O,'  H2SO4=',H2SO4
      write(16,*) 'Input Parameters:  T=',tmp,' Pr=',pr,'  H2O=',H2O,'  H2SO4=',H2SO4
      write(6,*) '     ncoag=',ncoag,'  ngrow=',ngrow,'  tstep=',tstep,'  time=',time
      write(16,*) '    ncoag=',ncoag,'  ngrow=',ngrow,'  tstep=',tstep,'  time=',time
      write(6,*)
      write(16,*)

      AIRVD = PR/(1.38E-19*TMP)
      H2Ond = H2O*AIRVD
      H2SO4nd = H2SO4*airvd
      H2SO40=H2SO4nd

      call size

      call wpden(h2o)
      write(6,*) 'wp=',wp,'  den=',den,'  bvp=',bvp,'  st=',st
      write(6,*) 'after wpdne, AN=',AN
      dnucl=0.
      dcond=0.

      call coagkernel

      write(16,*)
      write(16,*)
      write(16,210)
  210 format(' Time(s)',5x,'H2SO4',3X,'Nucl_Rt',3x,'Cond_Rt',3x,'Bins:')
      write(16,211) 0.,H2SO4nd,dnucl,dcond,an(1:nsize)
      nstep = nint(time/tstep)
      time = 0.
      do istep=1,nstep
         time = time+tstep
         do ist=1,ngrow
!!            CALL GROW(dcond,tstep/ngrow)
            CALL NUCLEAR(dnucl,TSTEP/ngrow)
         end do
         write(6,*) 'time=',time,' H2SO4=',h2so4nd,dnucl,dcond,' Bins: ',an(1:6)
         H2O=H2Ond/airvd
         H2SO4=H2SO4nd/airvd
         call coag(tstep,ncoag)
         iwrt = 0
         itime = nint(time)
         if(itime .le. 3600) then
            iwrt=1
         else 
            if(mod(itime,3600).eq.0) iwrt=1
         end if
         H2SO4 = H2SO4 + H2SO4rate*tstep
         H2SO4nd = H2SO4*airvd
         if(iwrt.eq.1) write(16,211) time,H2SO4nd,dnucl,dcond,an
  211    format(f12.3,1pE10.3,1pE11.3,1pE11.3,150g9.2)
      end do

      stop
      end
