      SUBROUTINE COAG(tstep0,ncoag)
!
!         $Id: coag.f,v 1.14 2011/12/28 14:48:06 dkweis Exp $
!
      include 'aerosol.i'
      DIMENSION DN(NSIZE),DNAL(NSZAL),DNM(NSZAL),DNSO4(NSZAL),IOUT(NSIZE)
      dimension dso4(nsize),dmso4(nszal)
      character*24 aerid,alumid,IDENT,almixid,also4id
      logical nomixed,alsalt
      SAVE IOUT
      COMMON /CAPTN/MODNO,NODAY,NSTOP,IYEAR,IDOFY
      COMMON /SIZE1/R(NSIZE),DR(NSIZE),VOL(NSIZE),SAD1(NSIZE)
      COMMON /SIZESOL/Rg(NSZAL),dRg(nszal),cores(NSZAL), &
     &    VOLAL(NSZAL),SADAL1(NSZAL)
      COMMON /SIZEWET/Rgwet(NSZAL),dRgwet(nszal),SAD1wet(NSZAL)
      COMMON /AIR/TMP,PR,AIRVD
      COMMON /INFOR/WP,DEN,BVP,ST
      COMMON /NUMBER1/ AN(NSIZE)
      COMMON /SOLIDS/ALN(NSZAL),ALM(NSZAL),ALSO4(NSZAL)
      common /AEROID/AERID(nsize),ALUMID(NSZAL),almixid(nszal),also4id(nszal)
      COMMON /KERNEL/AKSS(NSIZE,NSIZE),AKALAL(NSZAL,NSZAL), &
     &  AKSAL(NSIZE,NSZAL),AKSM(NSIZE,NSZAL),  &
     &  AKALM(NSZAL,NSZAL),AKMM(NSZAL,NSZAL)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                          cccccccc
!  This subroutine calculates the change of N due to       cccccccc
!      coagulation using explicit scheme.                  cccccccc
!    This routine may work correctly only for Vrat >=2     cccccccc
!                                                          cccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      data first/1/

      if (first.ne.0) then
         write(6,*)
         write(6,*) 'Coag:  tstep0=',tstep0,'  ncoag=',ncoag
         write(6,*)
!         WRITE(6,*) 'SULFATE COAG KERNELS'
!         write(6,410) (k,k=1,nsize)
  410    format(3x,40i8)
!         do j=1,nsize
!            WRITE(6,420) J,(AKSS(J,K),k=1,nsize)
!  420       format(i3,40(1pe9.2))
!         END DO
      end if
      first=0

!
      DO ISTEP=1,NCOAG
!!        CALL COAGKERNEL
!           CHECK FOR NEGATIVE NUMBER DENSITIES AND FILL FROM
!              NEXT LARGER BIN, WHICH HAS DOUBLE VOLUME
         DO K=1,NSIZE     
            IF(AN(K).LT.0.) THEN
               if(an(k).lt.-1.E-20) then
                  write(16,*) 'Warning:  negative aerosol. start coag'
                  write(16,*) '  lat=',ii,'  lev=',jj,'  ISTEP=',ISTEP
                  write(16,*) '  ISZ=',K,'  AN=',AN(K)
                  STOP 'NEGATIVE AEROSOL IN COAG'
               ELSE
                  an(k)=0.
               end if
            END IF
         END DO
         tstep=tstep0/ncoag
         ncount=0

  250    CONTINUE
         DO I=1,NSIZE
            DN(I)=0.
         END DO
         do i=1,nszal
            DNAL(I)=0.
            DNM(I)=0.
            DNSO4(I)=0.
         end do

!       dso4 is number of H2SO4 molecules per sulfate particle in bin k
      do k=1,nsize
         dso4(k)=vol(k)*den*WP*0.01/98.*6.02E23
      end do

!       dmso4 is number of H2SO4 molecules per mixed particle in bin k
         do k=1,nszal
            if(alm(k).gt.0.) then
               dmso4(k) = also4(k)/alm(k) ! molec H2SO4 per particle
            else
               dmso4(k)=0.
            end if
         end do
!         if(ii.eq.1 .and. jj.eq.1 .and. istep.eq.1) then
!            write(16,*) 'dso4=',dso4
!            write(16,*) 'dmso4=',dmso4
!         end if
!!         write(16,*) 'Coag
!!         write(16,*) 'alm=',(alm(isz),isz=1,nszal)
!!         write(16,*) 'als=',(also4(isz),isz=1,nszal)
!!         write(16,*) 'dnso4=',dnso4

!  sulfate-sulfate coagulation
         DO 200 I=1,NSIZE-1
            DO 300 J=1,I
               S=1.
               IF (I.EQ.J) S=0.5 
               F=VOL(J)/VOL(I)
!               if(an(i).gt.1.e-18.and. an(j).gt.1.e-18) then
                  AL=AN(I)*AN(J)*S*AKSS(I,J)
                  DN(J)=DN(J)-AL
                  if(i.lt.nsize) then
                     DN(I)=DN(I)-AL*F/(Vrat-1.)
                     DN(I+1)=DN(I+1)+AL*F/(Vrat-1.)
                  else
                     DN(I)=DN(I)+AL*F
                  end if
!               end if
  300       CONTINUE
  200    CONTINUE


!  alumina-alumina coagulation
         DO I=1,NSZAL-1
            DO J=1,I
               S=1.
               IF (I.EQ.J) S=0.5 
               F=VOLAL(J)/VOLAL(I)
!               if(aln(i).gt.1.e-18 .and. aln(j).gt.1.e-18) then
                  AL=ALN(I)*ALN(J)*S*AKALAL(I,J)
                  DNAL(J)=DNAL(J)-AL
                  if(i.lt.nszal) then
                     DNAL(I)=DNAL(I)-AL*F/(Vrat-1.)
                     DNAL(I+1)=DNAL(I+1)+AL*F/(Vrat-1.)
                  else
                     DNAL(I)=DNAL(I)+AL*F
                  END IF
!               END IF
            END DO
         END DO
!##         write(16,*) 'after al-al coag, alm=',(alm(isz),isz=1,nszal)
!##         write(16,*) 'after al-al coag, als=',(also4(isz),isz=1,nszal)
!##         write(16,*) 'after al-al coag, dnso4=',dnso4

         if(nomixed) go to 500
!##         write(16,*) 'Coag Point '


  500    continue
!##         do k=1,nszal
!##            if(dnm(k).gt.0.) stop 'coag: dnm>0'
!##            if(dnso4(k).gt.0.) stop 'coag: dnso4>0'
!##         end do
         tmin=tstep
         do k=1,nsize
            if(dn(k).lt.0.) then
               if((-dn(k)*tstep) .gt. 0.9*an(k)) then
                  tnew=-0.9*an(k)/dn(k)
                  if(tnew.lt.tmin) tmin=tnew
!                  write(16,*) 'time step too large for sulfate point ',k
!                  write(16,*) 'an=',an(k),'  dn=',dn(k)*tstep
!                  write(16,*) 'sulfate coag tstep=',tstep,'  tnew=',tnew
!                  call flush(16)
               end if
            end if
         end do
         do k=1,nszal
            if(dnal(k).lt.0.) then
               if((-dnal(k)*tstep) .gt. 0.9*aln(k)) then
                  tnew=-0.9*aln(k)/dnal(k)
                  if(tnew.lt.tmin) tmin=tnew
!                  write(16,*) 'time step too large for alumina point ',k
!                  write(16,*) 'aln=',aln(k),'  dnal=',dnal(k)*tstep
!                  write(16,*) 'alum coag tstep=',tstep,'  tnew=',tnew
!                  call flush(16)
               end if
            end if
            if(dnm(k).lt.0.) then
               if((-dnm(k)*tstep) .gt. 0.9*alm(k)) then
                  tnew=-0.9*alm(k)/dnm(k)
                  if(tnew.lt.tmin) tmin=tnew
!                  write(16,*) 'time step too large for mixedal point ',k
!                  write(16,*) 'alm=',alm(k),'  dnm=',dnm(k)*tstep
!                  write(16,*) 'mixedal coag tstep=',tstep,'  tnew=',tnew
!                  call flush(16)
               end if
            end if
            if(dnso4(k).lt.0.) then
               if((-dnso4(k)*tstep) .gt. 0.9*also4(k)) then
                  tnew=-0.9*also4(k)/dnso4(k)
                  if(tnew.lt.tmin) then
                     tmin=tnew
                     write(16,*) 'time step too large for mixed sulfate point ',k
                     write(16,*) 'also4=',also4(k),'  dnso4=',dnso4(k)*tstep
                     write(16,*) 'mixed sulfate coag tstep=',tstep,'  tnew=',tnew
                     call flush(16)
                  end if
               end if
            end if
         end do
!##         if(tmin.lt.tstep) then
!##            write(16,*) 'coag timestep reduced from ',tstep,' to ',tmin, &
!##      &       'for point ','  iday, istep=',idofy,istep
!##            write(16,*) 'time remaining = ',tstep-tmin
!##            if(tmin/tstep.lt.0.05) tmin=0.05*tstep
!##         endif
!#         write(6,*) 'coag, using tmin=',tmin,' iter# ',ncount
         DO K=1,NSIZE
            AN(K)=AN(K)+DN(K)*TMIN
            if(an(k).lt.0.) then
!   if negative sulfate aerosol concentration, borrow from mixed particles
!##               borrow = -an(k)*dso4(k)
!##               write(16,*) 'negative AN for ',k
!##               write(16,*) 'borrowing ',borrow
!##               do ial=nszal,1,-1
!##                  reduce=min(also4(ial),borrow)
!##                  also4(ial)=also4(ial)-reduce
!##                  borrow=borrow-reduce
!##                  if(reduce.gt.0.) write(16,*) 'borrowing from mixed bin ',ial,reduce,' remaining=',borrow
!##               end do
               an(k)=0.
            end if
         END DO
         do k=1,nszal
            aln(k)=aln(k)+dnal(k)*tmin
            alm(k)=alm(k)+dnm(k)*tmin
            also4(k)=also4(k)+dnso4(k)*tmin
!##            also4(k)=also4(k)+dnso4(k)*tmin*(1.-ratiodnso4)
!##            CNa2SO4=cNa2SO4+dnso4(k)*tmin*ratiodnso4
!##            CNa=CNa-2.*dnso4(k)*tmin*ratiodnso4
            if(aln(k).lt.0.) aln(k)=0.
            if(alm(k).lt.0.) alm(k)=0.
            if(also4(k).le. 0.0) then
               also4(k)=0.
               if(abs(Df-Dfwet).lt.0.001) then
                  aln(k)=aln(k)+alm(k)
                  alm(k)=0.
               end if
            end if
            if(also4(k).ge.0. .or. also4(k).le.0.) then
               continue
            else
               write(16,*) 'end of coag loop, istep=',istep
               write(16,*) 'also4=NaN in coag at ',k
               write(16,*) 'dnso4=',dnso4(k),'  also4=',also4(k)
               stop 'error in coag-explicit.f'
            end if
         end do
!##         write(16,*) 'after update for point '
!##         totalso4=0.
!##         do k=1,nszal
!##            totalso4=totalso4+also4(k)
!##         end do
!##         write(16,*) 'CNa=',CNa,'  CNa2SO4=',CNa2SO4,'  tot also4=',totalso4
!##         if(iNa.gt.0 .and. iNa.le.NTSP) then
!##            allssp(iNa-NFSP)=CNa/airvd
!##            allssp(iNa2SO4-NFSP)=CNa2SO4/airvd
!##         end if
         TSTEP=TSTEP-TMIN
         ncount=ncount+1
         if(ncount.gt.150) then
            write(16,*) 'Stopping at point ','  day ',idofy,istep
            write(16,*) 'AN=',AN(:)
            WRITE(16,*) 'ALN=',ALN(:)
            WRITE(16,*) 'ALM=',ALM(:)
            WRITE(16,*) 'ALSO4=',ALSO4(:)
            stop 'coag too many iterations'
         end if
         if(tstep.gt.1.e-3) write(6,*) 'coag iteration# ',ncount,' time step=',tmin,', remaining time: ',tstep
         IF(TSTEP.GT.1.e-3) GO TO 250
!         write(16,*) 'coag:  # 250 iterations = ',ncount
!           CHECK FOR NEGATIVE NUMBER DENSITIES AND FILL FROM
!              NEXT LARGER BIN, WHICH HAS DOUBLE VOLUME
         DO K=1,NSIZE     
            IF(AN(K).LT.0.) THEN
               if(an(k).lt.-1.E-20) then
                  write(16,*) 'Warning:  time step too large for coag'
                  write(16,*) '  lat=',ii,'  lev=',jj
                  write(16,*) '  ISZ=',K,'  AN=',AN(K)
                  STOP 'TIME STEP TOO LARGE IN COAG'
               ELSE
                  an(k)=0.
               end if
            END IF
         END DO

  100 CONTINUE             ! END OF LOOP OVER GRID POINTS
      END DO               ! end loop over NCOAG substeps
!##      DO ISZ=1,NSIZE
!##         IF(ISZ.LE.9) THEN
!##            WRITE(IDENT,808) ISZ
!##  808       FORMAT('COAG RATE FOR BIN',I1)
!##         ELSE
!##            WRITE(IDENT,809) ISZ
!##  809       FORMAT('COAG RATE FOR BIN',I2)
!##         END IF
!##         CALL WRTSSP(DN(1,1,isz),1,NLT,NHT,IDENT,iout,72)
!##      END DO
      RETURN
      END 

