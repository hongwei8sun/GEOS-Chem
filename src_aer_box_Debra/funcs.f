      FUNCTION LENGTH(STRING,MAX)
      INTEGER MAX
      CHARACTER*1 STRING(MAX)
      DO L=MAX,1,-1
         IF(STRING(L).NE.' ') GO TO 10
      END DO
      LENGTH=0
      RETURN
   10 LENGTH=L
      RETURN
      END

      logical function isnumber(strval,maxc)
      integer maxc
      character*1 strval(maxc),numb(16)
      do i=1,10
         write(numb(i),110) i-1
  110    format(i1)
      end do
      numb(11)=' '
      numb(12)='e'
      numb(13)='E'
      numb(14)='+'
      numb(15)='-'
      numb(16)='.'
!##      write(6,*) numb
      isnumber=.false.
      do 10 i=1,maxc
         do j=1,10
            if(strval(i).eq.numb(j)) isnumber=.true.
         end do
   10 continue
      if(.not.isnumber) return
      do 25 i=1,maxc
         do j=1,16
            if(strval(i).eq.numb(j)) go to 25
         end do
         isnumber=.false.
         return
   25 continue
      return
      end

      function atoi(arg,n)
      character*60 arg,junk
      junk=' '
      l=length(arg,n)
      junk(60-l+1:60)=arg(1:l)
      read(junk,'(i60)') number
      atoi=number
      return
      end

      function atof(arg,n)
      character*60 arg,junk
      junk=' '
      l=length(arg,n)
      junk(60-l+1:60)=arg(1:l)
      read(junk,'(f60.0)') anumber
      atof=anumber
      return
      end

      function ifroma(arg,n)
      character*60 arg,junk
      junk=' '
      l=length(arg,n)
      junk(60-l+1:60)=arg(1:l)
      read(junk,'(i60)') number
      ifroma=number
      return
      end

      function ffroma(arg,n)
      character*60 arg,junk
      junk=' '
      l=length(arg,n)
      junk(60-l+1:60)=arg(1:l)
      read(junk,'(f60.0)') anumber
      ffroma=anumber
      return
      end

      function instring(string,length,char)
      character*1 char,string(length)
!##  if character 'char' is contained in string 'string' of length 'length'
!##  then function returns position of first occurrence of 'char'
!##  otherwise function returns 0.
      do i=1,length
         if(string(i).eq.char) then
            instring=i
            return
         end if
      end do
      instring=0
      return
      end

      logical function insubstring(string,length,sub,L2)
      character*(*) sub,string
!##  if string 'sub' of length 'L2' is contained in string 
!##  'string' of length 'length'
!##  then function returns position of first occurrence of 'sub'
!##  otherwise function returns 0.
      do i=1,length-L2+1
         if(string(i:i+L2-1).eq.sub(1:L2)) then
            insubstring=.true.
            return
         end if
      end do
      insubstring=.false.
      return
      end
