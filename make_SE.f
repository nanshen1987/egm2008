      program make_SE
c
c-----------------------------------------------------------------------
c$
      implicit real*8(a-h,o-z)
      real*4 exclud
	character*128 input
c
      parameter(dlat     = 2.5d0/60.d0,
     &          dlon     = 2.5d0/60.d0,
     &          icell    = 1,
     &          nrows    = nint(180.d0/dlat)+(1-icell),
     &          ncols    = nint(360.d0/dlon),        
     &          exclud   = 9999.0,
     &          excl     = exclud*1.d0,
     &          ibox     = 1)  !  input grid is global? yes/no => 0/1
c-----------------------------------------------------------------------
c     parameter(dnorth   = 30.d0,    !     BOX
c    &          dsouth   = 10.d0,
c    &          dwest    =-10.d0,
c    &          deast    = 10.d0)
c-----------------------------------------------------------------------
c     parameter(dnorth   = 50.d0,    !     CONUS
c    &          dsouth   = 25.d0,
c    &          dwest    =235.d0,
c    &          deast    =295.d0)
c-----------------------------------------------------------------------
      parameter(dnorth   = 90.d0,    !     WORLD
     &          dsouth   =-90.d0,
     &          dwest    =000.d0,
     &          deast    =360.d0)
c-----------------------------------------------------------------------
c     parameter(dnorth   = 55.d0,    !     HIMALAYAS
c    &          dsouth   = 20.d0,
c    &          dwest    = 70.d0,
c    &          deast    =110.d0)
c-----------------------------------------------------------------------
c     parameter(dnorth   = 72.d0,    !     SANDWELL AND SMITH
c    &          dsouth   =-72.d0,
c    &          dwest    =000.d0,
c    &          deast    =360.d0)
c-----------------------------------------------------------------------
      parameter(iglob=1-(int(((deast-dwest)/360.d0)+1.d-10)),
     &          in=nint((90.d0-dnorth)/dlat)+1,
     &          is=nint((90.d0-dsouth)/dlat)+1-icell,
     &          jw=nint(dwest/dlon)+1,
     &          je=nint(deast/dlon)+(1-icell)*iglob,
     &          irow=is-in+1,
     &          jcol=je-jw+1)   
c-----------------------------------------------------------------------
      real*8 scrd(ncols,2),statd(22)  
      real*4 grd1(nrows,ncols),temp(ncols)
      integer*4 lead,tail
c-----------------------------------------------------------------------
c
c     INITIATIZE GRIDS
c
c-----------------------------------------------------------------------
      do i = 1, nrows
        do j = 1, ncols
          grd1(i,j) = exclud
        enddo  !  j
      enddo  !  i      
c-----------------------------------------------------------------------
c
c     READ INPUT GRIDS
c
c-----------------------------------------------------------------------
      if (ibox.eq.0) then
        imin = 1
        imax = nrows
        jmin = 1
        jmax = ncols
        jread = ncols
      else
        imin = in
        imax = is
        jmin = jw
        jmax = je
        jread = jcol
      endif  !  ibox
c-----------------------------------------------------------------------
      print *, ' '
	print *, 'Enter path/input file name:'
      read (5, '(a)') input
	len_str = len_trim(input)
c
      open(unit  = 1,
     &     file   = input,
     &     form   = 'unformatted',
     &     access = 'direct',
     &     recl   = (ncols+2)*4,
     &     status = 'old')
c
      do i = imin, imax
        read(1,rec=i) lead,(temp(jj),jj=1,jread),tail
        jj = 0
        do j = jmin, jmax
          j1 = j
          if (j.lt.1) j1 = j + ncols
          jj = jj+1
          grd1(i,j1) = temp(jj)
        enddo  !  j
      enddo  !  i
      close(1)
c-----------------------------------------------------------------------
c
c     OUTPUT
c
c-----------------------------------------------------------------------
      open(unit  = 1,
     &     file   = input(1:len_str)//'_SE',
     &     form   = 'unformatted',
     &     access = 'direct',
     &     recl   = (ncols+2)*4)
c-----------------------------------------------------------------------
      call endian4(lead)
      call endian4(tail)
      ii = 0
      do i = in, is
        ii = ii+1
        jj = 0
        do j = jw, je
          j1 = j
          if (j.lt.1) j1 = j + ncols
          jj = jj+1
          temp(jj) = grd1(i,j1)
          call endian4(temp(jj))
          grd1(i,j1) = temp(jj)
        enddo  !  j
        write(1,rec=ii) lead,(temp(j),j=1,jj),tail
      enddo  !  i
      close(1) 
c-----------------------------------------------------------------------
      write(6,5000)
      write(6,'(a)')  '  '//input(1:len_str)//'_SE'
      write(6,*) ' -------------------------------------------'
      write(6,*)
      ii = 0
      do i = in, is 
        ii = ii+1
        jj = 0
        do j = jw, je 
          jj = jj + 1
          j1 = j
          if (j.lt.1) j1 = j + ncols
          scrd(jj,1) = grd1(i,j1)
          scrd(jj,2) = 1.d0
        enddo  !  j
        call stats(dnorth,dwest,ii,dlat,dlon,irow,jcol,scrd,excl,0,
     &             statd,icell)    
      enddo  !  i
      write(6,*) ' -------------------------------------------'
      call flush(6)
c-----------------------------------------------------------------------          
 5000 format(//,'c',71('-'),//)
      stop
      end
c-----------------------------------------------------------------------
      subroutine endian4(byte)
c
      integer*1 byte(4),tmp(4)
c
      do j = 1, 4
        tmp(j) = byte(j)
      enddo  !  j
c
      do j = 1, 4
        j1 = 5-j
        byte(j1) = tmp(j)
      enddo  !  j                      
c
      return
      end
c-----------------------------------------------------------------------
c
      SUBROUTINE STATS(TOPLAT,WSTLON,I,GRDN,GRDE,is,NCOLS,DATA,
     $                 EXCLUD,ISIG,STAT,icent)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*20 DLABEL(14)
      CHARACTER*20 SLABEL( 8)
      DIMENSION DATA(NCOLS,2)
      DOUBLE PRECISION DEXCLUD,DG,SD,STAT(22)
      SAVE
      data ix/0/
      DATA PI/3.14159265358979323846D+00/
      DATA DLABEL/'    Number of Values','  Percentage of Area',
     $            '       Minimum Value',' Latitude of Minimum',
     $            'Longitude of Minimum','       Maximum Value',
     $            ' Latitude of Maximum','Longitude of Maximum',
     $            '     Arithmetic Mean','  Area-Weighted Mean',
     $            '      Arithmetic RMS','   Area-Weighted RMS',
     $            '   Arithmetic S.Dev.','Area-Weighted S.Dev.'/
      DATA SLABEL/'       Minimum Sigma',' Latitude of Minimum',
     $            'Longitude of Minimum','       Maximum Sigma',
     $            ' Latitude of Maximum','Longitude of Maximum',
     $            'Arithmetic RMS Sigma','Area-wghtd RMS Sigma'/
C-----------------------------------------------------------------------
      IF(I.EQ.1) THEN
      DEXCLUD=EXCLUD
      DTR=PI/180.D0
      FOURPI=4.D0*PI
      DPR=GRDN*DTR
      DLR=GRDE*DTR
      CAREA=2.D0*DLR*SIN(DPR/2.D0)
      DO 10 K=1,22
      STAT(K)=0.D0
   10 CONTINUE
      STAT( 3)= DEXCLUD
      STAT( 6)=-DEXCLUD
      STAT(15)= DEXCLUD
      STAT(18)= 0.D0
      ENDIF
C-----------------------------------------------------------------------
c
      dlat = toplat -(i-1.d0)*grdn - (grdn/2.d0)*icent
      COLATC=(90.D0-DLAT)*DTR
      AREA=CAREA*SIN(COLATC)
c
c   The following statement avoids over-estimating the total area
c   covered by the grid (of POINT values) when the pole(s) are part
c   of the grid.
c
      if(abs(dlat).eq.90.d0) area=2.d0*dlr*sin(dpr/4.d0)**2
c
C-----------------------------------------------------------------------
      DO 110 J=1,NCOLS
      DLON=WSTLON+(J-1.D0)*GRDE +(GRDE/2.D0)*icent
      DG=DATA(J,1)
      SD=DATA(J,2)
      IF(DG.LT.DEXCLUD) THEN
C-----------------------------------------------------------------------
      STAT( 1)=STAT( 1)+1.D0
      STAT( 2)=STAT( 2)+AREA
      IF(DG.LE.STAT( 3)) THEN
      STAT( 3)=DG
      STAT( 4)=DLAT
      STAT( 5)=DLON
      ENDIF
      IF(DG.GE.STAT( 6)) THEN
      STAT( 6)=DG
      STAT( 7)=DLAT
      STAT( 8)=DLON
      ENDIF
      STAT( 9)=STAT( 9)+DG
      STAT(10)=STAT(10)+DG*AREA
      STAT(11)=STAT(11)+DG**2
      STAT(12)=STAT(12)+DG**2*AREA
      IF(SD.LE.STAT(15)) THEN
      STAT(15)=SD
      STAT(16)=DLAT
      STAT(17)=DLON
      ENDIF
      IF(SD.GE.STAT(18)) THEN
      STAT(18)=SD
      STAT(19)=DLAT
      STAT(20)=DLON
      ENDIF
      STAT(21)=STAT(21)+SD**2
      STAT(22)=STAT(22)+SD**2*AREA
C-----------------------------------------------------------------------
      ENDIF
  110 CONTINUE
C-----------------------------------------------------------------------
      IF(I.NE.is) RETURN
      IF(STAT(1).GT.0.D0) THEN
      STAT( 9)=STAT( 9)/STAT( 1)
      STAT(10)=STAT(10)/STAT( 2)
      STAT(11)=SQRT(STAT(11)/STAT( 1))
      STAT(12)=SQRT(STAT(12)/STAT( 2))
c-----------------------------------------------------------------------
      STAT(13)=0.d0
      STAT(14)=0.d0
      temp1 = STAT(11)**2-STAT( 9)**2
      temp2 = STAT(12)**2-STAT(10)**2
      if (temp1.gt.0.d0) STAT(13)=SQRT(temp1)
      if (temp2.gt.0.d0) STAT(14)=SQRT(temp2)
c-----------------------------------------------------------------------
c     STAT(13)=SQRT(STAT(11)**2-STAT( 9)**2)
c     STAT(14)=SQRT(STAT(12)**2-STAT(10)**2)
c-----------------------------------------------------------------------
      STAT(21)=SQRT(STAT(21)/STAT( 1))
      STAT(22)=SQRT(STAT(22)/STAT( 2))
      STAT( 2)=STAT( 2)/FOURPI*100.D0
      ELSE
      DO 120 J=3,22
      STAT(J)=DEXCLUD
  120 CONTINUE
      ENDIF
C=======================================================================
      NUM=INT(STAT(1))
      WRITE(6,6001) DLABEL(1),NUM
 6001 FORMAT(5X,A20,3X,I11)
      DO 210 K=2,14
      WRITE(6,6002) DLABEL(K),STAT(K)
 6002 FORMAT(5X,A20,3X,F16.4)
  210 CONTINUE
      WRITE(6,6003)
 6003 FORMAT(' ')
      IF(ISIG.EQ.1) THEN
      DO 220 K=1,8
      WRITE(6,6002) SLABEL(K),STAT(K+14)
  220 CONTINUE
      ENDIF
C=======================================================================
      RETURN
      END