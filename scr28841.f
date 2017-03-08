      program readmin1
      implicit real*8(a-h,o-z)
      parameter(exclud=9999.d0,
     $   dlatg=2.5d0/60.d0,dlong=2.5d0/60.d0,nrowsg=04320,ncolsg=08640)
         real*4 data(ncolsg)
      dimension sdat(ncolsg,2),stati(22)
c-----------------------------------------------------------------------
c
c   Read-in Gravity Anomaly 2.5x2.5 minute point value file.
c
      write(6,6001)
 6001 format(///15x,'Statistics of Gravity Anomaly Values  (mGal)',/)
      do i = 1, nrowsg
        read(1) (  data(j),j=1,ncolsg)
        rlat = 90.d0 - (i-1.d0)*dlatg - dlatg*0.5d0
        do j = 1, ncolsg
          clon = (j-1.d0)*dlong + dlong*0.5d0
c
c   Print a few values.
c
          if(i.eq.2001.and.j.lt.20) then  !  top row
            write(6,6101) rlat,clon,data(j)
 6101       format(5x,2f15.10,f15.5)
          endif  !  top row
          sdat(j,1) = data(j)
          sdat(j,2) = 1.d0
        enddo  !  j
        call stats(90.d0,0.d0,i,dlatg,dlong,nrowsg,ncolsg,sdat,exclud,0,
     $             stati)
      enddo  !  i
c-----------------------------------------------------------------------
c
c   Read-in xi Deflection of Vertical 2.5x2.5 minute point value file.
c
      write(6,6002)
 6002 format(///15x,'Statistics of xi DOV Values  (arc-second)',/)
      do i = 1, nrowsg
        read(2) (  data(j),j=1,ncolsg)
        rlat = 90.d0 - (i-1.d0)*dlatg - dlatg*0.5d0
        do j = 1, ncolsg
          clon = (j-1.d0)*dlong + dlong*0.5d0
c
c   Print a few values.
c
c         if(i.eq.2001.and.j.lt.20) then  !  top row
          if( (rlat .gt. 39.0) .and. (rlat .lt. 42.0) ) then
            if(clon.gt.28.0 .and. clon.lt.35.0) then
              write(6,6101) rlat,clon,data(j)
            endif
          endif  !  top row
          sdat(j,1) = data(j)
          sdat(j,2) = 1.d0
        enddo  !  j
        call stats(90.d0,0.d0,i,dlatg,dlong,nrowsg,ncolsg,sdat,exclud,0,
     $             stati)
      enddo  !  i
c-----------------------------------------------------------------------
c
c   Read-in eta Deflection of Vertical 2.5x2.5 minute point value file.
c
      write(6,6003)
 6003 format(///15x,'Statistics of eta DOV Values  (arc-second)',/)
      do i = 1, nrowsg
        read(3) (  data(j),j=1,ncolsg)
        rlat = 90.d0 - (i-1.d0)*dlatg - dlatg*0.5d0
        do j = 1, ncolsg
          clon = (j-1.d0)*dlong + dlong*0.5d0
c
c   Print a few values.
c
          if(i.eq.   1.and.j.lt.20) then  !  top row
            write(6,6101) rlat,clon,data(j)
          endif  !  top row
          sdat(j,1) = data(j)
          sdat(j,2) = 1.d0
        enddo  !  j
        call stats(90.d0,0.d0,i,dlatg,dlong,nrowsg,ncolsg,sdat,exclud,0,
     $             stati)
      enddo  !  i
c-----------------------------------------------------------------------
      stop
      end
C
C
C
      SUBROUTINE STATS(TOPLAT,WSTLON,I,GRDN,GRDE,NROWS,NCOLS,DATA,
     $                 EXCLUD,ISIG,STAT)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*20 DLABEL(14)
      CHARACTER*20 SLABEL( 8)
      DIMENSION DATA(NCOLS,2)
      DOUBLE PRECISION DEXCLUD,DG,SD,STAT(22)
      SAVE
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
      DLAT=TOPLAT-(I-1.D0)*GRDN-GRDN/2.D0
      COLATC=(90.D0-DLAT)*DTR
      AREA=CAREA*SIN(COLATC)
C-----------------------------------------------------------------------
      DO 110 J=1,NCOLS
      DLON=WSTLON+(J-1.D0)*GRDE+GRDE/2.D0
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
      IF(I.NE.NROWS) RETURN
      IF(STAT(1).GT.0.D0) THEN
      STAT( 9)=STAT( 9)/STAT( 1)
      STAT(10)=STAT(10)/STAT( 2)
      STAT(11)=SQRT(STAT(11)/STAT( 1))
      STAT(12)=SQRT(STAT(12)/STAT( 2))
      STAT(13)=SQRT(STAT(11)**2-STAT( 9)**2)
      STAT(14)=SQRT(STAT(12)**2-STAT(10)**2)
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
 6001 FORMAT(/5X,A20,3X,I11)
      DO 210 K=2,14
      WRITE(6,6002) DLABEL(K),STAT(K)
 6002 FORMAT(5X,A20,3X,F15.3)
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
