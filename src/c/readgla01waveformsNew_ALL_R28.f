      program readgla01waveformsNew_ALL_R28
      implicit none
      character*80 gla01
      integer*4 i4gla01_1(1165),i4gla01_2(1165),i4gla01_3(1165)
      integer*2 i2gla01_1(2330),i2gla01_2(2330),i2gla01_3(2330)
      integer*1 i1gla01_1(4660),i1gla01_2(4660),i1gla01_3(4660)
      integer*4 i4gla06(1720)
      integer*2 i2gla06(3440)
      integer i1gla06(6880)
      integer*4 i4gla05(4350)
      integer*2 i2gla05(8700)
      integer i1gla05(17400)
c      integer*4 i_elvflg(40)
c      integer*1 i_attflg2(20),i_frameqf
c      integer*2 i_orbflg(2)
      integer*1 i_nPeaks1(40),i_surftype(1),i_cld1_mswf(1)
      integer*1 i_satNdx(40),i_satCorrFlg(40),i_msCorrFlg(40)
      integer*2 i_gdHt(2),i_SigmaElv(40)
      integer*2 i_kurt2(40),i_skew2(40)
      integer*2 i_srf_ruf(40),i_dTrop(40),i_srf_slope(40)
      integer*2 i_gla01_rectype1,i_gla01_rectype2,i_gla01_rectype3
      integer*2 i_gainSet1064_1(20),i_gainSet1064_2(20)
      integer*2 i_erElv(2),i_spElv(4),i_ldElv(4),i_ocElv(2),i_wTrop(2)
      integer*2 i_satElevCorr(40),i_satNrgCorr(40),i_satPwdCorr(40)
      integer*2 i_gval_rcv(40),i_RecNrgAll(40)
      integer*2 i_msRngCorr(40)
      integer*2 i_Surface_temp,i_Surface_pres,i_Surface_relh
      integer*2 i_wfFitSDev_2(40),i_RecNrgAll_05(40)
      integer*4 i_utctime1(2),i_utctime2(2),i_utctime3(2)
      integer*4 i_lat_pred,i_lon_pred
      integer*4 i_RecNrgLast_EU(40),i_RecNrgAll_EU(40),i_elev(40)
      integer*4 i_utctime(2),i_lat(40),i_lon(40)
      integer*4 i4undf,i_reflctUncorr(40),i_reflCor_atm(1)
      integer*4 i_SigBegOff(40),i_SigEndOff(40),i_cntRngOff(40)
      integer*4 i_isRngOff(40),i_siRngOff(40),i_ldRngOff(40),i_ocRngOff(40)
      integer*4 i_TrshRngOff(40),i_DEM_elv(40),i_refRng(40)
      integer*4 i_SolAng
      integer*4 i_PODFixedPos(6,40)
      real*4 sum,elev,predlat,predlon,gdHt
      real*4 recnrglast_eu,recnrgall_eu,drytrop 
      real*4 srfruf,kurt,skew,srfslope
      real*4 sigbeg,sigend,cntrng,reflctUncorr
      real*4 lat,lon,reflctCor,DEMelev
      real*4 latmin,latmax,sigmaelev
      real*4 isrngoff,sirngoff,ldrngoff,ocrngoff
      real*4 halfmax
      real*4 satElevCor,satNrgCor,satPwdCor
      real*4 recnrgall_06
      real*4 surf_temp,surf_pres,surf_rh
      real*4 SETide,SpTide,LdTide,OcTide,WetTrop_RngCor
      real*4 wfFit_SDev,recnrgall_05
      real*8 refrng,satalt_refrng_rad
      real*8 z_M,z_mm,x_M,x_mm,y_M,y_mm
      real*8 sat_alt,pod_x,pod_y,pod_z,pod_sum
      real*8 Re,pi,EF,jv,jv2,rad
      real*8 utc_time
      integer maxBinPos,fwhm,maxval
      integer irec,irec2,nsat,ishot,i,j,n,k,ngran,status
      integer*2 gain
      integer*2 gain_06
      integer*4 i_rec_ndx,i_rec_ndx1,i_rec_ndx2,i_rec_ndx3
      integer*4 nPeak1,srftype,cld_mswf
      integer*4 sat_ndx,satCorFlg,msCorFlg,msRngCor
      byte i_echo1(200,20),i_echo2(200,20)
      integer*2 echo1(200,40)
      byte i_txwf(48,40)
      integer*2 TXwf(48,40)
      character*72  filename1,filename2,f_read_cha,infile
      character*72  infile2,infile3,outfile,outfile2

      logical f_read_log,output_echoes

      equivalence (i2gla01_1(1),i4gla01_1(1))
      equivalence (i2gla01_2(1),i4gla01_2(1))
      equivalence (i2gla01_3(1),i4gla01_3(1))
      equivalence (i1gla01_1(1),i4gla01_1(1))
      equivalence (i1gla01_2(1),i4gla01_2(1))
      equivalence (i1gla01_3(1),i4gla01_3(1))
      equivalence (i_rec_ndx1,i4gla01_1(1))
      equivalence (i_rec_ndx2,i4gla01_2(1))
      equivalence (i_rec_ndx3,i4gla01_3(1))
      equivalence (i_utctime1(1),i4gla01_1(2))
      equivalence (i_utctime2(1),i4gla01_2(2))
      equivalence (i_utctime3(1),i4gla01_3(2))
      equivalence (i_gla01_rectype1,i2gla01_1(7))
      equivalence (i_gla01_rectype2,i2gla01_2(7))
      equivalence (i_gla01_rectype3,i2gla01_3(7))
      equivalence (i_lat_pred,i4gla01_1(44))
      equivalence (i_lon_pred,i4gla01_1(45))
      equivalence (i_RecNrgAll_EU(1),i4gla01_1(567))
      equivalence (i_RecNrgLast_EU(1),i4gla01_1(607))
      equivalence (i_txwf(1,1),i1gla01_1(2715))
      equivalence (i_gainSet1064_1(1),i2gla01_2(79))
      equivalence (i_gainSet1064_2(1),i2gla01_3(79))
      equivalence (i_echo1(1,1),i1gla01_2(417))
      equivalence (i_echo2(1,1),i1gla01_3(417))

      equivalence (i_wfFitSDev_2(1),i4gla05(3685))
      equivalence (i_RecNrgAll_05(1),i2gla05(8606))

      equivalence (i_rec_ndx,i4gla06(1))
      equivalence (i1gla06(1),i2gla06(1),i4gla06(1))
      equivalence (i_utctime(1),i4gla06(2))
      equivalence (i_lat(1),i4gla06(45))
      equivalence (i_lon(1),i4gla06(85))
      equivalence (i_elev(1),i4gla06(125))
      equivalence (i_SolAng,i4gla06(666))
      equivalence (i_gdHt(1),i2gla06(1339))
      equivalence (i_erElv(1),i2gla06(1341))
      equivalence (i_spElv(1),i2gla06(1343))
      equivalence (i_ldElv(1),i2gla06(1347))
      equivalence (i_ocElv(1),i2gla06(1351))
      equivalence (i_wTrop(1),i2gla06(1353))
      equivalence (i_dTrop(1),i2gla06(1355))
      equivalence (i_surftype,i2gla06(1395))
      equivalence (i_DEM_elv(1),i4gla06(699))
      equivalence (i_refRng(1),i4gla06(739))
      equivalence (i_TrshRngOff(1),i4gla06(779))
      equivalence (i_SigBegOff(1),i4gla06(819))
      equivalence (i_SigEndOff(1),i4gla06(859))
      equivalence (i_cntRngOff(1),i4gla06(899))
      equivalence (i_reflctUncorr(1),i4gla06(939))
      equivalence (i_reflCor_atm,i4gla06(979))
      equivalence (i_SigmaElv(1),i2gla06(1999))
      equivalence (i_kurt2(1),i2gla06(2059))
      equivalence (i_skew2(1),i2gla06(2099))
      equivalence (i_srf_ruf(1),i2gla06(2139))
      equivalence (i_isRngOff(1),i4gla06(1110))
      equivalence (i_siRngOff(1),i4gla06(1150))
      equivalence (i_ldRngOff(1),i4gla06(1190))
      equivalence (i_ocRngOff(1),i4gla06(1230))
      equivalence (i_nPeaks1(1),i4gla06(1270))
      equivalence (i_srf_slope(1),i2gla06(2179))
      equivalence (i_PODFixedPos(1,1),i4gla06(405))
      equivalence (i_cld1_mswf,i2gla06(2564))
      equivalence (i_satNdx(1),i4gla06(1368))
      equivalence (i_satElevCorr(1),i2gla06(2755))
      equivalence (i_satCorrFlg(1),i4gla06(1398))
      equivalence (i_satNrgCorr(1),i2gla06(2815))
      equivalence (i_satPwdCorr(1),i2gla06(2855))
      equivalence (i_gval_rcv(1),i2gla06(2895))
      equivalence (i_RecNrgAll(1),i2gla06(2935))
      equivalence (i_msRngCorr(1),i2gla06(3095))
      equivalence (i_msCorrFlg(1),i4gla06(1568))
      equivalence (i_Surface_temp,i2gla06(3155))
      equivalence (i_Surface_pres,i2gla06(3156))
      equivalence (i_Surface_relh,i2gla06(3157))
c      equivalence (i_elvflg(1),i4gla06(1293))
c      equivalence (i_attflg2(1),i1gla06(5309))
c      equivalence (i_orbflg(1),i2gla06(2666))
c      equivalence (i_frameqf,i1gla06(5330))

c      data irec/7/,irec2/3/,i4undf/2147483647/
      data irec/6/,irec2/3/,i4undf/2147483647/
      gla01='GLA06_in'
      infile = f_read_cha('GLA01_infile','Input filename',' ')
      open(10,file=infile,access='direct',form='unformatted',
     . recl=4660,status='old')
      outfile = f_read_cha('GLA01_outfile','Output filename',' ')
c      output_echoes=f_read_log('GLA01_output_echoes','Output echoes','n')
 
      write(*,*)'Opening :',outfile
      write(*,'(a4,a1,a50)')infile(1:4),'6',infile(6:55)
      infile2='GLA06_428_2115'//infile(15:55)
      write (*,*)'Infile 2: ',infile2
c      infile3='GLA05_428_2113'//infile(15:55)
c      write (*,*)'Infile 3: ',infile3
      open(11,file=infile2,access='direct',form='unformatted',
     . recl=6880,status='old')
      open(19,file=outfile,iostat=status,
     &status='unknown',form='formatted')

      i_gla01_rectype2=2
      latmax=-90.
      latmin=90.

100   continue
      ngran=0
      do while (i_gla01_rectype2.eq.2)
        irec=irec+1
        read(10,rec=irec,err=1999) i4gla01_1
        irec=irec+1
        read(10,rec=irec,err=1999) i4gla01_2
        irec=irec+1
        read(10,rec=irec,err=1999) i4gla01_3
      end do       

        irec=irec+1
        read(10,rec=irec,err=1999) i4gla01_1
*        write(*,*)'Read i4gla01_1 ',i_gla01_rectype1
        irec=irec+1
        read(10,rec=irec,err=1999) i4gla01_2
*        write(*,*)'Read i4gla01_2 ',i_gla01_rectype2
        irec=irec+1
        read(10,rec=irec,err=1999) i4gla01_3
*        write(*,*)'Read i4gla01_3 ',i_gla01_rectype3


      do while (i_utctime(1).lt.i_utctime1(1))
        irec2=irec2+1
        read(11,rec=irec2,err=1999) i4gla06
      end do

      do while (i_utctime(1).lt.i_utctime1(1))
        irec2=irec2+1
        read(11,rec=irec2,err=1999) i4gla05
      end do

*        write(*,*)'UTC : GLA01', i_utctime1(1),' GLA06 ',i_utctime(1)

      ngran=ngran+1
      ishot=0
      sum=0
        do j =1, 40
	 
	  wfFit_SDev=i_wfFitSDev_2(j)/100000.
          recnrgall_05=i_RecNrgAll_05(j)/100.
	  recnrglast_eu=i_RecNrgLast_EU(j)/1000.
	  recnrgall_eu=i_RecNrgAll_EU(j)/1000.
          elev=i_elev(j)/1000.
          srfruf=i_srf_ruf(j)/100.
          kurt=i_kurt2(j)/100.
          skew=i_skew2(j)/100.
          sigbeg=i_SigBegOff(j)/1000.
          sigend=i_SigEndOff(j)/1000.
          cntrng=i_cntRngOff(j)/1000.
          reflctUncorr=i_reflctUncorr(j)/1000000.
          reflctCor=i_reflCor_atm(1)/1000000.
	  lat=i_lat(j)/1000000.
	  lon=i_lon(j)/1000000.
          gdHt=(i_gdHt(1)+i_gdHt(2))/2./100.
	  drytrop=i_dTrop(j)/1000.
	  isrngoff=i_isRngOff(j)/1000.
	  sirngoff=i_siRngOff(j)/1000.
	  ldrngoff=i_ldRngOff(j)/1000.
	  ocrngoff=i_ocRngOff(j)/1000.
	  refrng=i_refRng(j)/1000.
	  srftype=i_surftype(1)
	  DEMelev=i_DEM_elv(j)/100.
	  sigmaelev=i_SigmaElv(j)/1000.
	  cld_mswf=i_cld1_mswf(1)
	  nPeak1=i_nPeaks1(j)  
	  srfslope=i_srf_slope(j)/1000.
          sat_ndx=i_satNdx(j)
	  satElevCor=i_satElevCorr(j)/1000.
          satCorFlg=i_satCorrFlg(j)
	  satNrgCor=i_satNrgCorr(j)/1000.
	  satPwdCor=i_satPwdCorr(j)/1000.
          gain_06=i_gval_rcv(j)
c          recnrgall_06=i_RecNrgAll(j)/100000000000000000.
          recnrgall_06=i_RecNrgAll(j)/100.
          msRngCor=i_msRngCorr(j)
	  msCorFlg=i_msCorrFlg(j)
	  surf_temp=i_Surface_temp/100.
	  surf_pres=i_Surface_pres/10.
	  surf_rh=i_Surface_relh/100.
	  x_M=i_PODFixedPos(1,j)
	  x_mm=i_PODFixedPos(2,j)/1000.
	  y_M=i_PODFixedPos(3,j)
	  y_mm=i_PODFixedPos(4,j)/1000.
	  z_M=i_PODFixedPos(5,j)
	  z_mm=i_PODFixedPos(6,j)/1000.
	  pod_x=(x_M+x_mm)
	  pod_y=(y_M+y_mm)
	  pod_z=(z_M+z_mm)
c         POD is given as a vector (x,y,z). So find magnitude to calculate Satellite Altitude
c	  sqrt:Ans**0.5
	  pod_sum=(pod_x*pod_x)+(pod_y*pod_y)+(pod_z*pod_z)
	  sat_alt=pod_sum**0.5
c	  CALCULATING THE Earth's Radius based on the T/P REFERNECE ELLIPSOID 
c         Earth's radius, rad, at a given geographic latitude (theta) in metres : r = Re(1-EF(sin2(theta)-EFsin2(2theta)))
c         Re is the equatorial radius and EF is Earth's Flattening
	  pi = 3.141592654
	  Re = 6378136.3000
          EF = 0.003352813178
          jv = sin((lat*pi)/180.)*sin((lat*pi)/180.)
          jv2 = sin((2*lat*pi)/180.)*sin((2*lat*pi)/180.)
	  rad = Re*(1-EF*(jv-(EF*(jv2))))
	  satalt_refrng_rad = sat_alt-refrng-rad

          utc_time=dble(i_utctime1(1))+dble(i_utctime1(2))/1000000.+dble(j-1)*0.025

	  SEtide=((i_erElv(1)+i_erElv(2))/2.)/1000.
          SpTide=((i_spElv(1)+i_spElv(2)+i_spElv(3)+i_spElv(4))/4.)/1000.
          LdTide=((i_ldElv(1)+i_ldElv(2)+i_ldElv(3)+i_ldElv(4))/4.)/1000.
          OcTide=((i_ocElv(1)+i_ocElv(2))/2.)/1000.
          WetTrop_RngCor=((i_wTrop(1)+i_wTrop(2))/2.)/1000.

          predlat = float(i_lat_pred)/1000000.
          predlon = float(i_lon_pred)/1000000.
	  if(j.le.20)then
               gain=int(i_gainSet1064_1(j))
             else if(j.gt.20)then 
               gain=int(i_gainSet1064_2(j-20))
           end if	  
	  do i=1,200
             if(j.le.20)then
               echo1(i,j)=int(i_echo1(i,j))
             else
               echo1(i,j)=int(i_echo2(i,(j-20)))
             end if
             if(echo1(i,j).lt.0)echo1(i,j)=echo1(i,j)+256
          end do
	  
	  maxval = -999
	  maxBinPos = -999
	  do i=1,200
             if(echo1(i,j).gt.maxval)then
               maxval=echo1(i,j)
	       maxBinPos = i
             end if
          end do 
	  
	  halfmax = ((dble(maxval))/2.0)
	  fwhm = 0
	  do i=1,200
             if(echo1(i,j).ge.halfmax)then
               fwhm=fwhm+1
             end if
          end do 

	  do k=1,48
             TXwf(k,j)=int(i_txwf(k,j))
 	     if(TXwf(k,j).lt.0)TXwf(k,j)=TXwf(k,j)+256
         end do 
          n=n+1
	
	  if(elev.gt.-100..and.elev.lt.100..and.lat.lt.90.)then
            if(lat.gt.latmax)latmax=lat
            if(lat.lt.latmin)latmin=lat
             write(19,'(2f11.6,i8,i4,2f18.6,18f11.6,4f16.4,2i8,6i4,f6.1,200i4,f18.6,48i4,i4,f9.4,i4,2f9.4,i8,i4,9f13.6)')
     &        lat,lon,irec,j,utc_time,elev,gdHt,
     +        kurt,skew,srfruf,srfslope,sigbeg,sigend,cntrng,
c     +        reflctUncorr,drytrop,recnrglast_eu,recnrgall_eu,recnrgall_06,recnrgall_06
     +        reflctUncorr,recnrglast_eu,recnrgall_eu,recnrgall_06,drytrop,sigmaelev,
     +        isrngoff,sirngoff,ldrngoff,ocrngoff,
     +        refrng,sat_alt,rad,satalt_refrng_rad,
     +        gain,gain_06,srftype,cld_mswf,nPeak1,
     +        maxBinPos,fwhm,maxval,halfmax,
     +        (echo1(k,j),k=1,200),utc_time,(TXwf(k,j),k=1,48),
     +        sat_ndx,satElevCor,satCorFlg,satNrgCor,satPwdCor,
c     +        gain_06,recnrgall_06,msRngCor,msCorFlg,
     +        msRngCor,msCorFlg,
     +        surf_temp,surf_pres,surf_rh,
c     +        SETide,SpTide,LdTide,OcTide,WetTrop_RngCor
     +        SETide,SpTide,LdTide,OcTide,WetTrop_RngCor,wfFit_SDev 
c             write(19,'(2f11.6,i8,i4,2f18.6)')
c     &        lat,lon,irec,j,utc_time,elev
         end if
	end do
      go to 100
c   elevuse flag 4247.5
c12400 format(' index',i10,'elvflg ',10(2x,z8))
 1999 print*,'execution terminating at i_rec_ndx',i_rec_ndx,'nsat',nsat
      write(*,*)'PROCESSING FINISHED ',n,' records output ', 'latmin ',latmin,' latmax ',latmax
      close(10)
      close(11)
      close(19)
      stop
      end
 
