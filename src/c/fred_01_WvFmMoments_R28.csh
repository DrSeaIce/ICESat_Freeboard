################################################################################################################################################################
# Sinead Farrell, CPOM, 17.09.05
# modified 10/10/05
# modified 28/02/07
#
# Extracting: Tx and Rx pulse data and specific parameters from GLAS text file constructed using FORTRAN code:
#  readgla01waveformsNew_ALL_R28.f
#
# Script reads in each individual GLAS track, Performs C calculations, Outputs data files, Moves onto next GLAS track. 
#
# Infile is a GLAS track text file 
# The GLAS input file has the form: lat, lon, irec, j,utc_time,elev,gdHt,kurt,skew,srfruf,srfslope,sigbeg,sigend,cntrng,reflctUncorr,
# 	recnrglast_eu,recnrgall_eu,recnrgall_06,drytrop,sigmaelev,isrngoff,sirngoff,ldrngoff,ocrngoff,refrng,sat_alt,rad,satalt_refrng_rad,
#       gain,gain_06,srftype,cld_mswf,nPeak1,maxBinPos,fwhm,maxval,halfmax,(echo1(k,j),k=1,200) ,utc_time,(TXwf(k,j),k=1,48)
#       sat_ndx,satElevCor,satCorFlg,satNrgCor,satPwdCor,
#       msRngCor,msCorFlg,surf_temp,surf_pres,surf_rh,SETide,SpTide,LdTide,OcTide,WetTrop_RngCor,wfFit_SDev 
#
# Infile has 54 parameter fields. Fields 38 - 237 incl. are the Rx waveform samples. Field 238 is utctime (again).
# Fields 239 - 286 are the Tx waveform samples. Fields 287 - 302 are the new R28 parameter fields.
# The RX waveform has NOT been re-ordered at this stage
#
# OUTPUT file from csh script: params.R28.dat
# INPUT files for C code: TXpulse.dat and RXpulse.dat
# OUTPUT files from C code: TXmoments.dat, RXmoments.dat, TX-nTHnorm_array.dat, RX-nTHnorm_array.dat  
#
################################################################################################################################################################
#!/bin/csh

echo
echo "setting constants "

#Arctic
set latmin = 65.0
set latmax = 90.0
set lonmin = 0.0
set lonmax = 360.0
##Antarctic
#set latmin = -90.0
#set latmax = -60.0
#set lonmin = 0.0
#set lonmax = 360.0
#echo 'Grid Box :',$lonmin,$lonmax,$latmin,$latmax
echo "-------------------------------------------------------------- "

echo
echo "reading the GLAS data file..."
echo

foreach f ( GLA01*.wvfm.gz )
sleep 2
gunzip $f

#foreach f ( GLA01*.wvfm.Z )
#sleep 2
#uncompress $f

echo "setting infile "
setenv infile $f:r
echo $f
echo $infile
echo
echo "lat(deg),lon(deg),irec,j,utctime(s),elev(m),gdHt(m),kurt,skew,srfruf(m),srfslope,sigbeg(m),sigend(m),cntrng(m),reflctUncorr,recnrglast_eu(fJ),recnrgall_eu(fJ),recnrgall_06(fJ),drytrop,sigmaelev,isrngoff,sirngoff,ldrngoff,ocrngoff,refrng,sat_alt,rad,satalt_refrng_rad,gain,gain_06,srftype,cld_mswf,nPeak1,maxBinPos,fwhm,maxval,halfmax,(echo1(k,j),k=1,200),utc_time,(TXwf(k,j),k=1,48),sat_ndx,satElevCor,satCorFlg,satNrgCor,satPwdCor,msRngCor,msCorFlg,surf_temp,surf_pres,surf_rh,SETide(m),SpTide(m),LdTide(m),OcTide(m),WetTrop_RngCor(m),wfFit_SDev "
#       1        2        3  4     5      6        7      8     9       10     11         12       13        14         15             16                17              18             19        20       21       22       23       24      25    26      27         28         29   30        31     32       33      34      35   36     37        38-237             238          239-286        287       288        289     290         291      292       293      294      295       296    297       298        299      300         301               302
echo "head "
head -1 $infile
echo 
echo "tail "
tail -1 $infile
echo
echo "word count "
wc $infile
echo "-------------------------------------------------------------- "


echo
echo "extracting params.R28.dat file for data within grid box "
echo "please wait ..."

# 1) select only lines which have 302 fields
# dataname of each field:
#     "lat(deg),lon(deg),irec,j,utctime(s),elev(m),gdHt(m),kurt,skew,srfruf(m),srfslope,sigbeg(m),sigend(m),cntrng(m),reflctUncorr,recnrglast_eu(fJ),recnrgall_eu(fJ),recnrgall_06(fJ),drytrop,sigmaelev,isrngoff,sirngoff,ldrngoff,ocrngoff,refrng,sat_alt,rad,satalt_refrng_rad,gain,gain_06,srftype,cld_mswf,nPeak1,maxBinPos,fwhm,maxval,halfmax,(echo1(k,j),k=1,200),utc_time,(TXwf(k,j),k=1,48),sat_ndx,satElevCor,satCorFlg,satNrgCor,satPwdCor,msRngCor,msCorFlg,surf_temp,surf_pres,surf_rh,SETide(m),SpTide(m),LdTide(m),OcTide(m),WetTrop_RngCor(m),wfFit_SDev "
#       1        2        3  4     5      6        7      8     9       10     11         12       13        14         15             16                17              18             19        20       21       22       23       24      25    26      27         28         29   30        31     32       33      34      35   36     37        38-237             238          239-286        287       288        289     290         291      292       293      294      295       296    297       298        299      300         301               302
awk '{ if (NF==302 && $1 >= '$latmin' && $1 <= '$latmax' && $2 >= '$lonmin' && $1 <= '$lonmax') {print NR,$1,$2,$5,$6,$18,$30,$15,$19,$7,$288,$289,$290,$295,$297,$298,$299,$300,$301,$302 } }' $infile > params.R28.dat
echo "line no.,lat,lon,utctime,elev,recnrgall_06,gain_06,reflctUncorr,drytrop,gdHt,satElevCor,satCorFlg,satNrgCor,surf_pres,SETide(m),SpTide(m),LdTide(m),OcTide(m),WetTrop_RngCor(m),wfFit_SDev "
head -1 params.R28.dat
echo
tail -1 params.R28.dat
echo
#minmax params.R28.dat
wc params.R28.dat
echo
echo
echo "------------------ "

echo
echo "extracting Rx pulse data within grid box "
echo "please wait ..."
awk '{ for (i=38; i<=237; i++) { if (NF==302 && $1 >= '$latmin' && $1 <= '$latmax' && $2 >=  '$lonmin' && $2 <= '$lonmax') {printf ("%d  ", $i) } } ; { if ($1 >= '$latmin' && $1 <= '$latmax' && $2 >=  '$lonmin' && $2 <= '$lonmax')  {print "" }} }' $infile > RXpulse.dat
echo "Rx_pulse((RXwf(k,j),k=1,200) "
echo "head "
head -1 RXpulse.dat
echo
echo "tail "
tail -1 RXpulse.dat
echo
#minmax RXpulse.dat
wc RXpulse.dat
echo
echo
echo "------------------ "

echo
echo "extracting Tx pulse data within grid box "
echo "please wait ..."
awk '{ for (i=239; i<=286; i++) { if (NF==302 && $1 >= '$latmin' && $1 <= '$latmax' && $2 >=  '$lonmin' && $2 <= '$lonmax') {printf ("%d  ", $i) } } ; { if ($1 >= '$latmin' && $1 <= '$latmax' && $2 >=  '$lonmin' && $2 <= '$lonmax') {print "" }} }' $infile > TXpulse.dat
echo "Tx_pulse((TXwf(k,j),k=1,48) "
echo "head "
head -1 TXpulse.dat
echo
echo "tail "
tail -1 TXpulse.dat
echo
#minmax TXpulse.dat
wc TXpulse.dat
echo
echo
echo "------------------ "
echo 'Finished parameter extract'
echo "-------------------------------------------------------------- "



echo "starting C program now... "
#
# can run this code on S5 or M12
#
# Cprogram name: Fred_moments_v3.c
# compile and run C code:
./makec_fred_moments
#makec_fred_moments
echo 
echo "-------------------------------------------------------------- "


echo
echo "checking output file from C code: "
echo

echo
echo "TXmoments.dat "
echo "head "
echo "line no., TX_maxval, TX_maxvalBinPos, TX_noise_mean, TX_noise_sd, TX_nTH, TX_numelems_abv_nTH, TX_mean, TX_sd, TX_skew, TX_kurt "
head -1 TXmoments.dat
echo

echo
echo "RXmoments.dat "
echo "head "
echo "line no., RX_maxval, RX_maxvalBinPos, RX_noise_mean, RX_noise_sd, RX_nTH, RX_numelems_abv_nTH, RX_mean, RX_sd, RX_skew, RX_kurt "
head -1 RXmoments.dat
echo

echo
echo "TX-nTHnorm_array.dat "
echo "head "
echo "ARRAY: Tx_norm[48] "
head -1 TX-nTHnorm_array.dat
echo

echo
echo "RX-nTHnorm_array.dat "
echo "head "
echo "ARRAY: Rx_norm[200] "
head -1 RX-nTHnorm_array.dat
echo

echo "-------------------------------------------------------------- "


echo
echo "renaming dat files"
mv params.R28.dat $infile.R28.params
mv TXpulse.dat $infile.TX.array
mv RXpulse.dat $infile.RX.array
mv TXmoments.dat $infile.TXmoments
mv RXmoments.dat $infile.RXmoments
mv TX-nTHnorm_array.dat $infile.TX-nTHnorm_array
mv RX-nTHnorm_array.dat $infile.RX-nTHnorm_array
mv infile_lines.flag $infile.infile_lines.flag
mv outfile_lines.flag $infile.outfile_lines.flag
echo "-------------------------------------------------------------- "

echo
echo "zipping dat files"
gzip $infile
gzip $infile.R28.params
gzip $infile.TX.array
gzip $infile.RX.array
gzip $infile.TXmoments
gzip $infile.RXmoments
gzip $infile.TX-nTHnorm_array
gzip $infile.RX-nTHnorm_array
echo "-------------------------------------------------------------- "

echo
echo "moving dat files"
mv $infile.R28.params.gz /cpnet/xserve3a/slf/GLA_PRO_TMP/L3E/R28/Fred/L3E_FM06_data/wvfm_params
mv $infile.TX.array.gz /cpnet/xserve3a/slf/GLA_PRO_TMP/L3E/R28/Fred/L3E_FM06_data/TX_arrays
mv $infile.RX.array.gz /cpnet/xserve3a/slf/GLA_PRO_TMP/L3E/R28/Fred/L3E_FM06_data/RX_arrays
mv $infile.TXmoments.gz /cpnet/xserve3a/slf/GLA_PRO_TMP/L3E/R28/Fred/L3E_FM06_data/TX_moment_stats
mv $infile.RXmoments.gz /cpnet/xserve3a/slf/GLA_PRO_TMP/L3E/R28/Fred/L3E_FM06_data/RX_moment_stats
cp $infile.TX-nTHnorm_array.gz /cpnet/xserve3a/slf/GLA_PRO_TMP/L3E/R28/Fred/L3E_FM06_data/TX-nTHnorm_arrays
cp $infile.RX-nTHnorm_array.gz /cpnet/xserve3a/slf/GLA_PRO_TMP/L3E/R28/Fred/L3E_FM06_data/RX-nTHnorm_arrays
echo "-------------------------------------------------------------- "


end
echo "-------------------------------------------------------------- "

