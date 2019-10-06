################################################################################################################################################################
# Sinead Farrell, CPOM, 11.10.05
# modified 
#
#
# Script reads in each individual GLAS track, Performs C calculations, Outputs data files, Moves onto next GLAS track. 
#
# Infile1 is a data file containing GLAS Tx pulses. Each has had its noise floor removed and has been normalised.
# Infile2 is a data file containing GLAS Rx pulses. Each has had its noise floor removed and has been normalised.
#
# these are read into a C code and the Cross-correlation between Tx and Rx is calculated.
#
# INPUT files for C code: TX-nTHnorm_array.dat, RX-nTHnorm_array.dat  
# OUTPUT files from C code: TX_padded_array.dat, RX_padded_array.dat, xcorrel_TX_RX_array.dat.old, xcorrel_TX_RX_rar_array.dat.old, xcorrel_TX_RX_stats.dat.old
#
################################################################################################################################################################
#!/bin/csh
echo
echo "reading GLAS pulse files..."
echo

#foreach f1 ( GLA01*.wvfm.TX-nTHnorm_array.gz )
foreach f1 ( GLA01*.wvfm.TX-nTHnorm_array )
echo 
echo " ---new echo pair--- "
#echo "unzipping "
echo $f1
#gunzip $f1


echo "setting infile "
#setenv infile $f1:r:r
setenv infile $f1:r
echo $infile
echo

#echo "unzipping "
#echo $infile.RX-nTHnorm_array.gz
#gunzip $infile.RX-nTHnorm_array.gz
echo $infile.RX-nTHnorm_array

echo "copying infiles "
mv $infile.TX-nTHnorm_array TX-nTHnorm.dat
mv $infile.RX-nTHnorm_array RX-nTHnorm.dat
#cp $infile.TX-nTHnorm_array TX-nTHnorm.dat
#cp $infile.RX-nTHnorm_array RX-nTHnorm.dat
#gzip $infile.TX-nTHnorm_array
#gzip $infile.RX-nTHnorm_array


echo
echo "TX-nTHnorm_array.dat "
echo "head "
echo "ARRAY: Tx_norm[48] "
head -1 TX-nTHnorm.dat
echo

echo
echo "RX-nTHnorm_array.dat "
echo "head "
echo "ARRAY: Rx_norm[200] "
head -1 RX-nTHnorm.dat
echo
echo "-------------------------------------------------------------- "

echo "starting C program now... "
#
# can run this code on S5 only - because of NAG routines!!
#
# Cprogram name: Fred_xcorrel_v3.c
# compile and run C code:
#./makecf_fred_xcorrel_v3
makecf_fred_xcorrel_v3
echo 
echo "-------------------------------------------------------------- "

echo
echo "checking output file from C code: "
echo

echo
echo "xcorrel_TX_RX_array.dat "
echo "head "
echo "ARRAY: Xcorrel_Tx_Rx[200] "
head -1 xcorrel_TX_RX_array.dat
echo

echo
echo "xcorrel_TX_RX_rar_array.dat "
echo "head "
echo "ARRAY: Xcorrel_Tx_Rx_rearranged[200] "
head -1 xcorrel_TX_RX_rar_array.dat
echo

echo
echo "xcorrel_TX_RX_stats.dat "
echo "head "
echo "line no.,std ratio,test stat,xcorrel@lag0,xcorrel_maxval,xcorrel_maxval_binpos,numXcorrelvals_abv0,Xcorrel_mean,Xcorrel_sd,Xcorrel_skew,Xcorrel_kurt "
head -1 xcorrel_TX_RX_stats.dat
echo
echo "-------------------------------------------------------------- "

echo
echo "renaming dat files"
mv xcorrel_TX_RX_array.dat $infile.xcorrel_TX_RX_array
#mv xcorrel_TX_RX_rar_array.dat $infile.xcorrel_TX_RX_rar_array
rm xcorrel_TX_RX_rar_array.dat 
mv xcorrel_TX_RX_stats.dat $infile.xcorrel_TX_RX_stats
#mv TX-nTHnorm.dat $infile.TX-nTHnorm_array
#mv RX-nTHnorm.dat $infile.RX-nTHnorm_array
mv infile_lines.flag $infile.infile_lines.flag
rm TX-nTHnorm.dat 
rm RX-nTHnorm.dat 
echo "-------------------------------------------------------------- "

echo
echo "zipping dat files"
gzip $infile.xcorrel_TX_RX_array
#gzip $infile.xcorrel_TX_RX_rar_array
gzip $infile.xcorrel_TX_RX_stats
#gzip $infile.TX-nTHnorm_array
#gzip $infile.RX-nTHnorm_array
echo "-------------------------------------------------------------- "

echo
echo "moving dat files"
#mv $infile.xcorrel_TX_RX_array /cpnet/xserve3a/slf/GLA_PRO_TMP/L3E/R28/Fred/L3E_FM06_data/TX_RX_xcorrel_arrays
##mv $infile.xcorrel_TX_RX_rar_array.gz /cpnet/xserve3a/slf/GLA_PRO_TMP/L3E/R28/Fred/L3E_FM06_data/TX_RX_xcorrel_rar_arrays
#mv $infile.xcorrel_TX_RX_stats /cpnet/xserve3a/slf/GLA_PRO_TMP/L3E/R28/Fred/L3E_FM06_data/TX_RX_xcorrel_stats
##mv $infile.TX-nTHnorm_array.gz /cpnet/xserve2c/slf/RETRACKER/FRED/L1_FM03_c4_data/TX-nTHnorm_arrays
##mv $infile.RX-nTHnorm_array.gz /cpnet/xserve2c/slf/RETRACKER/FRED/L1_FM03_c4_data/RX-nTHnorm_arrays
mv $infile.xcorrel_TX_RX_array.gz LXX_XXXX_data/TX_RX_xcorrel_arrays
mv $infile.xcorrel_TX_RX_stats.gz LXX_XXXX_data/TX_RX_xcorrel_stats
echo "-------------------------------------------------------------- "

end
echo "-------------------------------------------------------------- "


