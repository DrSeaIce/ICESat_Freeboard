# !/usr/bin/env python3
# encoding: utf-8

""" A module to read the IceSat GLAH files and extract data needed for Freeboard processing into an intermediate file
    
"""
import os

from netCDF4 import Dataset
import h5py
import numpy as np

def readGLA01File(filepath):
    ''' Initial function. This may be expanded or scrapped based on need for
        the algorithm and associated software
        
    '''

    # Determine file type from extension
    fname,dot,ext = filepath.rpartition('.')
    
    datadict = {}
    
    if ext.lower() == 'h5':
        h5f = h5py.File(filepath, mode='r')
        
        # Read everything in the file and store it in a dictionary
        datadict['rec_ndx_gla01'] = np.array( h5f.get('/Data_40HZ/Time/i_rec_ndx') )
        datadict['UTCTime_40_gla01'] = np.array( h5f.get('/Data_40HZ/DS_UTCTime_40') )
        #datadict['rectype'] = None
        datadict['lat_pred'] = np.array( h5f.get('/Data_1HZ/Geolocation/d1_pred_lat') ) # 1HZ variable
        datadict['lat_pred'] = expand40Hz(datadict['lat_pred'])
        datadict['lon_pred'] = np.array( h5f.get('/Data_1HZ/Geolocation/d1_pred_lon') ) # 1HZ variable
        datadict['lon_pred'] = expand40Hz(datadict['lon_pred'])
        datadict['RecNrgAll_EU'] = np.array( h5f.get('/Data_40HZ/Waveform/Characteristics/d_RecNrgAll_EU') )
        datadict['TxNrg_EU'] = np.array( h5f.get('/Data_1HZ/Transmit_Energy/d_TxNrg_EU') )
        datadict['TxNrg_EU'] = expand40Hz(datadict['TxNrg_EU'])
        datadict['gainset1064'] = np.array( h5f.get('/Data_40HZ/Waveform/Characteristics/i_gainSet1064') )
        #datadict['echo'] = None
        datadict['time_txWfPk'] = np.array( h5f.get('/Data_40HZ/Waveform/TransmitWaveform/i_time_txWfPk') )
        datadict['txWfPk_Flag'] = np.array (h5f.get('/Data_40HZ/Waveform/TransmitWaveform/i_txWfPk_Flag') )
        datadict['TxWfStart'] = np.array( h5f.get('/Data_40HZ/Waveform/TransmitWaveform/i_TxWfStart') )
        datadict['tx_wf'] = np.array( h5f.get('/Data_40HZ/Waveform/TransmitWaveform/r_tx_wf') )
        datadict['RespEndTime'] = np.array( h5f.get('/Data_40HZ/Waveform/RecWaveform/i_RespEndTime') )
        datadict['rec_wf_location_index'] = np.array( h5f.get('/Data_40HZ/Waveform/RecWaveform/i_rec_wf_location_index') )
        datadict['rng_wf'] = np.array( h5f.get('/Data_40HZ/Waveform/RecWaveform/r_rng_wf') )
        
        # Cleanup
        h5f.close()
        
    # May as well implement netCDF4 reading while we're at it
    elif ext.lower() == 'nc':
        nc_file = Dataset(filepath, 'r')
        
        # TODO - Attribute/variable reading
        
        nc_file.close()
    
    return datadict
    
    
def readGLA05File(filepath):
    ''' Function to read GLAH05. GLA reading functions may be merged later if a good way to handle the naming 
        differences is found
    
    ''' 
    
    # Determine file type from extension
    fname,dot,ext = filepath.rpartition('.')
    
    datadict = {}
    
    if ext.lower() == 'h5':
        h5f = h5py.File(filepath, mode='r')
        
        # Read everything in the file and store it in a dictionary
        datadict['wfFitSDev_2'] = np.array( h5f.get('/Data_40HZ/Waveform/d_wfFitSDev_2') )
        datadict['RecNrgAll_gla05'] = np.array( h5f.get('/Data_40HZ/Waveform/d_RecNrgAll') )
        
        # Cleanup
        h5f.close()
        
    # May as well implement netCDF4 reading while we're at it
    elif ext.lower() == 'nc':
        nc_file = Dataset(filepath, 'r')
        
        # TODO - Attribute/variable reading
        
        nc_file.close()
    
    return datadict


def readGLA06File(filepath):
    ''' Function to read GLAH06. GLA reading functions may be merged later if a good way to handle the naming 
        differences is found
    
    '''
    
    # Determine file type from extension
    fname,dot,ext = filepath.rpartition('.')
    
    datadict = {}
    
    if ext.lower() == 'h5':
        h5f = h5py.File(filepath, mode='r')
        
        # Read everything in the file and store it in a dictionary
        datadict['rec_ndx_gla06'] = np.array( h5f.get('/Data_40HZ/Time/i_rec_ndx') )
        datadict['UTCTime_40_gla06'] = np.array( h5f.get('/Data_40HZ/DS_UTCTime_40') )
        datadict['lat'] = np.array( h5f.get('/Data_40HZ/Geolocation/d_lat') )
        datadict['lon'] = np.array( h5f.get('/Data_40HZ/Geolocation/d_lon') )
        datadict['elev'] = np.array( h5f.get('/Data_40HZ/Elevation_Surfaces/d_elev') )
        datadict['SolAng'] = np.array( h5f.get('/Data_1HZ/Reflectivity/d_SolAng') ) # 1HZ variable
        datadict['SolAng'] = expand40Hz(datadict['SolAng'])
        datadict['gdHt'] = np.array( h5f.get('/Data_40HZ/Geophysical/d_gdHt') )
        datadict['erElv'] = np.array( h5f.get('/Data_40HZ/Geophysical/d_erElv') )
        #datadict['spElv'] = None
        datadict['poTide'] = np.array( h5f.get('/Data_40HZ/Geophysical/d_poTide') )
        datadict['ldElv'] = np.array( h5f.get('/Data_40HZ/Geophysical/d_ldElv') )
        datadict['ocElv'] = np.array( h5f.get('/Data_40HZ/Geophysical/d_ocElv') )
        datadict['wTrop'] = np.array( h5f.get('/Data_40HZ/Elevation_Corrections/d_wTrop') )
        datadict['dTrop'] = np.array( h5f.get('/Data_40HZ/Elevation_Corrections/d_dTrop') )
        #datadict['surftype'] = None
        datadict['surf_ld_flg'] = np.array( h5f.get('/Data_1HZ/Elevation_Flags/surf_ld_flg') ) # 1HZ variable
        datadict['surf_si_flg'] = np.array( h5f.get('/Data_1HZ/Elevation_Flags/surf_si_flg') ) # 1HZ variable
        datadict['surf_oc_flg'] = np.array( h5f.get('/Data_1HZ/Elevation_Flags/surf_oc_flg') ) # 1HZ variable
        datadict['surf_is_flg'] = np.array( h5f.get('/Data_1HZ/Elevation_Flags/surf_is_flg') ) # 1HZ variable
        datadict['surf_ld_flg'] = expand40Hz(datadict['surf_ld_flg'])
        datadict['surf_si_flg'] = expand40Hz(datadict['surf_si_flg'])
        datadict['surf_oc_flg'] = expand40Hz(datadict['surf_oc_flg'])
        datadict['surf_is_flg'] = expand40Hz(datadict['surf_is_flg'])
        datadict['DEM_elv'] = np.array( h5f.get('/Data_40HZ/Geophysical/d_DEM_elv') )
        datadict['refRng'] = np.array( h5f.get('/Data_40HZ/Elevation_Surfaces/d_refRng') )
        datadict['TrshRngOff'] = np.array( h5f.get('/Data_40HZ/Elevation_Offsets/d_TrshRngOff') )
        datadict['SigBegOff'] = np.array( h5f.get('/Data_40HZ/Elevation_Offsets/d_SigBegOff') )
        datadict['SigEndOff'] = np.array( h5f.get('/Data_40HZ/Elevation_Offsets/d_SigEndOff') )
        datadict['cntRngOff'] = np.array( h5f.get('/Data_40HZ/Elevation_Offsets/d_cntRngOff') )
        datadict['reflctUC'] = np.array( h5f.get('/Data_40HZ/Reflectivity/d_reflctUC') )
        datadict['reflCor_atm'] = np.array( h5f.get('/Data_1HZ/Reflectivity/d_reflCor_atm') ) # 1HZ variable
        datadict['reflCor_atm'] = expand40Hz(datadict['reflCor_atm'])
        #datadict['SigmaElv'] = None
        datadict['kurt2'] = np.array( h5f.get('/Data_40HZ/Waveform/d_kurt2') )
        datadict['skew2'] = np.array( h5f.get('/Data_40HZ/Waveform/d_skew2') )
        #datadict['srf_ruf'] = None
        datadict['isRngOff'] = np.array( h5f.get('/Data_40HZ/Elevation_Offsets/d_isRngOff') )
        datadict['siRngOff'] = np.array( h5f.get('/Data_40HZ/Elevation_Offsets/d_siRngOff') )
        datadict['ldRngOff'] = np.array( h5f.get('/Data_40HZ/Elevation_Offsets/d_ldRngOff') )
        datadict['ocRngOff'] = np.array( h5f.get('/Data_40HZ/Elevation_Offsets/d_ocRngOff') )
        datadict['nPeaks1'] = np.array( h5f.get('/Data_40HZ/Waveform/i_nPeaks1') )
        #datadict['srf_slope'] = None
        #datadict['PODFixedPos'] = None
        datadict['cldl_mswf_flg'] = np.array( h5f.get('/Data_1HZ/Atmosphere/cld1_mswf_flg') ) # 1HZ variable
        datadict['cldl_mswf_flg'] = expand40Hz(datadict['cldl_mswf_flg'])
        datadict['satNdx'] = np.array( h5f.get('/Data_40HZ/Quality/i_satNdx') )
        datadict['satElevCorr'] = np.array( h5f.get('/Data_40HZ/Elevation_Corrections/d_satElevCorr') )
        datadict['sat_corr_flg']  = np.array( h5f.get('/Data_40HZ/Quality/sat_corr_flg') )
        datadict['satNrgCorr'] = np.array( h5f.get('/Data_40HZ/Reflectivity/d_satNrgCorr') )
        #datadict['satPwdCorr'] = None
        datadict['gval_rcv'] = np.array( h5f.get('/Data_40HZ/Waveform/i_gval_rcv') )
        datadict['RecNrgAll_gla06'] = np.array( h5f.get('/Data_40HZ/Reflectivity/d_RecNrgAll') )
        #datadict['msRngCorr'] = None
        #datadict['msCorrFlg'] = None
        datadict['Surface_temp'] = np.array( h5f.get('/Data_1HZ/Atmosphere/d_Surface_temp') ) # 1HZ variable
        datadict['Surface_pres'] = np.array( h5f.get('/Data_1HZ/Atmosphere/d_Surface_pres') ) # 1HZ variable
        datadict['Surface_relh'] = np.array( h5f.get('/Data_1HZ/Atmosphere/d_Surface_relh') ) # 1HZ variable
        datadict['Surface_temp'] = expand40Hz(datadict['Surface_temp'])
        datadict['Surface_pres'] = expand40Hz(datadict['Surface_pres'])
        datadict['Surface_relh'] = expand40Hz(datadict['Surface_relh'])
        datadict['skew2'] = np.array( h5f.get('/Data_40HZ/Waveform/d_skew2') )
        datadict['kurt2'] = np.array( h5f.get('/Data_40HZ/Waveform/d_kurt2') )
        datadict['shot_count'] = np.array( h5f.get('/Data_40HZ/Time/i_shot_count') )
        datadict['elevBiasCorr'] = np.array( h5f.get('/Data_40HZ/Geophysical/d_ElevBiasCorr') )
        datadict['GmC'] = np.array( h5f.get('/Data_40HZ/Elevation_Corrections/d_GmC') )
        
        
        # Cleanup
        h5f.close()
        
    # May as well implement netCDF4 reading while we're at it
    elif ext.lower() == 'nc':
        nc_file = Dataset(filepath, 'r')
        
        # TODO - Attribute/variable reading
        
        nc_file.close()
    
    return datadict


def expand40Hz(one_hz_array):
    ''' Expand 1Hz variables to match the rest at 40Hz
    
    '''
    # For each element of the array, repeat it 40 times in the new array
    forty_hz_array = one_hz_array.take([x for x in range(one_hz_array.size) for n in range(40)])
    
    return forty_hz_array


def readGLAHeader(filepath):
    ''' Function to return a structure for GLA06's header attributes
        tobe included in the intermediate file
        
    '''
    # Determine file type from extension
    fname,dot,ext = filepath.rpartition('.')
    
    attr_dict = {}
    
    if ext.lower() == 'h5':
        h5f = h5py.File(filepath, mode='r')
        
        # Read the attribute information into a dictionary
        attr_dict = dict(h5f.attrs)
        
        # Decode ASCII byte strings into unicode
        for key, val in attr_dict.items():
            if isinstance(val, bytes):
                try:
                    attr_dict[key] = val.decode()
                    
                except AttributeError:
                    print("Bytes attribute {} can't be decoded.".format( str(key) ))
    
    return attr_dict


def writeIntermediateFile(filepath, datadict, gla6_dict, attributes=None):
    ''' Function to write a data dictionary out to a simple NetCDF4 file
    
    '''
    
    ncfile = Dataset(filepath, 'w')
    
    # Write header attributes if included
    if attributes:
        ncfile.setncatts(attributes)
    
    # Create dimensions for UTCTime
    ncfile.createDimension('UTCTime_GLA01', datadict['UTCTime_40_gla01'].size)
    # GLA(01/05) have a different UTCTime array size than GLA06
    ncfile.createDimension('UTCTime_GLA06', gla6_dict['UTCTime_40_gla06'].size)
    # Create dimension "tx_bin" for tx_wf variable
    ncfile.createDimension('tx_bin', datadict['tx_wf'].shape[1])
    # Create dimention "rx_bin" for rng_wf variable
    ncfile.createDimension('rx_bin', datadict['rng_wf'].shape[1])
    
    # Special case for two-dimensional variable tx_wf
    ds_tx_wf = ncfile.createVariable('tx_wf', datadict['tx_wf'].dtype, ('UTCTime_GLA01','tx_bin'))
    ds_tx_wf[:] = datadict['tx_wf']
    
    # Special case for two-dimensional variable rng_wf
    ds_rng_wf = ncfile.createVariable('rng_wf', datadict['rng_wf'].dtype, ('UTCTime_GLA01', 'rx_bin'))
    ds_rng_wf[:] = datadict['rng_wf']
    
    # Then delete special cases from the dictionary so it doesn't duplicate
    del datadict['tx_wf']
    del datadict['rng_wf']
    
    '''
    del datadict['UTCTime_40_gla06']
    del datadict['UTCTime_40_gla01']'''
    
    # Write GLA01 and GLA05 variables
    for varname, vardata in datadict.items():
        
        try:
            print(varname)
            ds_var = ncfile.createVariable(varname, vardata.dtype, 'UTCTime_GLA01')
            
        except TypeError:
            raise TypeError('Unsupported numpy dtype: {}'.format(vardata.dtype.str))
        
        ds_var[:] = vardata
    
    # Write GLA06 variables
    for varname, vardata in gla6_dict.items():
        
        try:
            print(varname)
            ds_var = ncfile.createVariable(varname, vardata.dtype, 'UTCTime_GLA06')
            
        except TypeError:
            raise TypeError('Unsupported numpy dtype: {}'.format(vardata.dtype.str))
        
        ds_var[:] = vardata
        
    ncfile.close()
    
    return filepath
    

def readGLA(track, segment, folder_path):
    
    gla1_filename = 'GLAH01_033_1102_004_{}_{}_02_0001.H5'.format(track, segment)
    gla5_filename = 'GLAH05_634_1102_004_{}_{}_01_0001.H5'.format(track, segment)
    gla6_filename = 'GLAH06_634_1102_004_{}_{}_01_0001.H5'.format(track, segment)
    
    gla1_dict = readGLA01File('/'.join([folder_path, gla1_filename]))
    gla5_dict = readGLA05File('/'.join([folder_path, gla5_filename]))
    gla6_dict = readGLA06File('/'.join([folder_path, gla6_filename]))
    
    # Combine GLA01 and GLA05 dictionaries to write to file
    # Keep GLA06 data separate as it has different dimensions
    datadict = dict( **gla1_dict, **gla5_dict)
    
    # Mask invalid values
    for key, val in datadict.items():
        if np.issubdtype(val.dtype, np.floating):
            datadict[key] = np.ma.masked_equal(val, np.finfo(val.dtype).max)
            
        elif np.issubdtype(val.dtype, np.integer):
            datadict[key] = np.ma.masked_equal(val, np.iinfo(val.dtype).max)
            
        else:
            # Needed?
            datadict[key] = np.ma.masked_equal(val, val.fill_value)
            
    # Repeat masking for the GLA06 dictionary
    for key, val in gla6_dict.items():
        if np.issubdtype(val.dtype, np.floating):
            gla6_dict[key] = np.ma.masked_equal(val, np.finfo(val.dtype).max)
            
        elif np.issubdtype(val.dtype, np.integer):
            gla6_dict[key] = np.ma.masked_equal(val, np.iinfo(val.dtype).max)
            
        else:
            # Needed?
            gla6_dict[key] = np.ma.masked_equal(val, val.fill_value)
    
    # Read the GLA06 attributes
    gla06_header = readGLAHeader('/'.join([folder_path, gla6_filename]))
    
    file = writeIntermediateFile('intermediate_test.nc', datadict, gla6_dict, gla06_header)
    
    print('Wrote intermediate file: {}'.format(file))
    

if __name__ == '__main__':
    import sys
    
    readGLA(sys.argv[1], sys.argv[2], os.getcwd())
    