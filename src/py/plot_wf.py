# !/usr/bin/env python3
# encoding: utf-8

""" A simple module to extract and plot the waveforms from a GLA01 file
    
"""
from collections.abc import Iterable

from netCDF4 import Dataset
import numpy as np
from scipy.stats import mode
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator

from processing import filterNoise

def plot_tx_wf(filepath, indices=None, output_path=''):
    ''' Function to plot tx_wf
    
        @param filepath: The path to the intermediate netCDF4 file to generate a plot for
        @param indices: The indices of the array to plot. This can be either an integer or
                        list of integers. If a list, the index will be appended to the end
                        of filenames for the output images. If not given, the entire array
                        from the given file will have plots created (this will create a lot
                        of images!)
        @param output_path: The path of the output imgge that will be produced. If left
                            blank, it will output to the local directory with a generated
                            filename
    
    '''
    # TODO - Tx-Rx - Separate functions? Parameter control? Way to detect?
    
    # Get data
    ds = Dataset(filepath, mode='r')
    txwf = ds.variables['tx_wf'][:]
    txwf = filterNoise(txwf)
    
    # Argument processing
    if indices is None:
        indices = range(txwf.shape[0])
    
    if not isinstance(indices, Iterable):
        # Make it a list to work with
        indices = [indices]
        
    if not output_path:
        output_path='txplot.png'
    
    
    for record in indices:
    
        plt.title("Transmit Waveform for Orbit {} Track {}-{}, Index {}".format(ds.OrbitNumber, ds.Track, 
                                                                                ds.Track_Segment, record))
        
        # Mask potential invalidly high values
        #txwf = np.ma.masked_greater(txwf, 1e+10)
        
        plt.minorticks_on()
        plt.xlim(0,48)
        
        plt.gca().xaxis.set_minor_locator( FixedLocator(np.arange(48)) )
        
        plt.xlabel('Bin Number')
        plt.ylabel('Volts (noise floor removed)')
        
        plt.plot(txwf[record])
        
        # Text statistics
        stats = '\n'.join([
            'Record: {}'.format(ds.variables['rec_ndx_gla01'][record]),
            'UTCTime: {:.5f}'.format(ds.variables['UTCTime_40_gla01'][record]),
            'Shot: {}'.format(expandGLA06(filepath, ds.variables['shot_count'][:])[record]),
            'Lat : {:.5f}'.format(expandGLA06(filepath, ds.variables['lat'][:])[record]),
            'Lon: {:.5f}'.format(expandGLA06(filepath, ds.variables['lon'][:])[record]),
            'Energy: {:.5f}'.format(ds.variables['TxNrg_EU'][record]),
            'Mean: {:.5f}'.format(np.mean(txwf[record])),
            'Max: {:.5f}'.format(np.max(txwf[record])),
            'Sigma: {:.5f}'.format(np.std(txwf[record]))
            
            ])
        
        plt.text(.70, .95, stats, transform=plt.gca().transAxes, verticalalignment='top', fontsize='x-small')
        
        # Preserve output path argument
        out_path = output_path
        
        if len(indices) > 1:
            # Append index number to image filename
            if out_path.split('/')[-1].find('.'):
                # Path has an extension
                name, sep, ext = out_path.rpartition('.')
                
                name += '_{}'.format(record)
                
                # Put it back together
                out_path = '.'.join([name, ext])
                
            else:
                out_path += '_{}'.format(record)
                
        # Provide the proper extension if none if none is specified
        if not out_path.split('/')[-1].find('.'):
            out_path += '.png'
        
        plt.savefig(out_path)
        
        # Cleanup
        plt.clf()
    

def plot_rx_wf(filepath, indices=None, output_path=''):
    ''' Function to plot rng_wf
    
        @param filepath: The path to the intermediate netCDF4 file to generate a plot for
        @param indices: The indices of the array to plot. This can be either an integer or
                        list of integers. If a list, the index will be appended to the end
                        of filenames for the output images.If not given, the entire array
                        from the given file will have plots created (this will create a lot
                        of images!)
        @param output_path: The path of the output imgge that will be produced. If left
                            blank, it will output to the local directory with a generated
                            filename
                            
    '''
    # Get data
    ds = Dataset(filepath, mode='r')
    rxwf = ds.variables['rng_wf'][:]
    rxwf = filterNoise(rxwf)
    
    # Argument processing
    if indices is None:
        indices = range(rxwf.shape[0])
    
    if not isinstance(indices, Iterable):
        # Make it a list to work with
        indices = [indices]
        
    if not output_path:
        output_path='rxplot.png'
    
    for record in indices:
    
        plt.title("Receive Waveform for Orbit {} Track {}-{}, Index {}".format(ds.OrbitNumber, ds.Track, 
                                                                                ds.Track_Segment, record))
        
        # Mask potential invalidly high values
        #rxwf = np.ma.masked_greater(rxwf, 1e+10)
        
        plt.minorticks_on()
        #plt.xlim(0,544)
        
        #plt.gca().xaxis.set_minor_locator( FixedLocator(np.arange(544)) )
        
        plt.xlabel('Bin Number')
        plt.ylabel('Volts (noise floor removed)')
        
        plt.plot(rxwf[record])
        
        # Text statistics
        stats = '\n'.join([
            'Record: {}'.format(ds.variables['rec_ndx_gla01'][record]),
            'UTCTime: {:.5f}'.format(ds.variables['UTCTime_40_gla01'][record]),
            'Shot: {}'.format(expandGLA06(filepath, ds.variables['shot_count'][:])[record]),
            'Lat : {:.5f}'.format(expandGLA06(filepath, ds.variables['lat'][:])[record]),
            'Lon: {:.5f}'.format(expandGLA06(filepath, ds.variables['lon'][:])[record]),
            'Energy: {:.5f}'.format(ds.variables['RecNrgAll_EU'][record]),
            'Reflectivity: {:.5f}'.format(expandGLA06(filepath, ds.variables['reflctUC'][:])[record]),
            'Gain: {:.5f}'.format(expandGLA06(filepath, ds.variables['gval_rcv'][:])[record]),
            'Mean: {:.5f}'.format(np.mean(rxwf[record])),
            'Max: {:.5f}'.format(np.max(rxwf[record])),
            'Sigma: {:.5f}'.format(np.std(rxwf[record]))
            
            ])
        
        plt.text(.70, .95, stats, transform=plt.gca().transAxes, verticalalignment='top', fontsize='x-small')
        
        # Preserve output path argument
        out_path = output_path
        
        if len(indices) > 1:
            # Append index number to image filename
            if out_path.split('/')[-1].find('.'):
                # Path has an extension
                name, sep, ext = out_path.rpartition('.')
                
                name += '_{}'.format(record)
                
                # Put it back together
                out_path = '.'.join([name, ext])
                
            else:
                out_path += '_{}'.format(record)
                
        # Provide the proper extension if none if none is specified
        if not out_path.split('/')[-1].find('.'):
            out_path += '.png'
        
        plt.savefig(out_path)
        
        # Cleanup
        plt.clf()
    
    
def expandGLA06(filepath, gla06_array):
    ''' Expand an array from GLA06 to cross-reference with GLA01. Keeps GLA06 elements
        in the proper place of the expanded array according to their record number
    
        @param filepath: The path of the intermediate NetCDF4 file the array comes from
        @param gla06_array: The array to expand to GLA01 size
        
    '''
    
    ds = Dataset(filepath, mode='r')
    
    rec_ndx_01 = ds['rec_ndx_gla01'][:]
    rec_ndx_06 = ds['rec_ndx_gla06'][:]
    
    expanded = np.zeros(rec_ndx_01.shape, dtype=gla06_array.dtype)
    
    expanded[np.isin(rec_ndx_01, rec_ndx_06)] = gla06_array
    
    if np.issubdtype(gla06_array.dtype, np.floating): 
        # Floating point, set fill to nan
        expanded[~np.isin(rec_ndx_01, rec_ndx_06)] = gla06_array.fill_value
        
    else:
        # Set to array's fill value
        expanded[~np.isin(rec_ndx_01, rec_ndx_06)] = gla06_array.fill_value
        
    expanded = np.ma.masked_equal(expanded, gla06_array.fill_value)    
    
    return expanded
    
    
if __name__ == '__main__':
    import sys
    
    plot_tx_wf(sys.argv[1], indices=range(10000,10010))
    plot_rx_wf(sys.argv[1], indices=range(10000,10010))
    
    