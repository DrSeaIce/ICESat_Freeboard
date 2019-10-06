# !/usr/bin/env python3
# encoding: utf-8

""" A module that will handle the processing operations to
    create a Freeboard product
    
"""

import os

import numpy as np
from scipy.stats import skew, kurtosis


def noiseThreshold(wave_array):
    ''' Compute a threshold below which waveform values will 
        be ignored
        
        Returns a noise threshold for each set of bins in the array
        
    '''

    # Compute mean noise using first 10 elements of each bin in TX array
    noise_mean = np.mean(wave_array[:, 0:10], axis=1)
    
    noise_std = np.std(wave_array[:, 0:10], axis=1)
    
    # Use robust standard deviation to determine noise threshold
    noise_threshold = noise_mean + 3*noise_std # Noise threshold
    
    return noise_threshold
        

def filterNoise(wave_array, return_thresh=False):
    ''' Filter out low waveform values using the noise threshold
    
    '''
    
    noise_thresh = noiseThreshold(wave_array)
    
    filtered = np.zeros_like(wave_array)
    
    # Setup rolling window parameters
    num_neighbors = 2 # Number of neighboring elements in each direction to check for being below threshold
    window = num_neighbors*2 + 1 # Size of window to check for noise threshold - equal to either side plus center element
    
    # Pad array bin axis to be able to check starting and ending elements in just one direction
    padded_array = np.pad(wave_array, [ (0,0), (num_neighbors,num_neighbors) ], mode='edge')
    
    # Create a rolling window view of the padded array
    rolling_view = np.lib.stride_tricks.as_strided(padded_array,
                                                   shape=padded_array.shape[:-1] + (padded_array.shape[-1] - window + 1, window), 
                                                   strides=padded_array.strides + (padded_array.strides[-1],),
                                                   writeable=False)
    
    # Create a mask of all windows where all values in window are greater than noise threshold
    mask = np.all(rolling_view > noise_thresh[:,np.newaxis,np.newaxis], axis=rolling_view.ndim-1)
    
    # Set filtered array values to wave array values, subtracting the noise threshold
    filtered[mask] = (wave_array - noise_thresh[:,np.newaxis])[mask]
    
    if return_thresh:
        return filtered, noise_thresh
    
    else:
        return filtered


def normalize(wave_array):
    ''' Normalize the waveform array
    
    '''
    # Compute max for each bin
    maxA = np.max(wave_array, axis=1, keepdims=True)
    
    # Divide each bin by its max (division by zero results in a masked element)
    normed = wave_array / maxA # [:, np.newaxis] - unneeded with keepdims?
    # Set masked entries back to their original values (should be zeroes)
    normed[normed.mask] = wave_array[normed.mask]
    
    ''' Alternative division
    normed = np.zeros_like(wave_array)
    # Divide each bin by its max where the max is not zero, leave rest as zeros
    normed[maxA != 0, :] = wave_array[maxA != 0] / maxA[maxA != 0, np.newaxis]
    '''
    
    return normed


def getStats(wave_array):
    ''' Compute the Mean, Standard deviation, Skewness, and Kurtosis
        on a wave array
        
    '''
    
    maxVal = np.max(wave_array, axis=1)
    maxPos = np.argmax(wave_array, axis=1)
    mean = np.mean(wave_array, axis=1)
    std_dev = np.std(wave_array, axis=1)
    skewness = skew(wave_array, axis=1)
    kurt = kurtosis(wave_array, axis=1)
    
    # TODO - return dict?
    return maxVal, maxPos, mean, std_dev, skewness, kurt


def writeStats(filepath, **variables):
    ''' Write the passed keyword parameters to the intermediate NetCDF4 file
    
        Variables passed in should be arrays matching either the transmit (tx_wf)
        or receive (rng_wf) waveform array dimensions
    
    '''
    from netCDF4 import Dataset
    
    nc_file = Dataset(filepath, mode='a')
    
    # Transmit and receive waveform array shapes
    tx_shape = nc_file['tx_wf'].shape
    rx_shape = nc_file['rng_wf'].shape
    
    # For each variable parameter, add it to the file using the correct dimension
    for varname, vararray in variables.items():
        # Check for pre-existence of variable names in NetCDF4 file
        if varname in nc_file.variables:
            # Overwrite it
            nc_var = nc_file.variables[varname]
            
        else:
            # Create it
            if vararray.shape == tx_shape:
                # Transmit statistics
                nc_var = nc_file.createVariable(varname, vararray.dtype, ('UTCTime_GLA01', 'tx_bin'))
                
            elif vararray.shape == rx_shape:
                # Receive statistics
                nc_var = nc_file.createVariable(varname, vararray.dtype, ('UTCTime_GLA01', 'rx_bin'))
            
            else:
                # 1-D statistic
                nc_var = nc_file.createVariable(varname, vararray.dtype, 'UTCTime_GLA01')
            
        
            
        ''' Unneeded with some stats being of 1 dimension   
        else:
            raise ValueError("Keyword argument value is not an array of the same shape as " +
                             "either the transmit or receive waveform array: {}".format(varname))
        '''
        nc_var[:] = vararray
        
    nc_file.close()
    
    
def xcorr_pad(wave_array):
    ''' Pad a wave array to a consistent value for cross-correlation of Rx and Tx arrays
        Center the peak within padded bins
    
    '''
    TARGET_BIN_SIZE = 1090
    
    # Pad array to the target bin size
    # Padded elements are zero
    padded_array = np.pad(wave_array, [ (0,0), (int(np.ceil((TARGET_BIN_SIZE - wave_array.shape[1]) / 2)),
                                                int(np.floor((TARGET_BIN_SIZE - wave_array.shape[1]) / 2))) ], 
                          mode='constant', constant_values=0 )

    # Find peak and position of peak within bins
    peaks = np.max(padded_array, axis=1, keepdims=True)
    peaks_index = np.argmax(padded_array, axis=1)
    
    # Center the peak
    # More complicated roll in order to move each bin by a different amount
    
    # Determine the shift amount for each bin
    shifts = int(np.ceil(padded_array.shape[1]/2)) - peaks_index
    # Make all shifts positive
    shifts[shifts < 0] = shifts[shifts < 0] + padded_array.shape[1]
    # Array of indices to the padded array 
    padded_indices = np.ogrid[:padded_array.shape[0], :padded_array.shape[1]]
    # Shift the indices by he intended amount
    padded_indices[1] = padded_indices[1] - shifts[:, np.newaxis]
    # Create the resulting array from indexing the original
    padded_centered = padded_array[tuple(padded_indices)]
    
    #padded_array = np.roll(padded_array, int(np.ceil(padded_array.shape[1]/2)) - peaks_index, axis=1)
    
    return padded_centered


def xcorr(tx_array, rx_array):
    ''' Perform cross correlation between transmit and receive waveforms
    
    '''
    from scipy.signal import correlate
    
    # Prep arrays by padding them
    padded1 = xcorr_pad(tx_array)
    padded2 = xcorr_pad(rx_array)
    
    # Perform the cross-correlation and return result
    result = correlate(padded1, padded2, mode='same')
    
    # Reverse array (deactivated until known to be necessary
    #result = np.flip(result)
    
    return result


if __name__ == '__main__':
    from netCDF4 import Dataset
    
    ds = Dataset('intermediate_test.nc', mode='r')
    
    tx_wave = ds['tx_wf'][:]
    rx_wave = ds['rng_wf'][:]
    
    ds.close()
    
    # Transmit wave
    tx_filtered, tx_noise = filterNoise(tx_wave, return_thresh=True)
    print("Tx Noise filtered")
    tx_norm = normalize(tx_filtered)
    print("Tx normalized")
    tx_max, tx_maxPos, tx_mean, tx_std, tx_skew, tx_kurt = getStats(tx_norm)
    print("Transmit stats processed")
    
    # Receive wave
    rx_filtered, rx_noise = filterNoise(rx_wave, return_thresh=True)
    print("Rx noise filtered")
    rx_norm = normalize(rx_filtered)
    print("Rx normalized")
    rx_max, rx_maxPos, rx_mean, rx_std, rx_skew, rx_kurt = getStats(rx_norm)
    print("Rx stats processed")
    
    writeStats('intermediate_test.nc', 
               tx_norm=tx_norm, rx_norm=rx_norm, tx_noise_thresh=tx_noise, rx_noise_thresh=rx_noise,
               rx_max=rx_max, rx_maxPos=rx_maxPos,
               rx_mean=rx_mean, rx_std_dev=rx_std, rx_skewness=rx_skew, rx_kurtosis=rx_kurt,
               tx_max=tx_max, tx_maxPos=tx_maxPos, 
               tx_mean=tx_mean, tx_std_dev=tx_std, tx_skewness=tx_skew, tx_kurtosis=tx_kurt)
    print("Stats written")
    
    x_corr = xcorr(tx_norm, rx_norm)
    print("Xcorrel calculated")
    xcorr_max, xcorr_maxPos, xcorr_mean, xcorr_std, xcorr_skew, xcorr_kurt = getStats(x_corr)
    print("Xcorrel stats processed")
    writeStats('intermediate_test.nc', 
               xcorr_max=xcorr_max, xcorr_maxPos=xcorr_maxPos, 
               xcorr_mean=xcorr_mean, xcorr_std=xcorr_std, xcorr_skew=xcorr_skew, xcorr_kurt=xcorr_kurt)
    print("Xcorrel stats written")
    print("Finished")
    