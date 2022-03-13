# TT-scna

## Usage:

## import ttscna
## ttscna.ttscna(file_name, on_start, off_start, sweep_len = 5.0, sample_rate = 20000, skip_rows = 9, sep = "\t", fmin = 20.0, fmax = 10000.0)


Takes as input a text file containing electrophysiological data in two conditions.
Computes the difference in single-channel noise between the two conditions.

Must contain one column, with time, expressed in seconds,
and a second column, with current, expressed in picoamps.

Inputs
----------
file_name: str
    Path and File_name to text data file.

on_start: float
    Start time of 'on' condition.

off_start: float
    Start time of 'off' condition.

sweep_len: float
    Length of time, in seconds, of each condition to analyze.

sample_rate: float
    Sampling rate of the time series, in Hz.

skip_rows: int
    Rows to skip in the header of the data file.

sep: str
    Delimiter used in each column for the data file.

fmin: float
    Minimum frequency to fit to, in Hz.

fmax: float
    Maximum frequency to fit to, in Hz.

Output
----------
Figure: .png
    A 4-panel plot containing:
        
    100 microsecond traces of the On and Off conditions.
    
    Overlaid power spectra from both conditions.
        
    The power spectra of the on condition minus the off condition, 
    overlaid with a Lorentzian model.
        
Table: .csv
    A .csv containing:

    Unitary conductance, in fS.
    
    Mean Dwell Time, in ms.
        
    Mean Bulk Current, in pA.
