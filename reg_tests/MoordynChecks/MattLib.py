#!/usr/bin/env python

# should combine this with other library I made at some point....
# Matt Hall (mtjhall@alumni.uvic.ca)
# updated April 1, 2015
# updated 2015-09-2X with smarter starting line detection features in read_output_file
# updated 2015-11-10 with RotMat for creating rotation matrix from Euler angles (need to check convention)
# updated 2015-11-19 with handling of extra data at ends of some rows
# updated 2016-07-30 with option to read FAST output files in binary form
# updated 2016-08-15 to add wave number iterative calculation function
# updated 2017-11    to add mass matrix calculation (I think it's right)
# updated 2018-04-17 to add utility functions for drawing mooring lines
# updated 2018-08-19 to add flexibility for loading csv files with missing entries for energy systems stuff
# updated 2020-02-24 to add function for mooring line 3D plotting (as was shared by Alex of Saitec)
# 2021-02-10: adding fix to load data files that have the weird fortran E missing: e.g. 3.23232-108
# 2021-05-12: adding dictionary output option for read_output_file

import os
import sys
import math
import re

#import matplotlib
import numpy as np

import matplotlib.colors as colors
import matplotlib.cm as cmx
import struct

import matplotlib.pyplot as plt
#from scipy import signal # for pwelch function (psd functionality)
import matplotlib.mlab as mlab       # for mlab's psd functionality
#import matplotlib.gridspec as gridspec

import shutil as shu
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

        
# ==============================================

def read_output_file(dirName,fileName, skiplines=-1, hasunits=1, chanlim=999, dictionary=True):

    # load data from FAST output file
    # looks for channel names, then units (if hasunits==1), then data lines after first skipping [skiplines] lines.
    # skiplines == -1 signals to search for first channel names line based on starting channel "Time".
    
#   print('attempting to load '+dirName+fileName)
    f = open(dirName+fileName, 'r')
    
    channels = []
    units = []
    data = []
    i=0
    
    for line in f:          # loop through lines in file
    
        if (skiplines == -1):               # special case signalling to search for "Time" at start of channel line
            entries = line.split()          # split elements by whitespace
            print(entries)
            if entries[0].count('Time') > 0 or entries[0].count('time') > 0:  # if we find the time keyword
                skiplines = i
                print("got skiplines="+str(i))
            else:
                pass
    
        if (i < skiplines or skiplines < 0):        # if we haven't gotten to the first channel line or we're in search mode, skip
            pass
            
        elif (i == skiplines):
            for entry in line.split():      # loop over the elemets, split by whitespace
                channels.append(entry)      # append to the last element of the list
                
        elif (i == skiplines+1 and hasunits == 1):
            for entry in line.split():      # loop over the elemets, split by whitespace
                if entry.count('kN') > 0 and entry.count('m') > 0:  # correct for a possible weird character
                    entry = '(kN-m)'
                    
                units.append(entry)         # append to the last element of the list
        
        elif len(line.split()) > 0:
            data.append([])  # add a new sublist to the data matrix
            
            r = re.compile(r"(?<=\d)\-(?=\d)")  # catch any instances where a large negative exponent has been written with the "E"
            line2 = r.sub("E-",line)           # and add in the E
            
            j=0
            for entry in line2.split():      # loop over the elements, split by whitespace
                if j > chanlim:
                    break
                j+=1    
                data[-1].append(entry)      # append to the last element of the list
    
        else:
            break
    
        i+=1
    
    f.close()  # close data file
    
    
    # use a dictionary for convenient access of channel columns (eg. data[t][ch['PtfmPitch'] )
    ch = dict(zip(channels, range(len(channels))))
    
    #print ch['WindVxi']
    data2 = np.array(data)
    
    data3 = data2.astype(float)
    
    if dictionary:
        dataDict = {}
        unitDict = {}
        for i in range(len(channels)):
            dataDict[channels[i]] = data3[:,i]
            unitDict[channels[i]] = units[i]
        return dataDict, unitDict
    else:
        return data3, ch, channels, units
    

    
def read_csv_file(dirName, fileName, skiplines=-1, hasunits=1, sep=",", chanlim=999):

    # load data from a CSV file that may have missing entries (what to do with them in that case?)
    
    f = open(dirName+fileName, 'r')
    
    channels = []
    units = []
    data = []
    i=0
    
    for line in f:          # loop through lines in file
    
        if (skiplines == -1):               # special case signalling to search for "Time" at start of channel line
            entries = line.split(sep)
            if entries[0].count('Time') > 0 or entries[0].count('time') > 0:  # if we find the time keyword
                skiplines = i
                print("got skiplines="+str(i))
            else:
                pass
    
        if (i < skiplines or skiplines < 0):        # if we haven't gotten to the first channel line or we're in search mode, skip
            pass
            
        elif (i == skiplines):
            for entry in line.split(sep):       # loop over the elemets, split 
                channels.append(entry.lstrip())         # append to the last element of the list
                
        elif (i == skiplines+1 and hasunits == 1):
            for entry in line.split(sep):       # loop over the elemets, split 
                units.append(entry)         # append to the last element of the list
        
        elif len(line.split(sep)) > 0:
            data.append([])  # add a new sublist to the data matrix
            
            j=0
            for entry in line.split(sep):       # loop over the elements, split
                if j > chanlim:
                    break
                j+=1    
                data[-1].append(entry)      # append to the last element of the list
    
        else:
            break
    
        i+=1
    
    f.close()  # close data file
    
    
    # use a dictionary for convenient access of channel columns (eg. data[t][ch['PtfmPitch'] )
    ch = dict(zip(channels, range(len(channels))))
    
    #print ch['WindVxi']

    data2 = np.array(data)  # make a numpy array of the read in data (string form for now)

    # now check for funny data / columns that can't be converted to floats
    print(data2.shape)
    badDataCounter = 0
    for i in range(len(ch)):        
        for j in range(data2.shape[0]):
        
            try:
                np.float(data2[j,i])
                
            except ValueError:
                #print("Detected non-floatable entry in channel "+channels[i]+" line "+str(j)+" so replacing with -1")
                #print(data2[j,i])
                if j>1:
                    data2[j,i] = data2[j-1,i]  # use previous value in channel
                else:
                    data2[j,i] = "-1"
                    
                badDataCounter += 1
                
    if badDataCounter>0:
        print("Warning: "+str(badDataCounter)+" instances of unrecognized numbers when loading "+fileName)
        
        #try:
        #   data2[:,i].astype(float)
        #   
        #except ValueError:
        #   print('\n')
        #   print("Detected non-floatable entries in channel "+channels[i]+" so replacing with -1")
        #   print(data2[:,i])
        #   
        #   #data2[:,i] = -1  # and replace them with -1s if necessary
        #   
        #   return data2[:,i], ch, channels, units
        #   
        #   #break
    
    data3 = data2.astype(float)
    
    return data3, ch, channels, units
    

    
    
def read_binary_output_file(dirName,fileName, LenName=10, LenUnit=10):

    # load data from FAST binary output file
    # LenName and LenUnit specify number of characters expected for channel names and units
    
    print('attempting to load '+dirName+fileName)
    
    #try:
    fid = open(dirName+fileName, 'rb')
    
    #print("fid is "+str(fid)+" I think")
    
    #------------------------  get the header information ----------------------------

    tempdata = fid.read(2)   # do this to avoid errors if the file is empty
    
    if len(tempdata) < 2:
        print("Unable to read data from file.")
        return [], [], [], []
    
    FileID       = struct.unpack('h', tempdata)[0]  # fread( fid, 1, 'int16');             # FAST output file format, INT(2)
    #FileID       = struct.unpack('h', fid.read(2))[0]  # fread( fid, 1, 'int16');             # FAST output file format, INT(2)
        
    NumOutChans  = struct.unpack('i', fid.read(4))[0]  # fread( fid, 1, 'int32');             # The number of output channels, INT(4)
    NT           = struct.unpack('i', fid.read(4))[0]  # fread( fid, 1, 'int32');             # The number of time steps, INT(4)

    if FileID == 1:
        TimeScl  = struct.unpack('d', fid.read(8))[0]  # fread( fid, 1, 'float64');           # The time slopes for scaling, REAL(8)
        TimeOff  = struct.unpack('d', fid.read(8))[0]  # fread( fid, 1, 'float64');           # The time offsets for scaling, REAL(8)
    else:
        TimeOut1 = struct.unpack('d', fid.read(8))[0]  # fread( fid, 1, 'float64');           # The first time in the time series, REAL(8)
        TimeIncr = struct.unpack('d', fid.read(8))[0]  # fread( fid, 1, 'float64');           # The time increment, REAL(8)

    ColScl       = np.fromfile(fid, count=NumOutChans, dtype=np.float32)  # fread( fid, NumOutChans, 'float32'); # The channel slopes for scaling, REAL(4)
    ColOff       = np.fromfile(fid, count=NumOutChans, dtype=np.float32)  # fread( fid, NumOutChans, 'float32'); # The channel offsets for scaling, REAL(4)

    LenDesc      = struct.unpack('i', fid.read(4))[0]  # fread( fid, 1,           'int32' );  # The number of characters in the description string, INT(4)

    DescStrASCII = fid.read(LenDesc)                # fread( fid, LenDesc,     'uint8' );  # DescStr converted to ASCII
    DescStr      = str( DescStrASCII)                     


    channels = []  #cell(NumOutChans+1,1);                   # initialize the ChanName cell array
    for iChan in range(NumOutChans+1):             #  1:NumOutChans+1 
        ChanNameASCII = fid.read(LenName)                     # fread( fid, LenName, 'uint8' ); # ChanName converted to numeric ASCII
        channels.append((str(ChanNameASCII).strip()[2:12]).strip())             #    = strtrim( char(ChanNameASCII') );

    units = [] #cell(NumOutChans+1,1);                   # initialize the ChanUnit cell array
    for iChan in range(NumOutChans+1):             #1:NumOutChans+1
        ChanUnitASCII = fid.read(LenUnit)          #fread( fid, LenUnit, 'uint8' ); # ChanUnit converted to numeric ASCII
        units.append(str(ChanUnitASCII).strip())          # {iChan}= strtrim( char(ChanUnitASCII') );
    
    #   if units[-1].count('kN') > 0 and units[-1].count('m') > 0:  # correct for a possible weird character
    #       units[-1] = '(kN-m)'
        
#   print('Reading from the file '+input_file+' with heading: ')
#   print(DescStr)

    #-------------------------  get the channel time series -------------------------

    nPts        = NT*NumOutChans;           # number of data points in the file   
    data        = np.zeros((NT,NumOutChans+1));  # output channels (including time in column 1)

    if FileID == 1:
        PackedTime = np.fromfile(fid, count=NT, dtype=np.int32)

        #[PackedTime cnt] = fread( fid, NT, 'int32' ); # read the time data
        #if ( cnt < NT )  :
        #   fid.close()
        #   print('Could not read entire file: read '+str( cnt )+' of '+str( NT )+' time values.')
        

    PackedData = np.fromfile(fid, count=nPts, dtype=np.int16)   
    #[PackedData cnt] = fread( fid, nPts, 'int16' ); # read the channel data
    #if ( cnt < nPts ) 
    #   fclose(fid);
    #   error(['Could not read entire ' FileName ' file: read ' num2str( cnt ) ' of ' num2str( nPts ) ' values.']);
    #end

    fid.close()

    #-------------------------
    # Scale the packed binary to real data
    #-------------------------

    ip=0                            #ip = 1;
    for it in range(NT):                            #for it = 1:NT
        for ic in range(NumOutChans):                       #   for ic = 1:NumOutChans
            data[it, ic+1] = ( PackedData[ip] - ColOff[ic] ) / ColScl[ic]   #       Channels(it,ic+1) = ( PackedData(ip) - ColOff(ic) ) / ColScl(ic) ;
            ip = ip + 1;
                                #   end # ic       
                                #end #it

    if FileID == 1:
        data[:,0] = ( PackedTime - TimeOff ) / TimeScl
    else:
        data[:,0] = TimeOut1 + TimeIncr*(np.arange(0,NT))
        
        
    # use a dictionary for convenient access of channel columns (eg. data[t][ch['PtfmPitch'] )
    ch = dict(zip(channels, range(len(channels))))
    

    return data, ch, channels, units

    
    
    

def read_marin_file(dirName,fileName, zRefShift=0):

    # load data from marin ascii file of test results
    # zRefShift (added Oct28,2014) should be >0 if data file is about a CG below water line
    
    print('attempting to load '+dirName+fileName)
    
    f = open(dirName+fileName, 'r')
    
    # dict of MARIN channel codes
    MARIN_channel_code = ['181',  '1',   
        '36',          '37',          '38',
        '158',         '159',         '160',
        '60',          '61',          '62',
        '21', '30', '34'] 
    MARIN_channel_name = ['VXWIND REF', 'WAVE CL', 
        'X COG SEMI',  'Y COG SEMI',  'Z COG SEMI', 
        'ROLL SEMI',   'PITCH SEMI',  'YAW SEMI',
        'FMOOR 1',     'FMOOR 2',     'FMOOR 3',
        'AX TOP', 'FX TURBINE', 'MY TURBINE']
    FAST_channel_name  = ['WindVxi',    'WaveElev', 
        'PtfmSurge',   'PtfmSway',    'PtfmHeave',  
        'PtfmRoll',    'PtfmPitch',   'PtfmYaw',
        'Fair1Ten',    'Fair2Ten',    'Fair3Ten',
        'YawBrAxp','YawBrFxp','YawBrMyp']
    
    M_to_F = dict(zip(MARIN_channel_code, FAST_channel_name))  # map codes to fast names
    
    
    channels = []
    Mchannels = []
    units = []
    data = []
    i=0
    
    for line in f:          # loop through lines in file
    
        if ((i==0) or (i==2)):  # all header lines but channel code
            pass
            
        elif (i == 1):                      # channel code line - use this to make map
            channels.append("Time")         # first column ISN'T labelled and is always time (note: it has a unit entry though)
            for entry in line.split():      # loop over the elemets, split by whitespace
                #print 'entry is'+entry
                if entry in MARIN_channel_code:
                    channels.append(M_to_F[entry])      # append FAST ch name to the last element of the list
                else:
                    channels.append(entry)
                    
                Mchannels.append(entry) # if we don't know what the code means, just append it as-is
                
        elif (i == 3):                      # units
            for entry in line.split():      # loop over the elemets, split by whitespace
                units.append(entry)         # append to the last element of the list
        
        elif len(line.split()) > 0:
            data.append(np.zeros(len(channels)))                # add a new sublist to the data matrix
            
            entries = line.split()
            for jj in range(len(entries)): # loop over the elemets, split by whitespace
                if (channels[jj].count('PtfmSurge') > 0 or channels[jj].count('PtfmSway') > 0 or
                    channels[jj].count('PtfmRoll') > 0 or channels[jj].count('PtfmPitch') > 0): # channels to flip
                    data[-1][jj] = -1.0*float(entries[jj]) 
                else:
                    data[-1][jj] = float(entries[jj]) 
    
        else:
            break
    
        i+=1
    
    f.close()  # close data file

    #print str(M_to_F)
    #print str(channels)
    
    # use a dictionary for convenient access of channel columns (eg. data[t][ch['PtfmPitch'] )
    ch = dict(zip(channels, range(len(channels))))
    Mch = dict(zip(Mchannels, range(1,len(channels))))
        
    data2 = np.array(data)
    
    data3 = data2.astype(float)
    
    # adjust displacements to account for changed reference point
    if zRefShift != 0:
        if ('PtfmSurge' in ch) and ('PtfmPitch' in ch):
            data3[:,ch['PtfmSurge']] = data3[:,ch['PtfmSurge']] + zRefShift*np.sin(np.pi/180.*data3[:,ch['PtfmPitch']])
            #print "adjusting MARIN surge data for zRef offset"
        if ('PtfmSway' in ch) and ('PtfmRoll' in ch):
            data3[:,ch['PtfmSway']] = data3[:,ch['PtfmSway']] - zRefShift*np.sin(np.pi/180.*data3[:,ch['PtfmRoll']])
            #print "adjusting MARIN roll data for zRef offset"
        if ('PtfmHeave' in ch) and ('PtfmRoll' in ch) and ('PtfmPitch' in ch):
            data3[:,ch['PtfmHeave']] = data3[:,ch['PtfmHeave']]*np.cos(np.pi/180.*data3[:,ch['PtfmRoll']])*np.cos(np.pi/180.*data3[:,ch['PtfmPitch']])
            #print "adjusting MARIN heave data for zRef offset"
    
    return data3, ch, Mch, channels, units

    
def read_marin_file2013(dirName,fileName, zRefShift=0):

    # this version is for the updated channels used in the 2013 UMaine/MARIN data
    # load data from marin ascii file of test results
    # zRefShift (added Oct28,2014) should be >0 if data file is about a CG below water line
    
    print('attempting to load '+dirName+fileName)
    
    f = open(dirName+fileName, 'r')
    
    # dict of MARIN channel codes
    MARIN_channel_code = ["1","158", "159","160", "26","805","806","807", 
                          "35   ","38", "304","305", "306","501", "502", "503", 
                          "504", "505", "506", "507", "508", "509", "510", "511", "512"]
    MARIN_channel_name = ["WAVE 180","ROLL SEMI","PITCH SEMI","YAW SEMI","AX TOP","FLINE BOWT",
                          "FLINE PSAT","FLINE SBAT","RPM","PITCHBLADE","X MWL","Y MWL",
                          "Z MWL","FX TOP","FY TOP","FZ TOP","MX TOP","MY TOP",
                          "MZ TOP","FX BOT","FY BOT","FZ BOT","MX BOT","MY BOT","MZ BOT"]
    FAST_channel_name  = ["WaveElev","PtfmRoll","PtfmPitch","PtfmYaw",
                           "YawBrAxp","Fair1Ten","Fair2Ten","Fair3Ten",
                           "RotorSpd","BldPitch","PtfmSurge","PtfmSway",
                           "PtfmYaw","YawBrFxp","YawBrFyp","YawBrFzp",
                           "YawBrMxp","YawBrMyp","YawBrmzp","TwrBsFxp",
                           "TwrBsFyp","TwrBsFzp","TwrBsMxp","TwrBsMyp",
                           "TwrBsMzp"]
    
    M_to_F = dict(zip(MARIN_channel_code, FAST_channel_name))  # map codes to fast names
    
    
    channels = []
    Mchannels = []
    units = []
    data = []
    i=0
    
    for line in f:          # loop through lines in file
    
        if ((i==0) or (i==2)):  # all header lines but channel code
            pass
            
        elif (i == 1):                      # channel code line - use this to make map
            channels.append("Time")         # first column ISN'T labelled and is always time (note: it has a unit entry though)
            for entry in line.split():      # loop over the elemets, split by whitespace
                #print 'entry is'+entry
                if entry in MARIN_channel_code:
                    channels.append(M_to_F[entry])      # append FAST ch name to the last element of the list
                else:
                    channels.append(entry)
                    
                Mchannels.append(entry) # if we don't know what the code means, just append it as-is
                
        elif (i == 3):                      # units
            for entry in line.split():      # loop over the elemets, split by whitespace
                units.append(entry)         # append to the last element of the list
        
        elif len(line.split()) > 0:
            data.append(np.zeros(len(channels)))                # add a new sublist to the data matrix
            
            entries = line.split()
            for jj in range(len(entries)): # loop over the elemets, split by whitespace
                if (channels[jj].count('PtfmSurge') > 0 or channels[jj].count('PtfmSway') > 0 or
                    channels[jj].count('PtfmRoll') > 0 or channels[jj].count('PtfmPitch') > 0): # channels to flip
                    data[-1][jj] = -1.0*float(entries[jj]) 
                else:
                    data[-1][jj] = float(entries[jj]) 
    
        else:
            break
    
        i+=1
    
    f.close()  # close data file

    #print str(M_to_F)
    #print str(channels)
    
    # use a dictionary for convenient access of channel columns (eg. data[t][ch['PtfmPitch'] )
    ch = dict(zip(channels, range(len(channels))))
    Mch = dict(zip(Mchannels, range(1,len(channels))))
        
    data2 = np.array(data)
    
    data3 = data2.astype(float)
    
    # adjust displacements to account for changed reference point
    if zRefShift != 0:
        if ('PtfmSurge' in ch) and ('PtfmPitch' in ch):
            data3[:,ch['PtfmSurge']] = data3[:,ch['PtfmSurge']] + zRefShift*np.sin(np.pi/180.*data3[:,ch['PtfmPitch']])
            #print "adjusting MARIN surge data for zRef offset"
        if ('PtfmSway' in ch) and ('PtfmRoll' in ch):
            data3[:,ch['PtfmSway']] = data3[:,ch['PtfmSway']] - zRefShift*np.sin(np.pi/180.*data3[:,ch['PtfmRoll']])
            #print "adjusting MARIN roll data for zRef offset"
        if ('PtfmHeave' in ch) and ('PtfmRoll' in ch) and ('PtfmPitch' in ch):
            data3[:,ch['PtfmHeave']] = data3[:,ch['PtfmHeave']]*np.cos(np.pi/180.*data3[:,ch['PtfmRoll']])*np.cos(np.pi/180.*data3[:,ch['PtfmPitch']])
            #print "adjusting MARIN heave data for zRef offset"
    
    return data3, ch, Mch, channels, units

    
    
# function to read a FAST-style input file
def read_input_file(dirName,fileName):
    
    # load parameters from FAST input file
    
    print('attempting to read '+dirName+fileName)
    
    f = open(dirName+fileName, 'r')
    
    values = []
    names = []
    comments = []
    
    for line in f:          # loop through lines in file
    
        if line.count('- OUTPUT -') > 0:
            break
            
        elif line.count('---') > 0:
            pass
            
        elif len(line.split()) >= 2:    # need at least two items per line: value and name
            entry = line.split(None, 3)     # split line by whitespace 3 times
            values.append(entry[0])         # append ...
            names.append(entry[1])      # append ...
            
            if len(entry) >= 4:
                comments.append(entry[3])       # append ...
                
            #else:
                #comments.append("")
    
        #else:
        #   break
            
    
    f.close()  # close data file
    
    #f = None
    
    #print channels
    
    #print data[2][4]
    
    
    # use a dictionary for convenient access of parameters (eg. values[p['HubRad'] )
    p = dict(zip(names, range(len(names))))
    
    return values, p
    


# function to read a MAP-style mooring line input file
def read_lines_file(dirName,fileName):
    
    # load parameters from FAST input file
    
    print('attempting to read '+dirName+fileName)
    
    f = open(dirName+fileName, 'r')
    
    LineNumb = [] 
    LineNodes = []  
    RodNumb = [] 
    RodNodes = []
    
    for line in f:          # loop through lines in file
    
        # get properties of each line
        if line.count('- LINE PROPERTIES -') > 0:
            line = next(f) # skip this header line, plus channel names and units lines
            line = next(f)
            line = next(f)
            while line.count('---') == 0:
                entry = line.split()
                LineNumb.append(entry[0])
                LineNodes.append(entry[3]);
                line = next(f)
            
        # get properties of each rod
        if line.count('- ROD PROPERTIES -') > 0:
            line = next(f) # skip this header line, plus channel names and units lines
            line = next(f)
            line = next(f)
            while line.count('---') == 0:
                entry = line.split()
                RodNumb.append(entry[0])
                RodNodes.append(entry[9])
                line = next(f)
    
    f.close()  # close data file
    
    
    
    
    print('read '+str(len(LineNumb))+' line properties.')
    
    return LineNumb, LineNodes, RodNumb, RodNodes

def read_mooring_file(dirName,fileName):

    # load data from time series for single mooring line
    
    print('attempting to load '+dirName+fileName)
    
    f = open(dirName+fileName, 'r')
    
    channels = []
    units = []
    data = []
    i=0
    
    for line in f:          # loop through lines in file
    
        if (i == 0):
            for entry in line.split():      # loop over the elemets, split by whitespace
                channels.append(entry)      # append to the last element of the list
                
        elif (i == 1):
            for entry in line.split():      # loop over the elemets, split by whitespace
                units.append(entry)         # append to the last element of the list
        
        elif len(line.split()) > 0:
            data.append([])  # add a new sublist to the data matrix
            
            r = re.compile(r"(?<=\d)\-(?=\d)")  # catch any instances where a large negative exponent has been written with the "E"
            line2 = r.sub("E-",line)            # and add in the E
            
            
            for entry in line2.split():      # loop over the elemets, split by whitespace
                data[-1].append(entry)      # append to the last element of the list
    
        else:
            break
    
        i+=1
    
    f.close()  # close data file
    
    
    # use a dictionary for convenient access of channel columns (eg. data[t][ch['PtfmPitch'] )
    ch = dict(zip(channels, range(len(channels))))
    
    data2 = np.array(data)
    
    data3 = data2.astype(float)
    
    return data3, ch, channels, units
    

# peak detection function
def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    Returns two arrays
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    % [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    % maxima and minima ("peaks") in the vector V.
    % MAXTAB and MINTAB consists of two columns. Column 1
    % contains indices in V, and column 2 the found values.
    %
    % With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    % in MAXTAB and MINTAB are replaced with the corresponding
    % X-values.
    %
    % A point is considered a maximum peak if it has the maximal
    % value, and was preceded (to the left) by a value lower by
    % DELTA.
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    """
    maxY = []
    maxX = []
    minY = []
    minX = []
    
    if x is None:
        x = np.arange(len(v))
        
    v = np.asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
        
    if not np.isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
        
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
        
    mn, mx = np.Inf, -np.Inf
    mnpos, mxpos = np.NaN, np.NaN
    
    lookformax = True
    
    for i in np.arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        if lookformax:
            if this < mx-delta:
                maxX.append(mxpos)
                maxY.append(mx)
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                minX.append(mnpos)
                minY.append(mn)
                mx = this
                mxpos = x[i]
                lookformax = True
     
    return np.array(maxX), np.array(minX), np.array(maxY), np.array(minY)
    

# my damping ratio calculation

def dampingRatio(Xdata, delta, Tdata):
    
    # call peak detection function
    maxX, minX, maxY, minY = peakdet(Xdata, delta, Tdata)   # find peaks greater than 0.1 m

    zetas = []
    iamps = []
    
    for ii in range(len(maxY)-1):
        x0 = maxY[ii]
        x1 = maxY[ii+1]
        zeta = 1./np.sqrt(1.0 + (2.*np.pi/np.log(x0/x1))**2.0 ) # damping ration from logarithmic decrement between adjacent peaks
        zetas.append(zeta)
        iamps.append(abs(x0))
    for ii in range(len(minY)-1):
        x0 = minY[ii]
        x1 = minY[ii+1]
        zeta = 1./np.sqrt(1.0 + (2.*np.pi/np.log(x0/x1))**2.0 ) # damping ration from logarithmic decrement between adjacent peaks
        zetas.append(zeta)
        iamps.append(abs(x0))
        
    zetas = np.array(zetas)
    iamps = np.array(iamps)
    sorted_inds = np.argsort(iamps)
    zetas = zetas[sorted_inds]
    iamps = iamps[sorted_inds]

    return iamps, zetas
    
    

def make_psd(xdata, ydata, ylabel='', figname='PSD_plot'):
    # this function creates a PSD of a selected time series
    
    #print 'plotting...'
    
    
    dt = xdata[2]-xdata[1]  
    
    #print 'length of time series is '+str(ydata.size)
    
    nfft = int(2*round(0.5*(ydata.size/2)))
    overlap = int(nfft/2)
    
    
    psd, freqs = mlab.psd(ydata, NFFT=nfft, Fs=1.0/dt, detrend=mlab.detrend_linear, window=mlab.window_hanning, noverlap=overlap)   # psd(x, NFFT=256, Fs=2, Fc=0, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=0, pad_to=None, sides='default', scale_by_freq=None, **kwargs)
    
    
    return psd, freqs
    
    
def plot_psd(xdata, ydata, ylabel='', figname='PSD_plot', xlim=[0,2.5]):
    
    fig = plt.figure()
    plt.plot(xdata, ydata)
    plt.xlim(xlim)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel(ylabel)
    
    plt.show()
    
    '''
    ax = plt.axes()
    #ax3 = plt.subplot(313)
    #ax.clear()
    plt.plot(xdata,ydata)
    ax.set_xlim([0, 2.5])
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel(ylabel)
    
    
    plt.savefig(figname+'.png', bbox_inches=0)
    print("Figure saved.")
    
    #plt.show()
    
    ax.clear()
    '''
    #print "Hs is "+str(4*np.sqrt( np.sum(psd) * (freqs[2]-freqs[1]) ) ) # check
    
    
    #fftP = dt/ydata.size * np.abs(np.fft.rfft(ydata))**2
    #fftX = np.linspace(0, 1./2./dt, fftP.size)  #np.fft.fftfreq(ydata.size, dt)
    #
    #print 'len P '+str(fftP.size)+' and len X '+str(fftX.size)
    #
    #plt.figure()
    #ax2 = plt.axes()
    #plt.plot(fftX, fftP)
    #smoothed = np.convolve(fftP, [0.333,0.333,0.333], mode='same')
    #plt.figure()
    #ax3 = plt.axes()
    #plt.plot(fftX, smoothed)
    #ax2.set_xlim([0, 0.4])
    ##ax3.set_xlim([0, 0.4])
    #print str(np.sum(fftP))
    #print str(np.sum(smoothed))
    
    
    #plt.show()

#def seperate_motions(data, ch, hHub)
    # this function examines the contribution of platform 
    # surge and pitch versus tower bending to nacelle motions
    
#   RMS_nac = np.rms(data[:,ch['...
    
    
# create rotation matrix
def RotMat( x2, x1, x3 ):

    # note above swapping of x1 and x2 to deal with weird coordinate system from FAST convention
    # ( x2 is roll, x1 is pitch, x3 is yaw )

    s1 = np.sin(x1) 
    c1 = np.cos(x1)
    s2 = np.sin(x2) 
    c2 = np.cos(x2)
    s3 = np.sin(x3) 
    c3 = np.cos(x3)
    
    #rmat = transpose([ 1 0 0; 0 c1 -s1; 0 s1 c1] * [c2 0 s2; 0 1 0; -s2 0 c2] * [c3 s3 0; -s3 c3 0; 0 0 1]);
    
    TransMat = np.matrix(np.zeros([3,3]))
    
    TransMat[0,0] =  c1*c3 + s1*s2*s3;   # just guessed the order (may be transposed!) <<<<<<<<<
    TransMat[0,1] =  c3*s1*s2-c1*s3;
    TransMat[0,2] =  c2*s1;
    TransMat[1,0] =  c2*s3;
    TransMat[1,1] =  c2*c3;
    TransMat[1,2] =  -s2;
    TransMat[2,0] =  c1*s2*s3 - c3*s1;
    TransMat[2,1] =  s1*s3 + c1*c3*s2;
    TransMat[2,2] =  c1*c2;
    
    return TransMat

# create a 6 by 6 mass matrix
def MassMat(ms, xs, ys, zs, Ix=0, Iy=0, Iz=0):

    m = np.sum(ms)
    n = len(ms)
    xs2 = xs*xs
    ys2 = ys*ys
    zs2 = zs*zs

    M = np.zeros([6,6])
    # mass
    M[0,0] = M[1,1] = M[2,2] = m
    # cross coupling
    M[0,4] = M[4,0] = -np.sum(ms*zs) # surge-pitch coupling!
    M[0,5] = M[5,0] =  np.sum(ms*ys) # surge-yaw
    M[1,3] = M[3,1] =  np.sum(ms*zs) # sway-roll
    M[1,5] = M[5,1] = -np.sum(ms*xs) # sway-yaw
    M[2,3] = M[3,2] = -np.sum(ms*ys) # heave-roll
    M[2,4] = M[4,2] =  np.sum(ms*xs) # heave-pitch
    # inertias (from https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor)
    M[3,3] = np.sum(ms*(      ys2 + zs2)) + Ix
    M[4,4] = np.sum(ms*(xs2       + zs2)) + Iy
    M[5,5] = np.sum(ms*(xs2 + ys2      )) + Iz
    M[3,4] = M[4,3] = -np.sum(ms*xs*ys   )
    M[3,5] = M[5,3] = -np.sum(ms*xs   *zs)
    M[4,5] = M[5,4] = -np.sum(ms*   ys*zs)

    return M

def PlotChannels(dat, ch, clist):
    plt.figure()
    for chan in clist:
        plt.plot(dat[:,ch["Time"]], dat[:,ch[chan]], label=chan)
        
    plt.legend()
    plt.show
    
    return
    
    
def PlotFilesChannels(directory, filenames, channels, skiplines=-1, hasunits=1, chanlim=999, makepsd=1, figName='', figsize=None, range=(0,-1)):
    # simply compares files


    fig, ax = plt.subplots(len(channels),1, figsize=figsize, sharex=True)
    if makepsd>0:
        fig2, ax2 = plt.subplots(len(channels),1, sharex=True)

    for filename in filenames:

        data, ch, chans, units = read_output_file(directory, filename, skiplines=skiplines, hasunits=hasunits,chanlim=chanlim, dictionary=False)
        
        for ic, chan in enumerate(channels):
            ax[ic].plot(data[range[0]:range[1],0], data[range[0]:range[1],ch[chan]], label=filename)
            
            if makepsd>0:
                SurgeQSpsd, fsQS = make_psd(data[:,0], data[:,ch[chan]])
                ax2[ic].loglog(fsQS, SurgeQSpsd, label=filename)
            
    #ax[0].legend()  
    if makepsd>0:
        ax2[ic].legend()
    
    for ic, chan in enumerate(channels):
        if 'Ptfm' in chan:
            chan = chan.replace('Ptfm','')
        if chan=='Roll' or chan=='Pitch' or chan=='Yaw':
            chan = chan+' (deg)'
        if chan=='Surge' or chan=='Sway' or chan=='Heave':
            chan = chan+' (m)'
        if 'Wave' in chan:
            chan = 'Wave Elevation (m)'
        if 'Wind' in chan:
            chan = f'Wind Velocity_{chan[-1]} (m/s)'
        if chan[0]=='L' and isinstance(int(chan[1]), int):
            chan = 'Line'+chan[1]+' Tension (N)'
            
        ax[ic].set_ylabel(chan)
        
        if makepsd>0:
            ax2[ic].set_ylabel(chan)
    
    
    ax[-1].set_xlabel('Time (s)')
    plt.tight_layout()
    plt.show()
    if len(figName) > 0:
        if os.path.isdir(directory+'figures'):
            pass
        else:
            os.mkdir(directory+'figures')
        fig.savefig(directory+'figures/'+figName)


def wavenumber(w,h,e=0.001):

    #takes wave frequency in rad/s and returns wave number%%%

    g = 9.81
    omega = w                       #angular frequency of waves
    k = 0                           #initialize  k outside of loop
                                    #error tolerance
    k1 = omega*omega/g                  #deep water approx of wave number
    k2 = omega*omega/(np.tanh(k1*h)*g)      #general formula for wave number
    while np.abs(k2 - k1)/k1 > e:           #repeate until converges
        k1 = k2
        k2 = omega*omega/(np.tanh(k1*h)*g)
    
    k = k2
    
    return k
    
    
    

# get lists with channel indices of line node positions for use in plotting
def getLineNodeCh(ch,n):
    
    Chx = []
    Chy = []
    Chz = []
    
    for i in range(n):
        Chx.append(ch['Node'+str(i)+'px'])
        Chy.append(ch['Node'+str(i)+'py'])
        Chz.append(ch['Node'+str(i)+'pz'])
        
    return Chx, Chy, Chz

'''
def lineVectors(data,ch,n,ts):

    xp = np.zeros([n])
    yp = np.zeros([n])
    zp = np.zeros([n])
    
    for i in range(n):
        xp[i] = data[ts][ch['Node'+str(i)+'px']]
        yp[i] = data[ts][ch['Node'+str(i)+'py']]
        zp[i] = data[ts][ch['Node'+str(i)+'pz']]    
    
    return xp,yp,zp
'''


    
def animateLines(moordirname, moorfilename):

        
    # class to hold data of each mooring line or rod
    class Mooring():
        
        #def __init__(self, vals, p, data, ch, channels, units):
        def __init__(self, number, length, diameter, nNodes, data, ch, channels, units):
            
            self.number = number
            self.length = length   # to be used in future for visualizing strain of elements!
            self.d = diameter
            self.Nnodes = int(nNodes) + 1 #int(vals[p["NumNodes"]])
            
            # get time info
            if ("Time" in ch):
                self.Tdata = data[:,ch["Time"]]
                self.dt = self.Tdata[1]-self.Tdata[0]
            else:
                print("ERROR: could not find Time channel for mooring line "+str(number))
        
            #self.data = data
            #self.ch = ch
            #self.channels = channels
            #self.units = units
            
            nT = len(self.Tdata)  # number of time steps
            
            self.xp = np.zeros([nT,self.Nnodes])
            self.yp = np.zeros([nT,self.Nnodes])
            self.zp = np.zeros([nT,self.Nnodes])
            
            
            for i in range(self.Nnodes):
                self.xp[:,i] = data[:, ch['Node'+str(i)+'px']]
                self.yp[:,i] = data[:, ch['Node'+str(i)+'py']]
                self.zp[:,i] = data[:, ch['Node'+str(i)+'pz']]


            self.xpi= self.xp[0,:]
            self.ypi= self.yp[0,:]
            self.zpi= self.zp[0,:]
            
    
            
            
    # load mooring info file
    LineNumb, nNodes, RodNumb, RodNodes = read_lines_file(moordirname, moorfilename)
    print('planning to read '+str(len(LineNumb))+" Lines with segments: "+str(nNodes))
    print('planning to read '+str(len(RodNumb))+" Rods with segments: "+str(RodNodes))
    
    #data, ch, channels, units = read_output_file(moordirname, "lines.out")     # load general mooring output file
    #self.extras.append( Extras(data, ch, channels, units) )  # append to list of "extra" data structures


    moorings = []   # list of Mooring objects containing line and rod data


    for i in range(int(len(LineNumb))):     #range(NumLines):
    
        # load mooring line timeseries
        data3, ch3, channels3, units3 = read_mooring_file(moordirname, "Line"+str(LineNumb[i])+".out") # starts on 1 rather than 0
        
        # make mooring line object
        moorings.append( Mooring(LineNumb[i], 0.0, 0.2, nNodes[i], data3, ch3, channels3, units3) )
                                                                                                                                  
    for i in range(int(len(RodNumb))):   
        data3, ch3, channels3, units3 = read_mooring_file(moordirname, "Rod"+str(RodNumb[i])+".out") # starts on 1 rather than 0
        moorings.append( Mooring(RodNumb[i], 0,  0.5, RodNodes[i], data3, ch3, channels3, units3) )
    
    
    ###############################################################################
    
        # Definition of the function to get the x,y and z coordinates of each node at each time interval "dt"
    def update_Coords(tStep, moorings, lines):
        #cont = 0; cont1 =0
        #for LinesFilePack in range(0,MoorLines.shape[0]):
        #for imooring in moorings:
        #
        #   imooring.xpi= imooring.xp[tStep,:]
        #   imooring.ypi= imooring.yp[tStep,:]
        #   imooring.zpi= imooring.zp[tStep,:]
        #
        #breakpoint()
        
        print ("update_Coords called with t step "+str(tStep))
        for line, imooring in zip(lines, moorings) :
            # NOTE: there is no .set_data() for 3 dim data...
            line[0].set_data(imooring.xp[tStep,:], imooring.yp[tStep,:])  # set x and y coordinates
            line[0].set_3d_properties(imooring.zp[tStep,:])               # set z coordinates
            
        return lines
    
    ###############################################################################
    
    # Prepare the figure      
    fig = plt.figure(figsize=(20/2.54,12/2.54))
    #ax = fig.gca(projection='3d')
    ax = Axes3D(fig)
    
    lines = [ax.plot(imooring.xpi, imooring.ypi, imooring.zpi, 'bo', markersize =0.7, markerfacecolor='green') for imooring in moorings]
    
    ax.set_xlim((-15,15)); ax.set_ylim((-15,15)); 
    ax.set_zlim((-40,0))
    ax.set_xlabel('x');    ax.set_ylabel('y');      ax.set_zlabel('z');
    #fig.tight_layout()   
    #text_time = ax.text(1, 1, 1, 'Time = 0.00', transform=ax.transAxes)


      
    #pdb.set_trace()
      
    # Animation: update the figure wit the updated coordinates from update_Coords function 
    line_ani = animation.FuncAnimation(fig, update_Coords, np.arange(1,len(moorings[0].Tdata)-1), fargs=(moorings, lines),
                                       interval=10, blit=False)#, repeat_delay=300)
    #plt.show()
    
    


    
    