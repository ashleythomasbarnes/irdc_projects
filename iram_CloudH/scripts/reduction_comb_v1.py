## Source gildas installation - outside of python
# source /vol/software/software/astro/gildas/initgildas-nov18a.sh
# source /vol/software/software/astro/gildas/initgildas-sep23a.sh
## Add gildas installation to python path
# export PYTHONPATH=$PYTHONPATH:/vol/software/software/astro/gildas/gildas-exe-sep23a/x86_64-ubuntu22.04-gfortran/python
# export PYTHONPATH=$PYTHONPATH:/vol/software/software/astro/gildas/gildas-exe-sep23a/x86_64-ubuntu22.04-gfortran/python
##

## Import  pygildas  modeules
import pygildas
import pyclass
##

## General imports
from glob import glob
import numpy as np
import os
import astropy.utils.console as console
##

####################################################
###### User inputs
# Dirs
inputdir = './../data_splitbaddata'
outputdir = './../data_processed'
inputfiles = glob('%s/*.30m' %inputdir)

# Source
source = 'CLOUDH'
source_out = 'CloudH'
ra0 = 18 + (57 + 7.8/60.) / 60. #30.830
dec0 = 2 + (10. + 39.1/60.) / 60. #3.53806
vlsr = 45.0

obs_param = {
    '13CS':             {'freq': 92494.30800, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'CCS_93':           {'freq': 93870.09800, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'CCD':              {'freq': 72107.72050, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'c-C3H':            {'freq': 91494.34900, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'HNCO_87_10K':      {'freq': 87925.23700, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'HNCO_87_53K':      {'freq': 87597.33000, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'DCOp':             {'freq': 72039.31220, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'HCOp':             {'freq': 89188.52470, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'H13COp':           {'freq': 86754.22840, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'DCN':              {'freq': 72414.90500, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'HCN':              {'freq': 88631.84750, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'H13CN':            {'freq': 86339.92140, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'HC15N':            {'freq': 86054.96100, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'H15NC':            {'freq': 88865.69200, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'HN13C':            {'freq': 87090.85000, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'DNC':              {'freq': 76305.72700, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'HNC':              {'freq': 90663.56800, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'H2CO':             {'freq': 72837.94800, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'HC3N_72':          {'freq': 72783.82200, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'HC3N_91':          {'freq': 90979.02300, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'N2Dp':             {'freq': 77109.24330, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'N2Hp':             {'freq': 93173.39770, 'vel_ext': '-30 130', 'vel_win':'30 58'},
    'SO2':              {'freq': 72758.24340, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'SiO':              {'freq': 86846.96000, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'SO':               {'freq': 86093.95000, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'CCH_873_a':        {'freq': 87284.10500, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'CCH_873_b':        {'freq': 87316.92500, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'CCH_873_c':        {'freq': 87328.62400, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'CH3CHO_935_a':     {'freq': 93580.90910, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'CH3CHO_935_b':     {'freq': 93595.23490, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'CH3OH_765':        {'freq': 76509.68400, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'CH3OH_924':        {'freq': 92409.49000, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    'H41a':             {'freq': 92034.43415, 'vel_ext': '-30 130', 'vel_win':'30 60'},
    }

lines = obs_param.keys()
#lines = ['hcop_10', 'dcn_10', 'n2hp_10', 'hcn_10', 'n2hp_10', ]
# lines = ['N2Dp', 'N2Hp']

# Output velocity axis
# https://publicwiki.iram.es/Iram30mEfficiencies
beff = 0.8
nchans = 2000  #number of channels
chan0 = 0      #reference channel
chan0_v = 0    #velocity at reference channel
delta_v = 0.05 #delta velocity
####################################################


####################################################
####################################################
### Initial
# get the gildas command line interpreter to issue commands quickly
sic = pyclass.comm
# get the pyclass.gdict object for easy Gildas variable access
g = pyclass.gdict
# sic('SIC MESSAGE GLOBAL ON') # or OFF for less noise
###
for line in lines:

    print('[INFO] Reducing line: %s' %line)

    #Get frequency
    freq = obs_param[line]['freq'] 
    vel_win = obs_param[line]['vel_win']
    vel_ext = obs_param[line]['vel_ext']

    ###Define output
    outputfile = '%s/%s_%s' %(outputdir, source_out, line)
    outputfile_fits = outputfile.replace('data_processed', 'cubes_processed')

    os.system('rm %s.30m' %outputfile)
    os.system('rm %s.lmv' %outputfile)
    os.system('rm %s.tab' %outputfile)
    os.system('rm %s.wei' %outputfile)
    os.system('rm %s.fits' %outputfile_fits)

    sic('file out %s s' %outputfile) #Do not use m here - m = multiple version in output, s = single, so s means you can not have two ids that are the same - e.g. if you merge files, a new ID is create instead of overwriting the old one.
    print('[INFO] Removing old output file')
    print('[INFO] Making new output file: %s' %outputfile)
    ####

    ###Loop through files - one file per date
    for inputfile in inputfiles:
        
        # if inputfile != inputfiles[2]:
        #    continue

        #Load file and det defaults
        sic('file in %s' %inputfile)
        sic('set source %s' %source)
        sic('set telescope %s' %'*')
        sic('set line %s' %'*')

        #only take data in range of map - x1 x2 y1 y2
        sic('set range %s %s %s %s' %(-200, 200, '*', '*'))

        #Open and check file
        sic('find /frequency %s' %freq) #only obs with refrenecy
        if g.found == 0:
            print('[INFO] No data found!')
            continue

        #Loop through spectral, baseline, and output
        inds = g.idx.ind.tolist()
        sic('sic message class s-i') #! Speed-up long loops by avoiding too many screen informational messages

        for ind in console.ProgressBar(inds):

            sic('set mo x t') #reset range to stop warning
            sic('get %s' %ind)
            sic('modify freq %s' %freq)
            sic('modify source %s' %source_out)
            sic('modify projection = {0} {1} ='.format(ra0, dec0))
            sic('modify telescope 30M-MRT')
            sic('modify linename %s' %line)

            # sic('modify beam_eff %s' %beff) #set manually
            sic('modify beam_eff /ruze') #auto beff determine
            sic('set unit v')
            sic('extract %s v' %vel_ext) #cut out spectra
            sic('resample %i %i %f %f v' %(nchans, chan0, chan0_v, delta_v))
            sic('set mo x 0 100') #only do baseline on this part
            sic('set win %s' %vel_win) #cut out spectra
            try:
                sic('base 1')
            except:
                print('[INFO] No basline fitted.')
                continue
            sic('write')
        sic('sic message class s+i') # Toggle back screen informational messages

    sic('file in %s.30m' %outputfile) # load 30m file
    sic('find /all') #regrid and output to fits file
    if g.found == 0: #some lines may not be covered
        print('[INFO] No data found!')
        continue
    sic('table %s new /resample %i %i %f %f velo /nocheck' %(outputfile, nchans, chan0, chan0_v, delta_v))
    sic('xy_map %s' %outputfile)
    sic('vector\\fits %s.fits from %s.lmv' %(outputfile_fits, outputfile))
    ###

