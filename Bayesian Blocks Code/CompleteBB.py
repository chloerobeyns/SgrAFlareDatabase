#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 21:22:24 2020

@author: Elisa

This is a modification of Hope Boyce's and Eli Bouffard's algorithm 
In turn, each of these was based on Peter William's work. 

In this version of the file, I contain all the necessary functions to 
get all the information we want from the bayseian block analysis 
for the x-ray flare database. 

Functions in this code are for plotting & creating/writing files.

Needs to be used in conjunction with xbblocks.py file. 

FUNCTIONS: 
    
1. Process ->> creates an 'output file' with bb blocks info. 
2. plot_bb ->> plots the  bayesian blocks on axes
3. plot_lc ->> plots the lightcurve on axes 
3. getInfo ->> creates an output file that contains all the
                information needed to populate the flare database.
 
"""

import matplotlib.pyplot as plt
import sys
import numpy as np
import numpy.ma as ma
import math



from astropy.io import fits
from astropy.table import Table

import xbblocks #this is a set of functions written by Eli Bouffard

def process(infile, outfile, p0=0.05):
    
    ''' 
    Generates the bayesian blocks and creates file with results 
    
    Parameters
    ----------
    
        infile   : name of input fits file (event file)
        
        p0       : default is 0.05
        
        outfile  : name of output text file into which results are 
                   written. 
                   
    Returns
    ----------
    Nothing, but creates a file with bayesian block results. 
    These results contain: 
        p0 value 
        timesys 
        start time of observation (MJD)
        end time of observation (MJD)
        time array length 
        
        Block info in order: 
            left edge of block mjd
            right edge of block  mjd
            counts in the block
            width of the block mjd
            count rate (avg) in block ct/day 
            count rate error ct/day 
    
    '''
    
    
    f = fits.open(infile)
    o = open(outfile, mode='w')


    timesys = f[0].header['timesys']
    mjdref = f[0].header['mjdref']
    print(mjdref)
    timeunit = f[0].header['timeunit']
    timezero = f[0].header['timezero']
    tstart = f[0].header['tstart']
    tstop = f[0].header['tstop']

    if timeunit == 's':
        tscale = 1. / 86400 #to convert s to days 
        #24*60*60 = 86400s/day so that timestamps can be added to mjdref
    else:
        die ('can\'t handle time unit "%s" in input "%s"', timeunit)

    eventhdu = None

    for hdu in f[1:]:
        if hdu.name == 'EVENTS':
            if eventhdu is None:
                eventhdu = hdu
            else:
                die('input "%s" has multiple EVENTS sections; don\'t know '
                     'which to use',)

    if eventhdu is None:
        die('input "%s" has no EVENTS sections')

    ccdid = eventhdu.data.ccd_id.min ()
    if eventhdu.data.ccd_id.max () != ccdid:
        die ('can\'t handle data from multiple CCDs in input "%s"')

    gtihdu = None

    for hdu in f[1:]:
        if hdu.name == 'GTI' and hdu.header['ccd_id'] == ccdid:
            if gtihdu is None:
                gtihdu = hdu
            else:
                die ('input "%s" has multiple matching GTI sections; don\'t know '
                     'which to use')

    if gtihdu is None:
        print (sys.stderr, 'warning: no GTI info for active CCD %d; trusting ' \
            'TSTART and TSTOP instead' % ccdid)
        tstarts = np.asarray ([tstart])
        tstops = np.asarray ([tstop])
    else:
        tstarts = gtihdu.data.START #defined in terms of good time interval
        tstops = gtihdu.data.STOP

    times = (eventhdu.data.time + timezero) * tscale + mjdref
    tstarts = (tstarts + timezero) * tscale + mjdref
    tstops = (tstops + timezero) * tscale + mjdref

    info = xbblocks.bsttbblock (times, tstarts, tstops, p0=p0,
                                nbootstrap=256)

    print('# p0 = %g' % p0, file=o)
    print('# timesys =', timesys, file=o)
    print('# tstarts =', ' '.join ('%.5f' % t for t in tstarts), file=o)
    print('#tstops  =', ' '.join ('%.5f' % t for t in tstops), file=o)
    print('# n = %d' % times.size, file=o)
    for i in range (info.nblocks):
        s = '%.5f %.5f %4d %g %g %g' % (info.ledges[i], info.redges[i],
                                        info.counts[i], info.widths[i],
                                        info.rates[i], info.bsrstds[i])
        
        print(s, file=o)
    o.close()



'''___________________________________________________________________
'''



def plot_bb(file):
	'''
	Plot Bayesian Blocks onto a plot

	Parameters
	----------   
    file     : string
		location and name of Bayesian Blocks results

	time_start : float
		MJD start of observation in days

	time_end : float
		MJD end of observation in days

	Returns
	-------
		nothing, plots the bayesian blocks. 

	'''

	### read Bayesian Blocks data in 
	ledges, redges, counts, widths, rates, bsrstds = np.transpose(np.loadtxt(file))

	bstart = ledges ### rename
	bend = redges ### rename
	rates = rates/86400.0 #convert from counts/day to counts/s

	#####
	### Recasting Bblocks results into arrays that are easier for plotting
	### x = time array
	### rates = the levels of blocks
	#####
	if hasattr(rates, "__len__"):
		### If there is more than one Bblock
		x = [bstart[0]]
		for j in range(0,len(rates)):
			x= np.concatenate([x,[bstart[j]]])
		x= np.concatenate([x,[max(bend)]])
		rates = np.concatenate([[0],rates,[0]])
	else:
		### If there's only one Bblock
		x = np.array([bstart,bstart,(bend)])

		rates = np.array([rates])
		rates = np.concatenate([[0],rates,[0]])

	##x_utc = Time(x, format='mjd', scale='utc') ## convert time array into UTC format

	###plotting
	plt.plot(x,rates,drawstyle='steps-post',color='#ff7b0f', lw=1.5,zorder=5) ### plot with MJD axis - timestart
	# plt.plot(x_utc.datetime,rates,drawstyle='steps-post',color='#ff7b0f', lw=1.5,zorder=5) ### plot with UTC axis
 


'''___________________________________________________________________
'''



def plot_lc(file):
    ''' Plot the lightcurve file 
    
    Parameters
    ----------
    
    file: fits file 
        File that contains the binned light curve data. 
        
    Since the bayesian blocks algorithm bins data by 300s, it is most convinent to use 300s 
    bins when plotting both on overlap. 
    '''
    
    f = fits.open(file)
    
    #obtain the mjdref from the primary header file : 
    #obtain the timezero correction from primary header: 
    primary = f[0].header
    mjdref = primary['MJDREF']
    timezero = primary['TIMEZERO']
    
    #obtain necessary data for plotting: 
    data = Table.read(f, hdu=1)
    time = mjdref + (timezero + data["TIME"] ) / 86400 #convert time data to mjd 
    counts = data["NET_RATE"]
    err = data["ERR_RATE"]
    
    
    
    #plot : 
    plt.errorbar(time, counts, yerr=err, marker=".", color="blue", mfc="black",mec="black", ecolor="navy")
    



'''___________________________________________________________________
'''



    
def getInfo(evtfile, lcfile, bbfile, outfile, countmin=8, amp_crit=3):
    ''' 
    Generate a document that contains all the data required for the flare
    table database. 
    
    Parameters
    ----------
    evtfile: fits file
        File on which the bayesian blocks was run. 
        
    lcfile: fits file 
        File with the lightcurve of the observation. Ideally, bin 300s
        
    bbfile: txt file 
        File that was output of the process function. Contains block info
    
    outfile: str 
        Name of the output file to generate. 
    
    countmin: int 
        Minimum number of counts in a block. If lower, it will be combined with another block
        Default value is 8 *(arbitrary)
    
    amp_crit: int 
        Sigma range above quiesence for a block to be considered a flare
        Default value is 3
        
    
    Returns
    -------
    Nothing, but generates output file with information needed for the 
    database. 
    
    '''
    
    #Open the event file to get the basic information about observation
    
    f = fits.open(evtfile)
    o = open(outfile, mode='w')
    primary = f[0].header
    data_info = f[1]
    evt = Table.read(f[1])
    gti = Table.read(f[2])
    
    g = fits.open(lcfile)
    dat = Table.read(g, hdu=1)
    count_rate_lc = dat["NET_RATE"]
    err = dat["ERR_RATE"]
    #convert times to mjd:
    
    mjdref_lc = primary['MJDREF']
    timezero_lc = primary['TIMEZERO']
    time_lc = mjdref_lc + (timezero_lc + dat["TIME"] ) / 86400
    
    mjdref_evt = primary['MJDREF']
    timezero_evt = primary['TIMEZERO']
    time_evt = mjdref_evt + (timezero_evt + evt["time"] ) / 86400
    
    ''' Get data from the primary header '''
    
    obsid = primary["OBS_ID"]
    date = primary["DATE-OBS"]
    instrument = primary["INSTRUME"]
    telescope = primary["TELESCOP"]
    exp = (gti["STOP"][0]-gti["START"][0])/1000 #the good time interval 
    
    #Add this data to the file: 
    print("-------------------------------------------------------------------------------------", "\n", file=o)
    print("BASIC INFORMATION", "\n", file=o)
    print("-------------------------------------------------------------------------------------", "\n", file=o)
    print("Obs_ID:", obsid, "\n", file=o)
    print("Obs Date:", date, "\n", file=o)
    print("Telescope:", telescope, "\n", file=o)
    print("Instrument:", instrument, "\n", file=o)
    print("Exposure (ks):", exp, "\n", file=o)
    
    '''Get data about the flares and Quiescent count rate'''
    
    #load bayesian blocks data 
    ledges, redges, counts, widths, rates, bsrstds = np.transpose(np.loadtxt(bbfile))
    #NOTE: the data that is given by 'process' function is in mjd 
    
    #Get the flare and block information:
    d,b,l = get_flare_bb_nobsnopcr(ledges, redges, counts, widths, rates, amp_crit, countmin)
    
    print("Quiescent Count Rate (10^-3 ct/s):", np.around((l[0][0]/86400.0)/(10**(-3)), 3), "+/-", np.around((l[0][1]/86400.0)/(10**(-3)), 3), "\n", file=o)
   # print("Quiescent Count Rate (10^-3 ct/s):", np.around((269.678/86400.0)/(10**(-3)), 3), "+/-", np.around((131.224/86400.0)/(10**(-3)), 3), "\n", file=o)
  
    #to find the number of flares: 
    if type(b) is type(None):
        print("No Flares in Observation \n", file=o)
      
    else:
        
        b = np.asarray(np.transpose(b))
        num_flare = b[0].size
    
        #for each flare, get the information 
        for i in range (0,num_flare):
            start = b[0][i] #start time of the flare
            end = b[1][i]  #end time of the flare 
        
            '''Picking out the mximum count rate from lightcurve:  '''
            #within this time region, need to find the max count rate:
            upper = np.where(time_lc <= b[1][i])
            lower = np.where(time_lc >= b[0][i]) #indices that say where 
            #within the light curve the lightcurves start/finish 
        
        
            #the maximum count rate can be found for times between indices min(lower) and max(upper)
            #then the maximum count rate is: 
            ct_max = np.max(count_rate_lc[np.min(lower):np.max(upper)])
            ct_max_err = err[np.where(count_rate_lc == ct_max)][0]
        
            '''Picking out Luminosity & Energy emitted by each flare: '''
            #According to Eli's paper there is a proportionality between 
            #counts in a block and the energy 
            #for non gratings: 0.013 ct / (10^34 erg) 
            #then the energy is given by: 
            ct_mean = b[4][i]/86400.0 #count rate in each block
            ct_mean_err = b[5][i]/86400.0 #err in mean ct rate
        
            ct_quies_sub = ct_mean - (l[0][0]/86400.0) #quiescence subtracted countrates
            dur = b[3][i]*86400.0 #duration of flare in s 
        
        
            #Determining values
            luminosity = (ct_quies_sub*(10**(34)/0.013)) # erg/s  
            energy = luminosity * dur  #erg
        
        
            #Error propagation
            lum_err = (np.sqrt((ct_mean_err)**2 + (l[0][1]/86400.0)**2))*(10**(34)/0.013)
            energy_err = energy * np.sqrt((lum_err/luminosity)**2)
            
        
            #Finding the Flux of each flare: 
            #flux is related to the distance between sgrA and the detector 
            #which is about radius = 2.45*10^22 cm 
            #the surface area over which we are dividing for the 
            #flux is 4*pi*radius^2 
            radius = 2.45*10**(22)
            err_rad = 4 *10**(19)
            flux = ((luminosity)/(4*np.pi*radius**2))# 10^-12 erg/s/cm2
            
            #Error propagation: 
            flux_err = flux*np.sqrt((lum_err/luminosity)**2 + (err_rad/radius)**2)
        
            #Finding the fluence of each flare: 
            #fluence is given in ct so I am assuming it is the total 
            #count number incident upon the detector 
            fluence = ct_mean*b[3][i]*86400.0
        
        
        
            print("-------------------------------------------------------------------------------------", "\n", file=o)
            print("FLARE NUMBER ", i+1,"\n", file=o)
            print("-------------------------------------------------------------------------------------", "\n", file=o)
            print("Start Time: ", start, "(MJD) \n", file=o)
            print("End Time: ", end, "(MJD) \n", file=o)
            print("Duration: ", dur, "(s) \n", file=o)
            print("Count Rate (mean): ",np.around(ct_mean,2) , "+/-", np.around(ct_mean_err, decimals=4), "(ct/s) \n", file=o)
            print("Count Rate (max): ",np.around(ct_max,2) ,"+/-", np.around(ct_max_err,2), "(ct/s) \n", file=o)
            print("Energy: ", np.around(energy/10**(37), 2), "+/-", np.around(energy_err/10**(37),2), "10^37 ergs \n", file =o)
            print("Luminosity: ", np.around(luminosity/10**(34),2), "+/-", np.around(lum_err/10**(34),2), "10^34 erg/s \n", file=o)
            print("Flux: ", np.around(flux/(10**(-12)),2) ,"+/-" ,np.around(flux_err/(10**(-12)),2) , "10^-12 erg/s/cm2 \n", file=o)
            print("Fluence: ", np.around(fluence,2), "ct \n", file=o)
        
    f.close()
    o.close()
    




'''___________________________________________________________________
'''



    
    
def get_flare_bb_nobsnopcr(ledges, redges, counts, widths, rates, amplitude_criteria = 3, minflu = 8):
    '''
    Outputs Flare parameters for a given set of parameters from a 
    Bayesian Blocks analysis without bootstrap. For subarray mode 
    (Not gratings!). 
    
    This function is taken from Eli Bouffard's code. 
    
    Helper function for getInfo. 
    
    Parameters
    ----------
    ledges: array of floats
        Time of the beginning of each block (s)
        
    redges: array of floats
        Time of the end of each block (s)
        
    counts: array of int
         Total counts in each block 
         
    widths: array of floats
         Total length of each block (s)
         
    rates: array of floats
        Mean count rate of each block (ct/s)

    amplitude_criteria: float:
        Sigma range above quiesence for a block to be considered a 
        flare
        Default value is 3
        
    minflu: int
        minimum number of counts in a block. If lower, it is combined 
        with a nearby block
        
    Returns
    ----------
    data:  array of floats
        Contains, in order, the time of the beginning of each block (s), 
        the time of the end of each block (s),
        the number of counts in each block, the total lenght of 
        each block (MJD), the mean count rate of each block (ct/MJD), the
        standard deviation in count rate
        in each block (ct/MJD) 
        
        (and the Poisson error in each block count
        rate (ct/MJD)) -> NOT
        
    block:  array of floats
        Same as data but for pile-up corrected values and for each 
        FLARE instead of each block (flares can be made of multiple 
        blocks)
        
    LoRate : array of floats
        Contains the mean count rate of the longest block and its 
        Poisson error
    '''    
    #Note how many blocks there are in the obsid
    num_blocks = np.size(redges)

    l = 0     #counts the number flaring blocks (the total, not the number of individual flares)
    
    counts = np.asarray(counts)
    widths = np.asarray(widths)
    rateserr = np.sqrt(counts)/widths
    block = None
    #If there are more than 1 block, then the ones significantly above the lowest one are flares
    #EXCEPT IF THEIR FLUENCE IS LESS THAN 8 COUNTS
    del_blocks = np.array([])
    if num_blocks > 1:    
        if (counts < 8).any():
            lowfluence = np.where(counts < 8)[0]
            #print('low counts:',counts[lowfluence])
            for i in lowfluence:
                if i == 0:
                    counts[i+1] = counts[i+1] + counts[i]
                    ledges[i+1] = ledges[i]
                    widths[i+1] = widths[i+1] + widths[i]
                    rates[i+1] = counts[i+1]/widths[i+1]
                    rateserr[i+1] = np.sqrt(counts[i+1])/widths[i+1]
                    del_blocks = np.append(del_blocks,i)
                else:
                    counts[i-1] = counts[i-1] + counts[i] 
                    redges[i-1] = redges[i]
                    widths[i-1] = widths[i-1] + widths[i]
                    rates[i-1] = counts[i-1]/widths[i-1]
                    rateserr[i-1] = np.sqrt(counts[i-1])/widths[i-1]
                    del_blocks = np.append(del_blocks,i)
                    
            counts = np.delete(counts,del_blocks.astype(int))
            ledges = np.delete(ledges,del_blocks.astype(int))
            redges = np.delete(redges,del_blocks.astype(int))
            widths = np.delete(widths,del_blocks.astype(int))
            rates = np.delete(rates,del_blocks.astype(int))
            rateserr = np.delete(rateserr,del_blocks.astype(int)) 
        #Note how many blocks there are in the obsid
        num_blocks = np.size(redges)
        quies_id = np.argmax(widths)#quiescent block is largest block
	#Identify each flaring block             
        flares_id = (rates-amplitude_criteria*rateserr)>(rates[quies_id] + amplitude_criteria*rateserr[quies_id])
        flares = ma.array(rates, mask = ~flares_id)
        data = np.ones(6).reshape(1,6)
        for h in range(num_blocks):
            data = np.concatenate((data, np.array([ledges[h], redges[h], counts[h], widths[h], rates[h], rateserr[h]]).reshape(1,6)))
            #print(ledges[h], redges[h], peakcr[h])
        data = np.delete(data, (0), axis = 0)
        LoRate = np.array([rates[quies_id], rateserr[quies_id]]).reshape(1,2)

        #UNPILE EACH FLARE BLOCK
        '''for p in range(num_blocks):
            if flares_id[p]:
                #print('block',p,'rate before pile-up corr ',rates[p], 'counts ', counts[p], 'rateserr ', rateserr[p])
                rates[p] = incident_cr[np.argmin(np.abs(rates[p] - observed_cr))]
                counts[p] = int(rates[p]*widths[p])
                rateserr[p] = np.sqrt(counts[p])/widths[p]
                #print('block',p,'rate after pile-up corr ',rates[p], 'counts ', counts[p], 'rateserr ', rateserr[p])
'''
        
        #If the obsid has only one flare
        if np.size(flares[~flares.mask]) == 1:
            indice = np.where(flares_id == True)[0][0]
            if l == 0:
                block = np.array([ledges[indice], redges[indice], counts[indice], widths[indice], 
                                                  rates[indice], rateserr[indice]]).reshape(1,6)
                l = l + 1
            else:
                block = np.concatenate((block, np.array([ledges[indice], redges[indice], counts[indice],
                                                  widths[indice], rates[indice], rateserr[indice]]).reshape(1,6)), axis = 0) 
        #If there are multiple blocks significantly above quiescence, then we need to figure out how many
        #flares there are
        else:
            k = 0        #Used to spot the first flare of the obsid
            j = 0        #Used to move through the blocks
            while j < num_blocks:
                if ~flares_id[j]:
                    j = j + 1
                    continue 
                else:
                    if k == 0:                  
                        k = k + 1
                        if (j < num_blocks - 1):
                            if ~flares_id[j + 1]:
                             #if the next block isnt a flare, then this flare is made of only one block
                                if l == 0: 
                                    block = np.array([ledges[j], redges[j], counts[j], widths[j], 
                                                  rates[j], rateserr[j]]).reshape(1,6)
                                    l = l + 1
                                    j = j + 1
                                else:
                                    block = np.concatenate((block, np.array([ledges[j], redges[j], counts[j], widths[j], 
                                                  rates[j], rateserr[j]]).reshape(1,6)), axis = 0)
                                    j = j + 1
                            else:
                           #But if the next block is also a flare then add the blocks until they end
                                flare_block = np.ones(6).reshape(1,6)
                                while(flares_id[j] and j < (num_blocks - 1)):
                                    flare_block = np.concatenate((flare_block, np.array([ledges[j], 
                                                     redges[j], counts[j], widths[j], 
                                                     rates[j], rateserr[j]]).reshape(1,6)), axis = 0)
                                    j = j + 1
                                    
                                if flares_id[j] and j == (num_blocks - 1):
                                    flare_block = np.concatenate((flare_block, np.array([ledges[j], 
                                                     redges[j], counts[j], widths[j], 
                                                     rates[j], rateserr[j]]).reshape(1,6)), axis = 0)
                                #Delete the first row that was used to initiate the array
                                flare_block = np.delete(flare_block, (0), axis = 0)
                                
                                #Finalize the block
                                if l == 0:
                                    l = l + 1
                                    block = np.array([flare_block[0,0],flare_block[-1,1],
                                                  np.sum(flare_block[:,2]), np.sum(flare_block[:,3]),
                                                  np.sum(flare_block[:,2])/float(np.sum(flare_block[:,3])),
                                                  math.sqrt(np.sum(flare_block[:,2]))/float(np.sum(flare_block[:,3]))]).reshape(1,6)
                                else:
                                    block = np.concatenate((block, np.array([flare_block[0,0],flare_block[-1,1],
                                                  np.sum(flare_block[:,2]), np.sum(flare_block[:,3]),
                                                  np.sum(flare_block[:,2])/float(np.sum(flare_block[:,3])), math.sqrt(np.sum(flare_block[:,2]))/
                                                     float(np.sum(flare_block[:,3]))]).reshape(1,6)), axis = 0)
                        else:
                            #If this is the last block, then this flare is also made up of only one block
                            if l == 0:
                                l = l + 1
                                block = np.array([ledges[j], redges[j], counts[j], widths[j], rates[j],
                                                        rateserr[j]]).reshape(1,6)
                            else:
                                block = np.concatenante((block,np.array([ledges[j], redges[j], counts[j], widths[j], rates[j],
                                                        rateserr[j]]).reshape(1,6)), axis = 0)
                            j = j + 1
                    
                    #If this isnt the first flare...
                    else:
                        if j < (num_blocks - 1):
                            #If this isnt the last block...
                            if ~flares_id[j + 1]:
                                #and if the next block isnt a flare, then this flare is made of only one block
                                block = np.concatenate((block, np.array([ledges[j], redges[j], counts[j], widths[j], 
                                                        rates[j], rateserr[j]]).reshape(1,6)), axis = 0)
                                j = j + 1
                            else:
                                #If the previous block wasnt a flare and this one  and the next are then 
                                #add the blocks until they end
                                flare_block = np.ones(6).reshape(1,6)
                                while(flares_id[j] and j < (num_blocks - 1)):
                                    flare_block = np.concatenate((flare_block, np.array([ledges[j], 
                                                     redges[j], counts[j], widths[j], 
                                                     rates[j], rateserr[j]]).reshape(1,6)), axis = 0)
                                    j = j + 1
                                    
                                if flares_id[j] and j == (num_blocks - 1):
                                    flare_block = np.concatenate((flare_block, np.array([ledges[j], 
                                                     redges[j], counts[j], widths[j], 
                                                     rates[j], rateserr[j]]).reshape(1,6)), axis = 0)
                                #Delete the first row that was use to initiate the array
                                flare_block = np.delete(flare_block, (0), axis = 0)
                                
                                #Finalize the block
                                block = np.concatenate((block,np.array([flare_block[0,0],flare_block[-1,1],
                                                  np.sum(flare_block[:,2]), np.sum(flare_block[:,3]),
                                                  np.sum(flare_block[:,2])/float(np.sum(flare_block[:,3])),math.sqrt(np.sum(flare_block[:,2]))/
                                                  float(np.sum(flare_block[:,3]))]).reshape(1,6)), axis = 0)
                        else:
                            if ~flares_id[j-1]:
                                #If this is the last block, then this flare is also made up of only one block
                                block = np.concatenate((block, np.array([ledges[j], redges[j], counts[j], widths[j], rates[j],
                                                        rateserr[j]]).reshape(1,6)), axis = 0)
                            j = j + 1
    else:
        data = np.array([ledges, redges, counts, widths, rates, rateserr]).reshape(1,6)
        LoRate = np.array([rates, rateserr]).reshape(1,2)
    
    #Make sure the arrays are sorted in time
    if block is not None:
        block = block[np.argsort(block[:,0])]
    data = data[np.argsort(data[:,0])]
    return data, block, LoRate

