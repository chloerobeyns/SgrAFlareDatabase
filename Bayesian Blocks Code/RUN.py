#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  5 16:02:49 2020

@author: Elisa
"""

#This script takes as input an ObsID number and runs the 
#flare detection algorithm on it. 

#IMPORTANT: 
#this code must be in the same directory as CompleteBayesianBlocks 
#and the observation files must all be in sub-folders within this 
#directory

import CompleteBB as bb
import matplotlib.pyplot as plt
import os
import sys

    
''' This script runs the flare detection on sgrA data (and magnetar, if True)
    
Expected notation: 
        
- all the required data files are within a directory ./obsid 
- The files names are: 
    obsid_sgra_2-8keV_evt.fits (for event files)
    obsid_sgrA_2-8keV_lc300.fits (for lightcurves)


------------------------------------------------------------------

ARGUMENTS : 

    obsid: int 
            a number that denotes the observation id for which we want 
            to write the code
    magnetar: bool 
        says whether we want to also run the analysis on the 
        magnetar region
    
RETURNS: 
    creates a folder in  directory ./obsid/results 
    with the results of the flare detection. 
    
    Also saves figure of lightcurve with blocks in the same directory
    '''


#read arguments from command line: 
#NOTE: this allows us to run the script for multiple observations in
#      one go. 
    
#for more than 1 observation the order of arguments : 
#       Obsid1, Magnetar1, ObsID2, Magnetar2 (...) 
    
    

maxi = len(sys.argv)
i = 0
while i < maxi-1: 
    obsid = sys.argv[i+1]
    magnetar = sys.argv[i+2]
    
    #create string for datafile directories: (input)
    lc = "./" +  str(obsid) + "/" + (str(obsid) + "_sgra_2-8keV_lc300.fits")
    evt = "./" + str(obsid) + "/" + (str(obsid) +  "_sgra_2-8keV_evt.fits")
    
    
    
    #OUTPUTs directories:
    
    filename = "./"  + str(obsid) + "/" + "Results/" 
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    
    bb_info = "./"  + str(obsid) + "/" + "Results/"  + str(obsid) + "_sgra_bayesianBlocks_info.txt" #block info 
    plot = "./" + str(obsid) + "/" + "Results/" + str(obsid) + "_PLOT_sgra.png" #plot 
    table_res = "./" + str(obsid) + "/" + "Results/"  + str(obsid) + "_SGRA_TABLE_RESULTS.txt" #info for flare table
    
    print("running code for ObsID " + str(obsid))
    #run bayesian block: 
    bb.process(evt, bb_info)
    
    #Create the plot: 
    fig = plt.figure()
    bb.plot_bb(bb_info) 
    bb.plot_lc(lc) 
    plt.xlabel("Time (days)")
    plt.ylabel("Count Rate")
    plt.title("ObsID " + str(obsid) + " - SgrA")
    fig.savefig(plot)
    
    #Get flare information for database: 
    bb.getInfo(evt , lc , bb_info, table_res)
    
    #run this on the magnetar, if magnetar = true: 
    
    if magnetar == True :
        #create string for datafile directories: (input)
        lc_m = "./"  + str(obsid) + "/" + (str(obsid) + "_magnetar_2-8keV_lc300.fits")
        evt_m = "./" + str(obsid) + "/" + (str(obsid) +  "_magnetar_2-8keV_evt.fits")
    
        #OUTPUTs directories:
        bb_info_m = "./" + str(obsid) + "/" + "Results/"  + str(obsid) + "_magnetar_bayesianBlocks_info.txt" #block info 
        plot_m = "./"  + str(obsid) + "/" + "Results/" + str(obsid) + "_PLOT_magnetar.png" #plot 
        table_res_m = "./" + str(obsid) + "/" + "Results/"  + str(obsid) + "_magnetar_TABLE_RESULTS.txt" #info for flare table
    
        #run bayesian block: 
        bb.process(evt_m, bb_info_m)
    
        #Create the plot: 
        fig = plt.figure()
        bb.plot_bb(bb_info_m) 
        bb.plot_lc(lc_m) 
        plt.xlabel("Time (days)")
        plt.ylabel("Count Rate")
        plt.title("ObsID " + str(obsid) + " - magnetar")
        fig.savefig(plot_m)
    
        #Get flare information for database: 
        bb.getInfo(evt_m , lc_m , bb_info_m, table_res_m)
        
    i = i + 2 #update value of i 
 
    

