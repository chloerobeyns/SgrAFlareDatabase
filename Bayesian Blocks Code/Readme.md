How to use the Bayesian Blocks Code: 

On terminal, go to the directory where all the python scripts are. To run the bayesian blocks analysis simply type a command of the form: 

>> python RUN.py 16218 False 

This assumes that there is a folder called "16218" in the same directory as RUN.py which contains two fits files: a level 2 events file named "16218_sgra_2-8keV_evt.fits" and a lightcurve file (here binned in 300s bins) named "16218_sgra_2-8keV_lc300.fits". The boolean command "False" indicates that the code will NOT be run on the magnetar region. 



If you also want to run the bayesian blocks analysis on the magnetar, run: 

>> python RUN.py 16218 True 

In this case, the folder "16218" should contain the files "16218_magnetar_2-8keV_evt.fits" and "16218_magnetar_2-8keV_lc300.fits". 



If you want to run the bayesian blocks code on multiple sets of data (for example 16218 and 2949) simply run the command: 

>> python RUN.py 16218 False 2949 False 

This can be done for any number of obsids given that they each have a folder containing the proper fits files. 

RETURN OF THE CODE: 
For each obsid specified, the code will create a "Results" directory within the folder "[obsid]/Results". In this "Results" directory there will be: a plot of the lightcurve and bayesian blocks, a file with the information about each block, and a file containing the data that is in the table (Flare duration, quiescent count rate, energy, luminosity etc.) 
