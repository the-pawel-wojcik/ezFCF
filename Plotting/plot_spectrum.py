#!/anaconda2/bin/python -tt

#Here I compute and plot XES spectra

import sys
import math

def main():

   #These needed for plotting
   import matplotlib.pyplot as plt # library to make plots
   from matplotlib.ticker import FormatStrFormatter
   import numpy as np # linear algebra library
   #These are general utility functions
   import spectral_tools as spt

   #To execute, run
   #python plot_spectrum.py trans.spectrum
   #It will produce a convoluted spectrum (in trans.spectrum.dat file)
   #and a pdf with the plot in trans.spectrum.pdf

   fname1=sys.argv[1]
   labels1=[]
   data1 = spt.read_data_with_labels(fname1,labels1)

   #now generate convoluted spectrum 
   x_0=data1[0,0]
   npts=len(data1[:,0])
   print "npts=", npts
   x_max=data1[npts-1,0]
   step=abs(x_max-x_0)/npts
   print "X-range:", x_0, x_max, npts, step  

   #gaussian width, in eV 
   #FWHM
   #width=0.025
   width=0.025
   spectrum1=spt.compute_spectrum(abs(x_0),abs(x_max),step,abs(data1[:,0]),data1[:,1],width)
   #print "Spectrum 1", spectrum1
   
   ticks_step=0.1 
   xtics=spt.compute_ticks(abs(x_0),abs(x_max),ticks_step)

   width = 8
   rect = (0.2, 0.2, 0.75, 0.75)

   fig = plt.figure(1,(width,width))
   ax = fig.add_axes(rect)
   # 2.0 for HCOH, 80 for thymine and 15 for adenine
   #ylim=2.0 
   ylim=40.0
   ax.set_xlim(abs(x_0),abs(x_max))
   ax.set_ylim(0.0,  ylim)  
   ax.set_xticks(xtics)
   ax.plot(spectrum1[:,0], spectrum1[:,1],'-', color='blue', label='FCFs')

   plt.ylabel('FCFs',fontsize=12)
   plt.xlabel('Energy, eV', fontsize=12)
   plt.legend()

   fname2=fname1+'.pdf'
   plt.savefig(fname2)

   

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
   
