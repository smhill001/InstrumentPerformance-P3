# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 21:25:30 2022

@author: smhil
"""
import sys
sys.path.append('c:/Astronomy/Python Play')
sys.path.append('c:/Astronomy/Python Play/Util_P3')
sys.path.append('c:/Astronomy/Python Play/SPLibraries_P3')
sys.path.append('c:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
import matplotlib.pyplot as pl
import numpy as np
import GeneralSpecUtils_P3 as GSU
###### Retrieve filter transmissions and convovle with disk integrated albedoes

FilterFile="c:/Astronomy/Projects/Stars/Capella/Spectral Data/1D Spectra/Capella20201125024213_1D_WVCal.txt"
FilterOPNC = np.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
FilterOPNC=np.reshape(FilterOPNC,[int(FilterOPNC.size/2),2])

FilterFile="c:/Astronomy/Projects/Stars/Capella/Spectral Data/1D Spectra/Capella20201125025745_1D_WVCal.txt"
Filter672 = np.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
Filter672=np.reshape(Filter672,[int(Filter672.size/2),2])
Transmission672=GSU.SpectrumMath(Filter672,FilterOPNC,"Divide")

pl.subplot(1,1,1)
x0=620.
x1=720.
xtks=6
y0=0.0
y1=1.2
ytks=7

# Set x limits
pl.xlim(x0,x1)
# Set x ticks
pl.xticks(np.linspace(x0,x1,xtks, endpoint=True))
# Set y limits
pl.ylim(y0,y1)
pl.yticks(np.linspace(y0,y1,ytks, endpoint=True))
# Set y ticks
pl.grid(linewidth=0.2)
pl.tick_params(axis='both', which='major', labelsize=8)
pl.ylabel("Albedo x Transmission",fontsize=8,color="black")
pl.xlabel("Wavelength (nm)",fontsize=8)
pl.plot(Transmission672[:,0],Transmission672[:,1],label='Filter Transmission',linewidth=1,color='b')

pl.legend(fontsize=7)

pl.savefig('c:/Astronomy/Projects/Techniques/InstrumentPerformance-P3/Filters/672SII/672SII_Transmission.png',dpi=320)
np.savetxt('c:/Astronomy/Projects/Techniques/InstrumentPerformance-P3/Filters/672SII/672SII_Transmission.txt',
           Transmission672,delimiter=" ",fmt="%10.3F %10.7F")
