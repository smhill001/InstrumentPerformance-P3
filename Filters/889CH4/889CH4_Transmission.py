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

SystemResponseFile="c:/Astronomy/Projects/Stars/Capella/Spectral Data/1D Spectra/Capella20201125031031_1D_WVCal.txt"
FilterResponseFile="c:/Astronomy/Projects/Stars/Capella/Spectral Data/1D Spectra/Capella20201125032159_1D_WVCal.txt"


###### Retrieve filter transmissions and convovle with disk integrated albedoes
SysResponse = np.fromfile(file=SystemResponseFile, dtype=float, count=-1, sep='\t')    
SysResponse=np.reshape(SysResponse,[int(SysResponse.size/2),2])

FilterResponse = np.fromfile(file=FilterResponseFile, dtype=float, count=-1, sep='\t')    
FilterResponse=np.reshape(FilterResponse,[int(FilterResponse.size/2),2])
Transmission889=GSU.SpectrumMath(FilterResponse,SysResponse,"Divide")

pl.subplot(1,1,1)
x0=840.
x1=940.
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
pl.plot(Transmission889[:,0],Transmission889[:,1],label='Filter Transmission',linewidth=1,color='b')

pl.legend(fontsize=7)

pl.savefig('c:/Astronomy/Projects/Techniques/InstrumentPerformance-P3/Filters/889CH4/889CH4_Transmission.png',dpi=320)
np.savetxt('c:/Astronomy/Projects/Techniques/InstrumentPerformance-P3/Filters/889CH4/889CH4_Transmission.txt',
           Transmission889,delimiter=" ",fmt="%10.3F %10.7F")
