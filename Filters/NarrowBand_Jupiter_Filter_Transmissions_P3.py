# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 12:38:27 2021

This code creates ....

@author: Steven Hill
"""
import sys
sys.path.append('c:/Astronomy/Python Play')
sys.path.append('c:/Astronomy/Python Play/Util_P3')
sys.path.append('c:/Astronomy/Python Play/SPLibraries_P3')
sys.path.append('c:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
import matplotlib.pyplot as pl
import numpy as np
import GeneralSpecUtils_P3 as GSU

###### Get reference disk-integrated albedo from Karkoschka, 1994


###### Plot disk-integrated reference albedo and simulated ammonia-free albedo


#Plot Layout Configuration
# Set y ticks

###### Retrieve filter transmissions and convovle with disk integrated albedoes

FilterFile="c:/Astronomy/Projects/Stars/Capella/Spectral Data/1D Spectra/Capella20201204033249_1D_WVCal.txt"
FilterOPNA = np.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
FilterOPNA=np.reshape(FilterOPNA,[int(FilterOPNA.size/2),2])

FilterFile="c:/Astronomy/Projects/Stars/Capella/Spectral Data/1D Spectra/Capella20201204034206_1D_WVCal.txt"
Filter730 = np.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
Filter730=np.reshape(Filter730,[int(Filter730.size/2),2])
Transmission730=GSU.SpectrumMath(Filter730,FilterOPNA,"Divide")

##########

FilterFile="c:/Astronomy/Projects/Stars/Capella/Spectral Data/1D Spectra/Capella20201125024213_1D_WVCal.txt"
FilterOPNC = np.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
FilterOPNC=np.reshape(FilterOPNC,[int(FilterOPNC.size/2),2])

FilterFile="c:/Astronomy/Projects/Stars/Capella/Spectral Data/1D Spectra/Capella20201125025745_1D_WVCal.txt"
Filter672 = np.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
Filter672=np.reshape(Filter672,[int(Filter672.size/2),2])
Transmission672=GSU.SpectrumMath(Filter672,FilterOPNC,"Divide")

##########

FilterFile="C:/Astronomy/Projects/Stars/Capella/Spectral Data/1D Spectra/Capella20201125031031_1D_WVCal.txt"
FilterNIR = np.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
FilterNIR=np.reshape(FilterNIR,[int(FilterNIR.size/2),2])

FilterFile="C:/Astronomy/Projects/Stars/Capella/Spectral Data/1D Spectra/Capella20201125032159_1D_WVCal.txt"
Filter889 = np.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
Filter889=np.reshape(Filter889,[int(Filter889.size/2),2])
Transmission889=GSU.SpectrumMath(Filter889,FilterNIR,"Divide")

FilterFile="C:/Astronomy/Projects/Stars/Capella/Spectral Data/1D Spectra/Capella20201125034035_1D_WVCal.txt"
Filter940 = np.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
Filter940=np.reshape(Filter940,[int(Filter940.size/2),2])
Transmission940=GSU.SpectrumMath(Filter940,FilterNIR,"Divide")




##########

FilterFile="c:/Astronomy/Projects/Planets/Mars/Spectral Data/1D Spectra/Mars20201122014325_1D_WVCal.txt"
FilterOPNM = np.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
FilterOPNM=np.reshape(FilterOPNM,[int(FilterOPNM.size/2),2])

FilterFile="c:/Astronomy/Projects/Planets/Mars/Spectral Data/1D Spectra/Mars20201122020503_1D_WVCal.txt"
Filter658 = np.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
Filter658=np.reshape(Filter658,[int(Filter658.size/2),2])
Transmission658=GSU.SpectrumMath(Filter658,FilterOPNM,"Divide")

FilterFile="c:/Astronomy/Projects/Planets/Mars/Spectral Data/1D Spectra/Mars20201122015240_1D_WVCal.txt"
Filter656 = np.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
Filter656=np.reshape(Filter656,[int(Filter656.size/2),2])
Transmission656=GSU.SpectrumMath(Filter656,FilterOPNM,"Divide")

FilterFile="c:/Astronomy/Projects/Planets/Mars/Spectral Data/1D Spectra/Mars20201122020035_1D_WVCal.txt"
Filter647 = np.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
Filter647=np.reshape(Filter647,[int(Filter647.size/2),2])
Transmission647=GSU.SpectrumMath(Filter647,FilterOPNM,"Divide")

##########

FilterFile="c:/Astronomy/Projects/Stars/Vega/Spectral Data/1D Spectra/Vega20210727051700_1D_WVCal.txt"
FilterOPNV = np.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
FilterOPNV=np.reshape(FilterOPNV,[int(FilterOPNV.size/2),2])

FilterFile="c:/Astronomy/Projects/Stars/Vega/Spectral Data/1D Spectra/Vega20210727051317_1D_WVCal.txt"
Filter632 = np.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
Filter632=np.reshape(Filter632,[int(Filter632.size/2),2])
Transmission632=GSU.SpectrumMath(Filter632,FilterOPNC,"Divide")

##########


###### Plot filter transmissions convolved with disk-integrated albedos
pl.figure(figsize=(6.5, 4.0), dpi=150, facecolor="white")

x0=600.
x1=1000.
xtks=9
y0=0.0
y1=1.0
ytks=11

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
pl.plot(Transmission940[:,0],Transmission940[:,1],linewidth=1,color='C7',label='940NIR')
pl.plot(Transmission889[:,0],Transmission889[:,1],linewidth=1,color='C6',label='889CH4')
pl.plot(Transmission730[:,0],Transmission730[:,1],linewidth=1,color='C5',label='730OII')
pl.plot(Transmission672[:,0],Transmission672[:,1],linewidth=1,color='C0',label='672SII')
pl.plot(Transmission658[:,0],Transmission658[:,1],linewidth=1,color='C1',label='658NII')
pl.plot(Transmission656[:,0],Transmission656[:,1],linewidth=1,color='C2',label='656HIA')
pl.plot(Transmission647[:,0],Transmission647[:,1],linewidth=1,color='C3',label='647CNT')
pl.plot(Transmission632[:,0],Transmission632[:,1],linewidth=1,color='C4',label='632OI')

pl.legend(fontsize=7)


pl.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/NarrowBand_Jupiter_Filter_Transmissions.png',dpi=320)

path='c:/Astronomy/Projects/Techniques/InstrumentPerformance-P3/Filters/'

np.savetxt(path+'632OI/632OI_Transmission.txt',Transmission632,delimiter=" ",
           fmt="%10.3F %10.7F")
np.savetxt(path+'647CNT/647CNT_Transmission.txt',Transmission647,delimiter=" ",
           fmt="%10.3F %10.7F")
np.savetxt(path+'656HIA/656HIA_Transmission.txt',Transmission656,delimiter=" ",
           fmt="%10.3F %10.7F")
np.savetxt(path+'658NII/658NII_Transmission.txt',Transmission658,delimiter=" ",
           fmt="%10.3F %10.7F")
np.savetxt(path+'672SII/672SII_Transmission.txt',Transmission672,delimiter=" ",
           fmt="%10.3F %10.7F")
np.savetxt(path+'730OII/730OII_Transmission.txt',Transmission730,delimiter=" ",
           fmt="%10.3F %10.7F")
np.savetxt(path+'940NIR/940NIR_Transmission.txt',Transmission940,delimiter=" ",
           fmt="%10.3F %10.7F")
np.savetxt(path+'889CH4/889CH4_Transmission.txt',Transmission889,delimiter=" ",
           fmt="%10.3F %10.7F")

