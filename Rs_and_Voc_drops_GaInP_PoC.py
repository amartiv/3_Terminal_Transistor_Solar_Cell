# -*- coding: utf-8 -*-
"""
Created on Mon Oct  12 12:41:40 2020

@author: Marius H. Zehender
https://orcid.org/0000-0002-2263-4560

Calculate series resistance of an Heterojuntion Bipolar Transistor Solar Cell (HBTSC)

method adapted from: 
    Luque, A. & Araújo, G. L. Solar cells and optics for photovoltaic concentration. (A. Hilger, 1989)
    
details on HBTSC: 
concept: Martí, A. & Luque, A. Three-terminal heterojunction bipolar transistor solar cell for high-efficiency photovoltaic conversion. Nature Communications 6, 6902 (2015)
design:  Antolín, E., Zehender, M. H., García-Linares, P., Svatek, S. A. & Martí, A. Considerations for the Design of a Heterojunction Bipolar Transistor Solar Cell. IEEE Journal of Photovoltaics 10, 2–7 (2020)
cell:    Zehender, M. H. et al. GaInP/GaAs three-terminal heterojunction bipolar transistor solar cell. in submission
"""
# constants
q = 16021766208e-29 #[C]
# ionized doping atoms densities
Ne = 4.2e17     #[cm³]
Nb = 6.3e18     #[cm³]
NbE= 6.4e18     #[cm³]
Nc = 4.7e17     #[cm³]
#material parameters
ne = Ne;
mue = 3500/(1+(ne/1e16)**(32e-2))      # emitter majority carrier mobility [cm²/V/s]   n-GaInP Parameters fit from:[DOI:10.1016/S0022-0248(01)01743-2],
# in accordance with [DOI:10.1063/1.112373] and [DOI:10.1063/1.343718]
mub = 40     #base majority carrier mobility [cm²/V/s]   #p-GaInP modified from 30:[DOI:10.1063/1.343718]
muc = 3000   #collector majority carrier mobility [cm²/V/s] #n-GaAs [DOI:10.1063/1.372274],[DOI: 10.1016/S0080-8784(08)60331-2]
#Semiconductor resistivities
rhoe    = 1 / (Ne * mue * q)    #emitter resistivity
rhob    = 1 / (Nb * mub * q)    #base resistivity
rhobE   = 1 / (NbE * mub * q)   #base resistivity
rhoc    = 1 / (Nc * muc * q)    #collector resistivity
#Semiconductor structure
we  = 485e-7                    #emitter thickness [cm]
wb   = 565e-7                   #base thickness [cm]
wbE = 373e-7                    #etched base thickness [cm]
wc  = 3000e-7;                  #collector thickness [cm]
#Semiconductor sheet resistances
Rse = rhoe  / we
Rsb = rhob  / wb
RsbE= rhobE / wbE
Rsc = rhoc  / wc
#Contact resistances
# n-GaAs:
ROce = 6e-5     #specific emitter contact resistance (emitter area) [Ohm*cm2]  #p-GaInP [DOI: 10.1002/1521-396X(200103)184:1<139::AID-PSSA139>3.0.CO;2-M]
ROcb = 3e-1     #specific base contact resistance (base area) [Ohm*cm2]
rcc  = 6e-5     #specific collector rear contact resistance [Ohm*cm2]
#CELL PARAMETERS
Lm = 1381e-4    #Mesa Length [cm]
Tm = 1422e-4    #Mesa Width [cm]
te = 15e-4      #emitter metal fingers width [cm]
le = 1114e-4    #emitter mal fingers length [cm]
de = 195e-4     #distance between emitter metal fingers [cm]
Te = Tm         #emitter width [cm]
Ler = 819e-4    #emitter length [cm]
Le = Ler / 2;   #two busbars, so we have to calculate only half emitter length
Rsme = 0.5      #emitter metal grid sheet resistance [sheetOhm]
letch = 100e-4  #distance between base bus and onset of emitter [cm]
Tb = Tm         #base width [cm]
Lb = (Lm - 200e-4) / 2    #base length [cm] # rest outer distance from base contact, two busbars
ROc = rhoc;     #collector resistivity [Ohm*cm]
ROm = 17e-6     #rear metal resistivity [Ohm*cm]
wmc = 600e-7    #rear metal thickness [cm]
Lc = Lm         #collector length [cm]
Tc = Tm         #collector width [cm]
at = 0.01026    #top junction area [cm^2]
ab = 0.01570    #bottom junction area [cm^2]
Jt = 14.2       #top junction current density at maximum power point [mA cm^-2]
Jb = 16.3       #top junction current density at maximum power point [mA cm^-2]


#EMITTER
Fse = 1-(te/de)                             #Emitter Transparency Factor
re  = (Rse * te**2) / (12 * ((1-Fse)**2))   #Specific Emitter Resistance (emitter area) in [Ohm*cm2]
rme = (Rsme * le**2) / (3 * (1 - Fse))      #Specific Emitter Metal Grid Resistance (emitter area) in [Ohm*cm2]
rce = ROce / (1 - Fse)                      #Specific Emitter Contact Resistance (emitter area) in [Ohm*cm2]
rE = re + rme + rce                         #Total Emitter Series Resistance in [Ohm*cm2]
print('\n')
print('Emitter series resistance = {:1.2e} '       .format(rE)+'\u03A9 cm^2')
print('Emitter resistance = {:1.2e} '              .format(re)+'\u03A9 cm^2')
print('Emitter contact resistance = {:1.2e} '      .format(rce)+'\u03A9 cm^2')
print('Emitter metal grid resistance = {:1.2e} '   .format(rme)+'\u03A9 cm^2')

#BASE
rBt =(1/3*Rsb*Le**2)+(RsbE*letch*Le)                #Base contribution to specific series resistance to the top jucntion current path
rBb =(1/3*(Le*Rsb+letch*RsbE)/(Le+letch)*Lb**2)     #Base contribution to specific series resistance to the bottom jucntion current path
print('\n')
print('Base series resistance (top junction) = {:1.2e} '   .format(rBt)+'\u03A9 cm^2')
print('Base series resistance (bottom junction) = {:1.2e} '.format(rBb)+'\u03A9 cm^2')

#COLLECTOR
rc = ROc*wc         # Specific collector semiconductor resistance
rmc= ROm*wmc        # Specific rear metal resistance
rC = rc+rmc+rcc     # Specific collector resistance
print('\n')
print('Collector series resistance = {:1.2e} '         .format(rC)+'\u03A9 cm^2')
print('Collector semiconductor resistance = {:1.2e} '  .format(rc)+'\u03A9 cm^2')
print('Collector contact resistance = {:1.2e} '        .format(rcc)+'\u03A9 cm^2')
print('Collector metal resistance = {:1.2e} '          .format(rmc)+'\u03A9 cm^2')

# whole jucntions
Rstop       = rE + rBt
Rsbottom    = rBb + rC
print('\n')
print('Top junction series resistance = {:1.2e} '    .format(Rstop)+'\u03A9 cm^2')
print('Bottom junction series resistance = {:1.2e} ' .format(Rsbottom)+'\u03A9 cm^2')

# voltage drops
Vt = rBb*Jb
Vb = rBt*Jt*at/ab
print('\n')
print('Voltage drop at top junction = {:1.2e} '    .format(Vt)+'mV')
print('Voltage drop at bottom junction = {:1.2e} ' .format(Vb)+'mV')