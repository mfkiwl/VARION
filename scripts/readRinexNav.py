# \file readRinexNav.py
# \author Michela Ravanelli. 2020  michela.ravanelli(at)uniroma1.it
#
# -------------------------------------------------------------------------
#
# Copyright (C) 2015-2020  (see AUTHORS file for a list of contributors)
#
# VARION is a open source software for GNSS processing
#
# This file is part of VARION.
#
# VARION is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# VARION is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with VARION. If not, see <http://www.gnu.org/licenses/>.
#
# -------------------------------------------------------------------------
import numpy as np
import sys
import pdb
from datetime import datetime


def read_rinex_v2_new(rinexnav):
    
    f = open(rinexnav, "r")
    #pdb.set_trace()     

    lns = f.readlines() 
    count = 0
    # Skip header
    while "END OF HEADER" not in lns[count]:
				count += 1
    count += 1
    #diz = {}
    keys = ['sat','year','month','day','hour','minute','seconds','a0','a1','a2',
            'IODE','Crs','Deltan','M0',
            'Cuc','e','Cus','sqrtA',
            'TOE', 'Cic','OMEGA','Cis',
            'i0','Crc','omega0','OMEGA_DOT',
            'idot','Codes','GPSweek','L2Pdata',
            'SVa','SVh','TGD','IODC',
            'TTom','HFI']
            
    #diz = diz.fromkeys(keys, A)  
    diz = {k:[] for k in keys }     
    for i in np.arange(count,len(lns),8):
        
        ln = lns[i]
        ln = ln.replace("D","E")   
        #pdb.set_trace()

        diz['sat'] .append( float(ln[0:3]) )
        diz['year'].append( int(ln[3:5]) )
        diz['month'] .append(  int(ln[6:8]) )
        diz['day']  .append(  int(ln[9:11]) )
        diz['hour'] .append(  int(ln[12:14]) )
        diz['minute'] .append(  int(ln[15:17])   )          
        diz['seconds'] .append(  int(ln[18:20]) )
        diz['a0'] .append(  float(ln[22:41]) )
        diz['a1']  .append(  float(ln[41:60]) )
        diz['a2'] .append(  float(ln[60:80]) )
        
        ln = lns[i+1]
        ln = ln.replace("D","E")           
        diz['IODE'] .append(  float(ln[3:22]) )
        diz['Crs'] .append(  float(ln[22:41]) )
        diz['Deltan'] .append(  float(ln[41:60]) )
        diz['M0']  .append(  float(ln[60:80]) )
        
        ln = lns[i+2]
        ln = ln.replace("D","E")           
        diz['Cuc'] .append(  float(ln[3:22]) )
        diz['e'] .append(  float(ln[22:41]))
        diz['Cus'] .append(  float(ln[41:60]))
        diz['sqrtA']  .append(  float(ln[60:80]))
        
        ln = lns[i+3]
        ln = ln.replace("D","E")           
        diz['TOE'] .append(  float(ln[3:22]))
        diz['Cic'] .append(  float(ln[22:41]))
        diz['OMEGA'] .append(  float(ln[41:60]))
        diz['Cis']  .append(  float(ln[60:80]))
 
        ln = lns[i+4]
        ln = ln.replace("D","E")           
        diz['i0'] .append(  float(ln[3:22]))
        diz['Crc'] .append(  float(ln[22:41]))
        diz['omega0'] .append(  float(ln[41:60]))
        diz['OMEGA_DOT']  .append(  float(ln[60:80])   )    
        
        ln = lns[i+5]
        ln = ln.replace("D","E")           
        diz['idot'] .append(  float(ln[3:22]))
        diz['Codes'] .append(  float(ln[22:41]))
        diz['GPSweek'] .append(  float(ln[41:60]))
        diz['L2Pdata']  .append(  float(ln[60:80])  ) 

 
        ln = lns[i+6]
        ln = ln.replace("D","E")           
        diz['SVa'] .append(  float(ln[3:22]))
        diz['SVh'] .append(  float(ln[22:41]))
        diz['TGD'] .append(  float(ln[41:60]))
        diz['IODC']  .append(  float(ln[60:80]) ) 
        
        ln = lns[i+7]
        ln = ln.replace("D","E")           
        diz['TTom'] .append(  float(ln[3:22]))
        if (len(ln)>24):
                if (ln[24]) != " ":
                    diz['HFI'] .append(  float(ln[22:41]))
     
    return diz

###########################################################################


def nav_gps_v2_new(rinexnav):

        #pdb.set_trace()

        diz = read_rinex_v2_new(rinexnav)

        prn=[] ; sow=[];epoch=[]; SVclockBias=[];SVclockDrift=[];SVclockDriftRate=[];IODE=[];
        Crs=[];DeltaN=[];M0=[];Cuc=[];Eccentricity=[];Cus=[];a=[];TimeEph=[];
        Cic=[];OMEGA=[];CIS=[];Io=[]; Crc=[];omega=[];OMEGA_DOT=[];IDOT=[];
        CodesL2=[];GPSWeek=[];L2Pflag=[];SVacc=[];SVhealth=[];TGD=[];IODC=[];
        TransTime=[];FitIntvl=[];

        new_sat=[]
        #pdb.set_trace()
        for i in xrange(0,len(diz['sat'])):
          
            new_sat.append(diz['sat'][i])
            epoch.append(datetime(year =int(diz['year'][i]),
                                          month   =diz['month'][i],
                                          day     =diz['day'][i],
                                          hour    =diz['hour'][i],
                                          minute  =diz['minute'][i],
                                          second  =int(diz['seconds'][i])   )   )
            sow.append((diz['hour'][i])*60*60.+(diz['minute'][i])*60.+diz['seconds'][i])  
            prn.append( diz['sat'][i])             # PRN number
            TimeEph.append( diz['TOE'][i])             # Time of Ephemeris (seconds into GPS week) 
            a.append( (diz['sqrtA'][i])**2)            # Semi-major axis (meters)
            M0.append( diz['M0'][i])                   #Mean anomaly at reference time  (radians)
            Eccentricity.append( diz['e'][i])          # Eccentricity (unitless)
            Cuc.append( diz['Cuc'][i])                 # Amplitude of the cosine harmonic correction term to the argument of latitude (radians) 
            Cus.append( diz['Cus'][i])                 #Amplitude of the sine harmonic correction term to the argument of latitude (radians)
            Crc.append( diz['Crc'][i])                 #Amplitude of the cosine harmonic correction term to the orbit radius (meters)
            Crs.append( diz['Crs'][i])                 # Amplitude of the sine harmonic correction term to the orbit radius (meters)
            Cic.append( diz['Cic'][i])                 # Amplitude of the cosine harmonic correction term to the angle of inclination (radians)
            CIS.append( diz['Cis'][i])                 #Mean motion difference from computed value (radians/sec)
            DeltaN.append( diz['Deltan'][i])           #Mean motion difference from computed value (radians/sec)
            OMEGA.append( diz['OMEGA'][i])             # longitude of ascending node of orbit plane at weekly epoch
            omega.append( diz['omega0'][i])            # argument of perigee (radians)
            OMEGA_DOT.append( diz['OMEGA_DOT'][i])     # Rate of right ascension (radians/sec)
            Io.append( diz['i0'][i])                    # Inclination angle at reference time (radians)
            IDOT.append( diz['idot'][i])               #  Rate of inclination angle (radians/sec)
            GPSWeek.append( diz['GPSweek'][i])



        prn=(np.asarray(prn))
        sow=(np.asarray(sow))
        epoch=(np.asarray(epoch))
        te  = (np.asarray(TimeEph))      # Time of Ephemeris (seconds into GPS week) 
        a   = (np.asarray(a))               # Semi-major axis (meters)
        mo  = (np.asarray(M0))           # Mean anomaly at reference time  (radians)
        e   = (np.asarray(Eccentricity)) # Eccentricity (unitless)
        cuc = (np.asarray(Cuc))          # Amplitude of the cosine harmonic correction term to the argument of latitude (radians) 
        cus = (np.asarray(Cus))          # Amplitude of the sine harmonic correction term to the argument of latitude (radians)
        crc = (np.asarray(Crc))          # Amplitude of the cosine harmonic correction term to the orbit radius (meters)
        crs = (np.asarray(Crs))          # Amplitude of the sine harmonic correction term to the orbit radius (meters)
        cic = (np.asarray(Cic))          # Amplitude of the cosine harmonic correction term to the angle of inclination (radians)
        cis = (np.asarray(CIS))          # Amplitude of the sine harmonic correction term to the angle of inclination (radians)
        dn  = (np.asarray(DeltaN))       # Mean motion difference from computed value (radians/sec)
        lo = (np.asarray(OMEGA))         # longitude of ascending node of orbit plane at weekly epoch
        ome = (np.asarray(omega))        # argument of perigee (radians)
        ome_dot = (np.asarray(OMEGA_DOT))# Rate of right ascension (radians/sec)
        i_init  = (np.asarray(Io))       # Inclination angle at reference time (radians)
        i_rate  = (np.asarray(IDOT))     #  Rate of inclination angle (radians/sec)
        gps_week_sat = np.asarray(GPSWeek)


        return epoch, sow,prn, te,a,mo,e,cuc,cus,crc,crs,cic,cis,dn,lo,ome,ome_dot,i_init,i_rate,gps_week_sat

##################################################################################

