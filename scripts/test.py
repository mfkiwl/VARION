# -*- coding: utf-8 -*-
import numpy as np

import readRinexNav as RN
import myRead as mR
import myObs as mO

##CONSTANTS
mu = 3.986005e14
omegae = 7.2921151467e-5
pigreco = 3.1415926535898

## LETTURA FILE NAVIGAZIONALE
data = RN.readRinexNav( 'ahup3020.12n' )

prn = (np.asarray(data['sv']))           # PRN number
te  = (np.asarray(data['TimeEph']))      # Time of Ephemeris (seconds into GPS week) 
a   = (np.asarray(data['sqrtA'])) **2    # Semi-major axis (meters)
mo  = (np.asarray(data['M0']))           # Mean anomaly at reference time  (radians)
e   = (np.asarray(data['Eccentricity'])) # Eccentricity (unitless)
cuc = (np.asarray(data['Cuc']))          # Amplitude of the cosine harmonic correction term to the argument of latitude (radians) 
cus = (np.asarray(data['Cus']))          # Amplitude of the sine harmonic correction term to the argument of latitude (radians)
crc = (np.asarray(data['Crc']))          # Amplitude of the cosine harmonic correction term to the orbit radius (meters)
crs = (np.asarray(data['Crs']))          # Amplitude of the sine harmonic correction term to the orbit radius (meters)
cic = (np.asarray(data['Cic']))          # Amplitude of the cosine harmonic correction term to the angle of inclination (radians)
cis = (np.asarray(data['CIS']))          # Amplitude of the sine harmonic correction term to the angle of inclination (radians)
dn  = (np.asarray(data['DeltaN']))       # Mean motion difference from computed value (radians/sec)
lo = (np.asarray(data['OMEGA']))         # longitude of ascending node of orbit plane at weekly epoch
ome = (np.asarray(data['omega']))        # argument of perigee (radians)
ome_dot = (np.asarray(data['OMEGA DOT']))# Rate of right ascension (radians/sec)
i_init  = (np.asarray(data['Io']))       # Inclination angle at reference time (radians)
i_rate  = (np.asarray(data['IDOT']))     #  Rate of inclination angle (radians/sec)
no = (mu/(a**3))**0.5
n  = no + dn

data_obs   = mR.read_rinex( 'ahup3020.12o' )

prn_sat_list = []
Xk_list = []
Yk_list = []
Zk_list = []
sod_list = []
toe_list = []
for s in xrange(1,33):
    prn_sat = s
    if prn_sat < 10:
        varion, ora, sod = mO.obs_sat( data_obs[0], data_obs[1], data_obs[2], data_obs[3], data_obs[4], 'G0'+str(prn_sat))
    else:
        varion, ora, sod = mO.obs_sat( data_obs[0], data_obs[1], data_obs[2], data_obs[3], data_obs[4], 'G'+str(prn_sat))
    for i in xrange(len(sod)):
        tk_arr = ( sod[i] - te[ prn==prn_sat ] )                                           # tk = t - te
        val, idx = min((val, idx) for (idx, val) in enumerate( np.abs(tk_arr) ))   # vedo quale tk è min
        tk = tk_arr[idx]
        mk = mo[ prn == prn_sat][ idx] + n[ prn == prn_sat][ idx]*tk        
        ## iterazioni per calcolare Ek
        eko = mk
        for j in xrange(50):   
            ek = mk + e[ prn == prn_sat][idx]*np.sin(eko)
            if abs( ek-eko ) <= 10e-10: break
            else:
                eko = ek
            if j == 49:
                print "Iterations didn't go through in 30 steps"
                break
        ##
        vk = np.arctan2(   ( (( 1- ( e[ prn == prn_sat][idx] )**2 )**0.5)*np.sin(ek) ) , ( np.cos(ek)-e[ prn == prn_sat][idx] )  )
        uk = ome[ prn == prn_sat ][idx] + vk
        duk = (  cuc[ prn == prn_sat][idx]*np.cos(2*uk) + (cus[ prn == prn_sat][idx] * np.sin(2*uk)) )
        drk = (  crc[ prn == prn_sat][idx]*np.cos(2*uk) + (crs[ prn == prn_sat][idx] * np.sin(2*uk)) )
        dik = (  cic[ prn == prn_sat][idx]*np.cos(2*uk) + (cis[ prn == prn_sat][idx] * np.sin(2*uk)) )
        omek = ome[ prn == prn_sat ][idx] + duk
        rk   = a[ prn == prn_sat ][idx]*( 1- e[ prn == prn_sat ][idx]*np.cos(ek)) + drk
        ik   = i_init[ prn == prn_sat ][idx]  +  i_rate[ prn == prn_sat ][idx] * tk + dik
        xk   = rk*np.cos(omek+vk)
        yk   = rk*np.sin(omek+vk)
        lk   = lo[ prn == prn_sat ][idx]  + ome_dot[ prn == prn_sat ][idx] * tk - omegae * (sod[i])               
        Xk   = xk*np.cos(lk) - yk * np.sin(lk) * np.cos(ik)
        Yk   = xk*np.sin(lk) + yk * np.cos(lk) * np.cos(ik)
        Zk   = yk*np.sin(ik)
        Xk_list.append( Xk )
        Yk_list.append( Yk )
        Zk_list.append( Zk )
        prn_sat_list.append( prn_sat )
        sod_list.append( sod[i])
        toe_list.append( te[ prn==prn_sat ][ idx] )

prn_sat_arr = np.asarray( prn_sat_list )
sod_arr = np.asarray( sod_list )
toe_arr = np.asarray( toe_list )
Xk_arr = np.asarray( Xk_list )
Yk_arr = np.asarray( Yk_list )
Zk_arr = np.asarray( Zk_list )

with open('ahup3020.12o_VARION.sky', 'w') as f:
    f.write( 'PRN' + '\t' + 'sod' + '\t' + 'Xs' + '\t' + 'Ys' + '\t' + 'Zs'+ '\t' + 'toe' + '\n' )
    for s in xrange(1,33):
        prn_sat = s
        for i in xrange(len(prn_sat_arr[ prn_sat_arr == s ])):
            f.write( str(prn_sat_arr[ prn_sat_arr == s ][i])+'\t'+str(sod_arr[ prn_sat_arr == s ][i])+ '\t'+\
                    str(Xk_arr[ prn_sat_arr == s ][i])+'\t'+str(Yk_arr[ prn_sat_arr == s ][i])+'\t'+str(Zk_arr[ prn_sat_arr == s ][i])+ \
                     '\t'+str(toe_arr[ prn_sat_arr == s ][i])+'\n' )