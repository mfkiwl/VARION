# -*- coding: utf-8 -*-
import numpy as np

import readRinexNav as RN
import myRead as mR
import myObs as mO
import myFunc as mF


def coord_satellite( rinex_nav, rinex_obs ):
    '''
        inputs:
            - rinex navigation file (anche brdc)
            - rinex obs file
        outputs:
            - array of the prn
            - array of the sod
            - array of the X of the satellites
            - array of the Y of the satellites
            - array of the Z of the satellites
            - array of the time of Time of Ephemeris (seconds into GPS week)
    '''
    ##CONSTANTS
    mu = 3.986005e14
    omegae = 7.2921151467e-5
    pigreco = 3.1415926535898

    ## READING THE NAVIGATION RINEX 
    data = RN.readRinexNav( rinex_nav )

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
    ## READING THE NAVIGATION RINEX -- > in order to compute the Xs,Ys,Zs for the epochs we want
    data_obs   = mR.read_rinex( rinex_obs )

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
                    print "Iterations didn't go through in 50 steps"
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

    return prn_sat_arr, sod_arr, Xk_arr, Yk_arr, Zk_arr, toe_arr
#####
def track_sat(sIP, sat_name):
    '''
    Function that returns the geodetic coordinats of a particular satellites.
    input:
        - sIP list
        - name of the satellite (e.g. 'G13')
    output:
        - sod
        - names
        - X
        - Y
        - Z
        - toe
    '''    
    mask = ( sIP[0] == sat_name )
    sod   =  sIP[1][mask]
    names =  sIP[0][mask]
    x     =  sIP[2][mask] 
    y     =  sIP[3][mask] 
    z     =  sIP[4][mask] 
    toe   =  sIP[5][mask] 
        
    return sod, names, x, y, z, toe
#####
def coord_ipps( xr,yr,zr,xs_arr,ys_arr,zs_arr,h_iono):
    X1, Y1, Z1 = xr, yr, zr

    phi_list   = []
    lamda_list = []
    h_list     = []

    for i in xrange(len(xs_arr)):
        X2=xs_arr[i]
        Y2=ys_arr[i]
        Z2=zs_arr[i]

        while True:
            X3=(X1+X2)/2
            Y3=(Y1+Y2)/2
            Z3=(Z1+Z2)/2
            phi3,lamda3,h3=mF.coord_geog(X3,Y3,Z3)
            if abs(h3-h_iono)<1000:
                phi_list.append(phi3)
                lamda_list.append(lamda3)
                h_list.append(h3)
                break
            if h3>h_iono:
                X2=X3 
                Y2=Y3
                Z2=Z3
            if h3<h_iono:
                X1=X3 
                Y1=Y3
                Z1=Z3
    return np.asarray(phi_list), np.asarray(lamda_list), np.asarray(h_list)        


# with open('ahup3020.12n_VARION.sky', 'w') as f:
#     f.write( 'PRN' + '\t' + 'sod' + '\t' + 'Xs' + '\t' + 'Ys' + '\t' + 'Zs'+ '\t' + 'toe' + '\n' )
#     for s in xrange(1,33):
#         prn_sat = s
#         for i in xrange(len(prn_sat_arr[ prn_sat_arr == s ])):
#             f.write( str(prn_sat_arr[ prn_sat_arr == s ][i])+'\t'+str(sod_arr[ prn_sat_arr == s ][i])+ '\t'+\
#                     str(Xk_arr[ prn_sat_arr == s ][i])+'\t'+str(Yk_arr[ prn_sat_arr == s ][i])+'\t'+str(Zk_arr[ prn_sat_arr == s ][i])+ \
#                      '\t'+str(toe_arr[ prn_sat_arr == s ][i])+'\n' )