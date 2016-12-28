# -*- coding: utf-8 -*-
# \author Giorgio Savastano, 2015. giorgio.savastano(at)uniroma1.it
#
# -------------------------------------------------------------------------
#
# Copyright (C) 2015-2016  (see AUTHORS file for a list of contributors)
#
# VARION is a opean source software for GNSS processing
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
import readRinexNav as RN
import myObs as mO
import myFunc as mF
##
def sat_selection( rinex_obj, sats_list, start, stop ):
    """
    Functions that select the satellites in view in a time period
    """
    sats_list_new = []
    for sa in sats_list:
        mask_sate = (rinex_obj.data[0] == sa)
        mask_time = (rinex_obj.data[2][mask_sate] >= start) & (rinex_obj.data[2][mask_sate] <= stop)
        length = (stop - start) / (rinex_obj.int)
        if len(rinex_obj.data[2][mask_sate][mask_time]) <= ( length / 4.0 ):
            continue
        else:
            sats_list_new.append(sa)
    return np.asarray( sats_list_new )
##
def coord_satellite( rinex_nav, rinex_obs, sats_write ):
    '''
        inputs:
            - rinex navigation file (also brdc)
            - rinex obs file
        outputs:
            - array of the prn
            - array of the sod
            - array of the X of the satellites
            - array of the Y of the satellites
            - array of the Z of the satellites
            - array of the time of Time of Ephemeris (seconds into GPS week)
    '''
    ##### CONSTANTS

    # WGS 84 value of earth's univ. grav. par.
    mu = 3.986005e14
    # WGS 84 value of earth's rotation rate
    omegae = 7.2921151467e-5

    pigreco = 3.1415926535898
    speedOfLight = 299792458.0

    # relativistic correction term constant (not used so far)
    F = -4.442807633E-10

    # WGS-84 earth rotation rate
    we = 7.292115e-5

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
    gps_week_sat = np.asarray(data['GPSWeek'])

    no = (mu/(a**3))**0.5
    n  = no + dn
    ## READING THE OBSERVATION RINEX -- > in order to compute the Xs,Ys,Zs at the epochs we want
    data_obs   = rinex_obs.data

    prn_sat_list = []
    Xk_list = []; Yk_list = []; Zk_list = []
    sod_list = []
    toe_list = []
    for s in sats_write:
        prn_sat = int( s[1:3] )
        varion, ora, sod = mO.obs_sat( data_obs[0], data_obs[1], data_obs[2], data_obs[3], data_obs[4], s )

        for i in xrange( len(sod) ):
            
            # Time of flight
            tof = data_obs[-1][1:][i] / speedOfLight

            ### correction due to the earth rotation
            # alpha = tof * we

            trasmitTime = (sod[i] + rinex_obs.gps_sow_ref) - tof
            tk_arr = ( trasmitTime - te[ prn==prn_sat ] )                              # tk = t - te
            val, idx = min((val, idx) for (idx, val) in enumerate( np.abs(tk_arr) ))   # vedo quale tk è min
            tk = tk_arr[idx] + (rinex_obs.gps_week_ref - gps_week_sat[ prn==prn_sat ][idx])*604800

            ## TEST
            if tk > 302400:
                tk = tk - 604800
            if tk < -302400:
                tk = tk + 604800

            mk = mo[ prn == prn_sat][ idx] + n[ prn == prn_sat][ idx]*tk        
            ## iterazioni per calcolare Ek
            eko = mk
            for j in xrange(50):   
                ek = mk + e[ prn == prn_sat][idx]*np.sin(eko)
                if abs( ek-eko ) <= 10e-10: break
                else:
                    eko = ek
                if j == 49:
                    print "WARNING: Kepler Eqn didn't converge for sat " + str(s)
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

            lk   = lo[ prn == prn_sat ][idx]  + ome_dot[ prn == prn_sat ][idx] * tk - omegae * (sod[i]+rinex_obs.gps_sow_ref)  

            Xk   = xk*np.cos(lk) - yk * np.sin(lk) * np.cos(ik)
            Yk   = xk*np.sin(lk) + yk * np.cos(lk) * np.cos(ik)
            Zk   = yk*np.sin(ik)

            # Correction the satellite position for the time it took the message to get to the reciver
            #Xk   =  Xk * np.cos(alpha) + Yk * np.sin(alpha)
            #Yk   = -Xk * np.sin(alpha) + Yk * np.cos(alpha)

            Xk_list.append( Xk )
            Yk_list.append( Yk )
            Zk_list.append( Zk )
            prn_sat_list.append( s )
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
def track_sat(sIP, sat_name, start, stop ):
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
    mask_st = ( sIP[0] == sat_name )
    mask_tm = ( sIP[1][mask_st] >= start ) & ( sIP[1][mask_st] <= stop )

    sod   =  sIP[1][mask_st][mask_tm]
    names =  sIP[0][mask_st][mask_tm]
    x     =  sIP[2][mask_st][mask_tm] 
    y     =  sIP[3][mask_st][mask_tm] 
    z     =  sIP[4][mask_st][mask_tm] 
    toe   =  sIP[5][mask_st][mask_tm] 
        
    return sod, names, x, y, z, toe
#####
def coord_ipps( xr,yr,zr,xs_arr,ys_arr,zs_arr,h_iono):
    '''
        Function that computes the phi,lambda and h of the IPP. It needs as input the Xr,Yr,Zr of the receiver and the
        Xs,Ys,Zs whose are the array with the position of the satellite (it is moving).
        input:
            - Xr, Yr, Zr --> 1 position
            - Xs, Ys, Zs --> arrays of position
        output:
            - PHIipp, LAMBDAipp, Hipp --> arrays of the same size of the one in input
    '''

    phi_list   = []
    lamda_list = []
    h_list     = []

    for i in xrange(len(xs_arr)):
        X1, Y1, Z1 = xr, yr, zr
        
        X2=xs_arr[i]
        Y2=ys_arr[i]
        Z2=zs_arr[i]

        while True:
            X3=(X1+X2)/2
            Y3=(Y1+Y2)/2
            Z3=(Z1+Z2)/2
            phi3,lamda3,h3=mF.coord_geog(X3,Y3,Z3)
            if abs(h3-h_iono)<(50):
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
